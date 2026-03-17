import sympy as sp
import networkx as nx
import heapq


class IdPool:
    '''
    编号池: 分配和回收整数编号

    分配时自动选择最小的可用编号 (从 0 开始),
    回收后编号可被后续分配重新使用。
    '''
    def __init__(self):
        self._in_use = set()       # 当前已分配的编号
        self._recycled = []        # 最小堆, 存放已回收的编号
        self._next_new = 0         # 下一个全新编号 (未被分配过的最小值)

    def allocate(self) -> int:
        '''分配最小的可用编号'''
        # 优先从回收池中取最小编号
        while self._recycled:
            candidate = heapq.heappop(self._recycled)
            if candidate not in self._in_use:
                self._in_use.add(candidate)
                return candidate
        # 回收池为空, 分配全新编号
        new_id = self._next_new
        self._next_new += 1
        self._in_use.add(new_id)
        return new_id

    def release(self, id_val: int):
        '''回收编号, 使其可被再次分配'''
        if id_val not in self._in_use:
            raise KeyError(f"编号 {id_val} 未被分配, 无法回收")
        self._in_use.discard(id_val)
        heapq.heappush(self._recycled, id_val)

    @property
    def used(self) -> set:
        '''当前已分配的编号集合'''
        return set(self._in_use)


class Component:
    '''定义元件, 用于校验输入'''
    def __init__(self, u: int, v: int, tp: str, value):
        self.u = min(u, v)
        self.v = max(u, v)
        self.type = tp.upper()
        if self.type not in {'C', 'L', 'JJ'}:
            raise ValueError(f"元件类型 '{self.type}' 不在支持列表 {{'C', 'L', 'JJ'}} 中")
        self.value = sp.symbols(value) if isinstance(value, str) else value

    @property
    def weight(self):
        '''生成树权重: JJ=0 < L=1 < C=2, 最小权重优先进入树枝 (电容倾向连支)'''
        return 2 if self.type == 'C' else 1 if self.type == 'L' else 0


class MutualInductance:
    '''定义互感系数'''
    def __init__(self, comp_id1: int, comp_id2: int, value):
        self.comp_id1 = comp_id1
        self.comp_id2 = comp_id2
        self.value = sp.symbols(value) if isinstance(value, str) else value


class PhysicalFlux:
    '''
    定义物理外磁通 (与物理回路绑定)

    物理回路由一组元件 ID 按顺序围成, 顺序决定磁通正方向 (右手定则)
    '''
    def __init__(self, comp_loop: list[int], flux):
        self.comp_loop = list(comp_loop)
        self.flux = sp.symbols(flux) if isinstance(flux, str) else flux


class Circuit:
    '''
    超导电路量子化对象 

    使用流程:
        1. 创建空电路:      circuit = Circuit()
        2. 添加元件:        cid = circuit.add_component(u, v, type, value)
        3. 添加互感 (可选): mid = circuit.add_mutual(cid1, cid2, value)
        4. 设置外磁通 (二选一):
           方式 A — 物理回路磁通 (需提供恰好 num_loops 个线性无关回路):
                            fid = circuit.add_physical_flux(comp_loop, flux)
                            磁通为 0 的回路也需添加: add_physical_flux(..., 0)
           方式 B — 直接指定基本回路磁通 (可指定任意个):
                            circuit.set_external_flux(loop_index, flux)
                            未指定的回路默认磁通为 0
        5. 查看信息:        circuit.print_edges() / circuit.print_loops()
        6. 计算哈密顿量:    H, info = circuit.hamiltonian()

    元件/互感/磁通可随时增删, 所有计算按需自动重建图结构
    两种磁通设置方式不可混用, 同时存在时会报错。

    外磁通计算原理 (物理回路模式):
        每个物理回路可分解为基本回路的线性组合, 系数由回路在各连支上的分量决定。
        将所有物理回路排成矩阵方程 A · φ_ext = b, 解出各基本回路的外磁通。
        要求物理回路数 = 基本回路数 (连支数), 且矩阵 A 可逆 (回路线性无关)。
    '''

    def __init__(self):
        self._components = {}        # comp_id -> Component
        self._mutuals = {}           # mutual_id -> MutualInductance
        self._physical_fluxes = {}   # flux_id -> PhysicalFlux
        self._fundamental_fluxes = {}  # chord_key -> flux (直接指定基本回路磁通)
        self._comp_pool = IdPool()     # 元件编号池
        self._mutual_pool = IdPool()   # 互感编号池
        self._flux_pool = IdPool()     # 物理磁通编号池
        self._node_pool = IdPool()     # 节点编号池
        self._node_refcount = {}       # 节点引用计数: node_id -> 引用该节点的元件数
        self._dirty = True
        self._clear_cache()

    def _clear_cache(self):
        '''清除所有缓存的计算结果'''
        self._G = None
        self._T = None
        self._edge_map = None
        self._loops = None
        self._nt = None
        self._m = None
        self._mutual_dict = None
        self._F_C_cache = None
        self._comp_to_edge = None   # comp_id -> edge_key (重编号后)

    def _mark_dirty(self):
        '''标记图结构已变更, 需要重新构建'''
        self._dirty = True
        self._clear_cache()
        # 图结构变更后 chord_key 可能改变, 清除直接指定的基本回路磁通
        if self._fundamental_fluxes:
            import warnings
            warnings.warn(
                "图结构已变更，之前通过 set_external_flux 设置的基本回路磁通已被清除。"
                "请在图结构稳定后重新设置。",
                UserWarning, stacklevel=3
            )
            self._fundamental_fluxes.clear()

    def _ensure_built(self):
        '''确保图结构已构建, 按需重建'''
        if self._dirty:
            if not self._components:
                raise ValueError("电路中没有元件，无法构建")
            self._rebuild()
            self._dirty = False

    # 节点管理

    def _ensure_node(self, node_id: int):
        '''确保节点已注册, 未注册则自动注册 (兼容手动指定节点编号的旧用法)'''
        if node_id not in self._node_refcount:
            # 手动指定的节点: 在编号池中标记为已使用
            if node_id not in self._node_pool._in_use:
                # 将跳过的编号加入回收池, 使 add_node() 能分配到它们
                for skipped in range(self._node_pool._next_new, node_id):
                    if skipped not in self._node_pool._in_use:
                        heapq.heappush(self._node_pool._recycled, skipped)
                self._node_pool._in_use.add(node_id)
                if node_id >= self._node_pool._next_new:
                    self._node_pool._next_new = node_id + 1
            self._node_refcount[node_id] = 0

    def _ref_node(self, node_id: int):
        '''增加节点引用计数'''
        self._ensure_node(node_id)
        self._node_refcount[node_id] += 1

    def _unref_node(self, node_id: int):
        '''减少节点引用计数, 归零时自动回收'''
        if node_id in self._node_refcount:
            self._node_refcount[node_id] -= 1
            if self._node_refcount[node_id] <= 0:
                del self._node_refcount[node_id]
                if node_id in self._node_pool._in_use:
                    self._node_pool.release(node_id)

    @property
    def nodes(self) -> set:
        '''当前活跃的节点编号集合'''
        return set(self._node_refcount.keys())

    # 元件管理

    def add_component(self, u: int, v: int, tp: str, value) -> int:
        '''
        添加元件

        参数:
            u, v: 节点编号 (可使用 add_node() 自动分配, 或手动指定)
            tp: 元件类型 ('C', 'L', 'JJ')
            value: 参数值 (字符串自动转为 sympy 符号)

        返回: 元件 ID (用于后续引用)
        '''
        comp_id = self._comp_pool.allocate()
        self._components[comp_id] = Component(u, v, tp, value)
        # 更新节点引用计数
        comp = self._components[comp_id]
        self._ref_node(comp.u)
        self._ref_node(comp.v)
        self._mark_dirty()
        return comp_id

    def remove_component(self, comp_id: int):
        '''删除元件 (关联的互感会一并删除, 孤立节点自动回收)'''
        if comp_id not in self._components:
            raise KeyError(f"元件 {comp_id} 不存在")
        comp = self._components[comp_id]
        del self._components[comp_id]
        self._comp_pool.release(comp_id)
        # 减少节点引用计数 (可能触发节点回收)
        self._unref_node(comp.u)
        self._unref_node(comp.v)
        # 移除引用了该元件的互感
        to_remove = [k for k, m in self._mutuals.items()
                     if m.comp_id1 == comp_id or m.comp_id2 == comp_id]
        for k in to_remove:
            del self._mutuals[k]
            self._mutual_pool.release(k)
        self._mark_dirty()

    # 互感管理 

    def add_mutual(self, comp_id1: int, comp_id2: int, value) -> int:
        '''
        添加互感

        参数:
            comp_id1, comp_id2: 两个电感元件的 ID
            value: 互感值

        返回: 互感 ID
        '''
        if comp_id1 not in self._components or comp_id2 not in self._components:
            raise KeyError("互感引用的元件不存在")
        for cid in (comp_id1, comp_id2):
            if self._components[cid].type != 'L':
                raise TypeError(f"元件 {cid} 不是电感（type='{self._components[cid].type}'），互感只能施加在电感之间")
        mutual_id = self._mutual_pool.allocate()
        self._mutuals[mutual_id] = MutualInductance(comp_id1, comp_id2, value)
        self._mark_dirty()
        return mutual_id

    def remove_mutual(self, mutual_id: int):
        '''删除互感'''
        if mutual_id not in self._mutuals:
            raise KeyError(f"互感 {mutual_id} 不存在")
        del self._mutuals[mutual_id]
        self._mutual_pool.release(mutual_id)
        self._mark_dirty()

    # 物理磁通管理

    @staticmethod
    def _normalize_comp_loop(comp_loop):
        '''
        标准化元件回路: 返回 (canonical_tuple, direction)

        canonical_tuple: 规范化的元件 ID 元组 (所有旋转/反转中字典序最小的)
        direction: +1 (同向) 或 -1 (反向, 即原始回路与规范形式方向相反)
        '''
        n = len(comp_loop)
        if n == 2:
            a, b = comp_loop
            if a <= b:
                return (a, b), 1
            else:
                return (b, a), -1

        best = None
        best_dir = 1
        for i in range(n):
            fwd = tuple(comp_loop[(i + j) % n] for j in range(n))
            if best is None or fwd < best:
                best = fwd
                best_dir = 1
            bwd = tuple(comp_loop[(i - j) % n] for j in range(n))
            if bwd < best:
                best = bwd
                best_dir = -1
        return best, best_dir

    def add_physical_flux(self, comp_loop: list[int], flux, direction: int = None) -> int:
        '''
        添加物理外磁通 (相同物理回路自动合并)

        参数:
            comp_loop: 围成物理孔洞的元件 ID 列表 (按顺序排列)
                       顺序决定磁通正方向 (右手定则)
                       元件首尾相连形成闭合回路
            flux: 磁通值 (字符串自动转为 sympy 符号, 0 表示无磁通)
            direction: 仅在两元件回路时可用, 取值 +1 或 -1。
                       +1 (默认) 表示该回路的正方向为第一个元件从小序号节点到大序号节点;
                       -1 表示第一个元件从大序号节点到小序号节点 (即整个回路反向)。
                       若回路恰好由两个元件组成且未指定此参数, 将发出警告并默认使用 +1。

        返回: 磁通 ID (若合并到已有条目, 返回已有 ID)
        '''
        import warnings

        # 校验 direction 参数合法性
        if direction is not None:
            if direction not in (1, -1):
                raise ValueError(f"direction 参数必须为 +1 或 -1，当前值: {direction}")
            if len(comp_loop) != 2:
                raise ValueError("direction 参数仅适用于两元件回路")

        flux_val = sp.symbols(flux) if isinstance(flux, str) else flux

        # 校验: 元件不可重复
        if len(comp_loop) != len(set(comp_loop)):
            seen = set()
            dupes = [c for c in comp_loop if c in seen or seen.add(c)]
            raise ValueError(f"物理回路中元件不可重复, 重复的元件 ID: {dupes}")

        # 校验: 节点不可重复 (不可经过同一节点两次)
        if len(comp_loop) >= 2:
            dirs = self._trace_comp_loop_directions(comp_loop)
            if dirs is None:
                raise ValueError("物理回路元件不首尾相连或无法构成闭合回路")
            visited_nodes = [d[0] for d in dirs]  # 每条边的起始节点
            if len(visited_nodes) != len(set(visited_nodes)):
                seen = set()
                dupes = [n for n in visited_nodes if n in seen or seen.add(n)]
                raise ValueError(f"物理回路经过了重复节点: {dupes}")

            # 两元件回路: 方向不确定时发出警告并默认 +1
            if len(comp_loop) == 2 and direction is None:
                warnings.warn(
                    "两元件物理回路方向不确定，已默认为 +1 方向"
                    "（第一个元件从小序号节点到大序号节点）。"
                    "如需明确指定方向，请传入 direction=1 或 direction=-1。",
                    UserWarning, stacklevel=2
                )
                direction = 1

        user_dir = direction if direction is not None else 1
        canonical, norm_dir = self._normalize_comp_loop(comp_loop)
        canonical_flux = norm_dir * user_dir * flux_val  # 转换到规范方向

        # 检查是否已有相同物理回路, 若有则合并
        for fid, pf in self._physical_fluxes.items():
            if tuple(pf.comp_loop) == canonical:
                pf.flux = sp.nsimplify(pf.flux + canonical_flux)
                return fid

        # 新回路: 以规范形式存储
        flux_id = self._flux_pool.allocate()
        self._physical_fluxes[flux_id] = PhysicalFlux(list(canonical), canonical_flux)
        return flux_id

    def remove_physical_flux(self, flux_id: int):
        '''删除物理磁通'''
        if flux_id not in self._physical_fluxes:
            raise KeyError(f"物理磁通 {flux_id} 不存在")
        del self._physical_fluxes[flux_id]
        self._flux_pool.release(flux_id)

    # 基本回路磁通管理 (直接指定模式)

    def _validate_loop_index(self, loop_index: int) -> int:
        '''验证回路编号有效性, 返回对应的 chord_key'''
        self._ensure_built()
        num = self._m - self._nt
        if loop_index < 0 or loop_index >= num:
            raise IndexError(f"回路编号 {loop_index} 超出范围 [0, {num - 1}]")
        return self._loops[loop_index]['chord_key']

    def set_external_flux(self, loop_index: int, flux):
        '''
        直接设置某个基本回路的外磁通

        此方式与 add_physical_flux 互斥, 不可混用。
        调用此方法后即进入"直接指定模式", 需用 clear_external_fluxes() 清除后
        才能改用 add_physical_flux。

        参数:
            loop_index: 回路编号 (0-based, 对应 print_loops 中显示的回路编号)
            flux: 外磁通值, 可以是 sympy 表达式、字符串(自动转符号) 或数值
                  设为 0 表示该回路无外磁通
        '''
        chord_key = self._validate_loop_index(loop_index)
        if isinstance(flux, str):
            flux = sp.symbols(flux)
        self._fundamental_fluxes[chord_key] = flux

    def get_external_flux(self, loop_index: int):
        '''获取某个基本回路当前直接设置的外磁通'''
        chord_key = self._validate_loop_index(loop_index)
        return self._fundamental_fluxes.get(chord_key, 0)

    def clear_external_fluxes(self):
        '''清除所有直接设置的基本回路磁通'''
        self._fundamental_fluxes.clear()

    # 图构建

    def _rebuild(self):
        '''重新构建电路图、生成树、基本回路'''
        comp_ids = sorted(self._components.keys())
        components = [self._components[cid] for cid in comp_ids]

        # 1. 构建原始图
        G_orig = nx.MultiGraph()
        for i, comp in enumerate(components):
            G_orig.add_edge(comp.u, comp.v, key=i, type=comp.type,
                            value=comp.value, weight=comp.weight)
        G_orig.add_nodes_from(range(max(G_orig.nodes) + 1 if G_orig.nodes else 0))

        # 2. 最小生成树
        T_orig = nx.minimum_spanning_tree(G_orig, weight='weight')
        tree_orig_keys = set(k for _, _, k in T_orig.edges(keys=True))

        # 3. 区分树枝和连支
        tree_items = []
        chord_items = []
        for i, cid in enumerate(comp_ids):
            item = {'comp_id': cid, 'comp': components[i]}
            if i in tree_orig_keys:
                tree_items.append(item)
            else:
                chord_items.append(item)

        # 4. 重新编号: 树枝 0..nt-1, 连支 nt..m-1
        edge_map = {}
        comp_to_edge = {}

        def register_edges(items, graph, start_key):
            '''将一组元件注册为图的边, 返回下一个可用 key'''
            key = start_key
            for item in items:
                comp = item['comp']
                graph.add_edge(comp.u, comp.v, key=key,
                               type=comp.type, value=comp.value)
                edge_map[key] = {
                    'u': comp.u, 'v': comp.v,
                    'type': comp.type, 'value': comp.value,
                    'comp_id': item['comp_id'],
                }
                comp_to_edge[item['comp_id']] = key
                key += 1
            return key

        T_new = nx.MultiGraph()
        T_new.add_nodes_from(G_orig.nodes)
        nt = register_edges(tree_items, T_new, 0)

        G_new = T_new.copy()
        m = register_edges(chord_items, G_new, nt)

        # 5. 互感映射
        mutual_dict = {}
        for mut in self._mutuals.values():
            ek1 = comp_to_edge[mut.comp_id1]
            ek2 = comp_to_edge[mut.comp_id2]
            mutual_dict[tuple(sorted((ek1, ek2)))] = mut.value

        # 6. 计算基本回路
        tree_simple = nx.Graph()
        for u, v, k in T_new.edges(keys=True):
            tree_simple.add_edge(u, v, key=k)

        loops = []
        for chord_key in range(nt, m):
            chord_info = edge_map[chord_key]
            u_chord, v_chord = chord_info['u'], chord_info['v']
            path_nodes = nx.shortest_path(tree_simple, u_chord, v_chord)

            tree_branch_keys = []
            for j in range(len(path_nodes) - 1):
                n1, n2 = path_nodes[j], path_nodes[j + 1]
                tree_branch_keys.append(tree_simple.edges[n1, n2]['key'])

            loops.append({
                'chord_key': chord_key,
                'chord_info': chord_info,
                'tree_branch_keys': tree_branch_keys,
                'path_nodes': path_nodes,
            })

        # 保存
        self._T = T_new
        self._G = G_new
        self._edge_map = edge_map
        self._nt = nt
        self._m = m
        self._mutual_dict = mutual_dict
        self._loops = loops
        self._comp_to_edge = comp_to_edge

    #  物理磁通分解 

    def _physical_loop_to_signed_vector(self, phys_flux: PhysicalFlux) -> list:
        '''
        将物理回路 (元件 ID 序列) 转换为有符号边向量

        复用 _trace_comp_loop_directions 追踪节点路径,
        再根据遍历方向标记 +1 (u→v) 或 -1 (v→u)
        '''
        comp_loop = phys_flux.comp_loop
        dirs = self._trace_comp_loop_directions(comp_loop)
        if dirs is None:
            raise ValueError("物理回路无法构成闭合回路")

        p = [0] * self._m
        for cid, (from_n, _) in zip(comp_loop, dirs):
            ek = self._comp_to_edge[cid]
            ei = self._edge_map[ek]
            p[ek] = 1 if from_n == ei['u'] else -1
        return p

    def _compute_external_fluxes(self) -> dict:
        '''
        计算各基本回路的外磁通

        支持两种模式 (不可混用):
        1. 物理回路模式: 从 add_physical_flux 提供的物理回路分解得到
        2. 直接指定模式: 从 set_external_flux 直接设置, 未指定的回路磁通为 0
        '''
        num_chords = self._m - self._nt

        if num_chords == 0:
            return {}

        has_physical = bool(self._physical_fluxes)
        has_fundamental = bool(self._fundamental_fluxes)

        if has_physical and has_fundamental:
            raise ValueError(
                "不可同时使用物理回路磁通 (add_physical_flux) 和"
                "直接指定基本回路磁通 (set_external_flux)。"
                "请只使用其中一种方式。"
            )

        chord_keys = [loop['chord_key'] for loop in self._loops]

        # 模式 B: 直接指定基本回路磁通
        if has_fundamental:
            return {ck: self._fundamental_fluxes.get(ck, 0) for ck in chord_keys}

        # 模式 A: 物理回路分解
        phys_list = list(self._physical_fluxes.values())
        num_phys = len(phys_list)

        if num_phys != num_chords:
            raise ValueError(
                f"需要恰好 {num_chords} 个线性无关的物理回路来确定外磁通，"
                f"当前已提供 {num_phys} 个。"
                f"即使磁通为 0 的回路也必须添加 (使用 add_physical_flux(..., 0))。"
            )

        # 构建分解矩阵 A 和磁通向量 b
        A = sp.zeros(num_phys, num_chords)
        b = sp.zeros(num_phys, 1)

        for i, pf in enumerate(phys_list):
            p = self._physical_loop_to_signed_vector(pf)
            for j, ck in enumerate(chord_keys):
                A[i, j] = p[ck]
            b[i] = pf.flux

        # 检查物理回路是否线性无关
        if A.rank() < num_chords:
            raise ValueError(
                "提供的物理回路线性相关，无法唯一确定各基本回路的外磁通。"
                "请确保每个物理回路对应不同的独立孔洞。"
            )

        # 解线性方程组: A · φ_ext = b
        x = A.LUsolve(b)

        return {ck: sp.nsimplify(x[j]) for j, ck in enumerate(chord_keys)}

    #  基本割集矩阵 

    def _compute_fundamental_cut_matrix(self):
        '''计算基本割集矩阵 Q_f, 行对应树枝编号'''
        T, edge_map = self._T, self._edge_map
        m, nt = self._m, self._nt
        Q_list = [[0] * m for _ in range(nt)]

        tree_graph = nx.Graph()
        for u, v, k in T.edges(keys=True):
            tree_graph.add_edge(u, v, key=k)

        for row_key in range(nt):
            Q_list[row_key][row_key] = 1
            u_tree, v_tree = edge_map[row_key]['u'], edge_map[row_key]['v']

            T_temp = tree_graph.copy()
            T_temp.remove_edge(u_tree, v_tree)
            comp_u = set(nx.node_connected_component(T_temp, u_tree))

            for col_key in range(nt, m):
                u_chord = edge_map[col_key]['u']
                v_chord = edge_map[col_key]['v']
                u_in_U = u_chord in comp_u
                v_in_U = v_chord in comp_u

                if u_in_U and not v_in_U:
                    Q_list[row_key][col_key] = 1
                elif not u_in_U and v_in_U:
                    Q_list[row_key][col_key] = -1

        return sp.Matrix(Q_list)

    #  参数矩阵 

    def _build_parameter_matrices(self):
        '''构建 D_C (电容), L_plus (电感逆), D_J (约瑟夫森) 矩阵'''
        edge_map = self._edge_map
        mutual_dict = self._mutual_dict
        num_edges = len(edge_map)

        # 电容矩阵 D_C
        Dc_list = [0] * num_edges
        for k, info in edge_map.items():
            if info['type'] == 'C':
                Dc_list[k] = info['value']
        D_C = sp.diag(*Dc_list)

        # 约瑟夫森矩阵 D_J
        Dj_list = [0] * num_edges
        for k, info in edge_map.items():
            if info['type'] == 'JJ':
                Dj_list[k] = info['value']
        D_J = sp.diag(*Dj_list)

        # 电感逆矩阵 L_plus
        L_indices = [k for k, info in edge_map.items() if info['type'] == 'L']

        if not L_indices:
            L_plus = sp.zeros(num_edges, num_edges)
        else:
            n_L = len(L_indices)
            L_sub = sp.zeros(n_L, n_L)
            local_to_global = {i: k for i, k in enumerate(L_indices)}

            for i in range(n_L):
                k_i = local_to_global[i]
                L_sub[i, i] = edge_map[k_i]['value']
                for j in range(i + 1, n_L):
                    k_j = local_to_global[j]
                    key_tuple = tuple(sorted((k_i, k_j)))
                    if key_tuple in mutual_dict:
                        val = mutual_dict[key_tuple]
                        L_sub[i, j] = val
                        L_sub[j, i] = val

            try:
                Gamma_sub = L_sub.inv()
            except ValueError:
                raise ValueError("电感矩阵奇异，无法求逆。请检查互感参数或电路拓扑。")

            L_plus = sp.zeros(num_edges, num_edges)
            for i in range(n_L):
                for j in range(n_L):
                    L_plus[local_to_global[i], local_to_global[j]] = Gamma_sub[i, j]

        return D_C, L_plus, D_J

    # 公开属性 

    @property
    def nt(self):
        '''树枝数 (自由度)'''
        self._ensure_built()
        return self._nt

    @property
    def m(self):
        '''总支路数'''
        self._ensure_built()
        return self._m

    @property
    def num_loops(self):
        '''基本回路数 (= 连支数 = m - nt)'''
        self._ensure_built()
        return self._m - self._nt

    @property
    def edge_map(self):
        '''边映射表 {key: {u, v, type, value, comp_id}}'''
        self._ensure_built()
        return self._edge_map

    @property
    def loops(self):
        '''基本回路信息列表'''
        self._ensure_built()
        return self._loops

    @property
    def F_C(self):
        '''基本割集矩阵 (带缓存)'''
        self._ensure_built()
        if self._F_C_cache is None:
            self._F_C_cache = self._compute_fundamental_cut_matrix()
        return self._F_C_cache

    @property
    def components(self):
        '''当前元件字典 {comp_id: Component}'''
        return dict(self._components)

    @property
    def physical_fluxes(self):
        '''当前物理磁通字典 {flux_id: PhysicalFlux}'''
        return dict(self._physical_fluxes)

    @property
    def fundamental_fluxes(self):
        '''当前直接设置的基本回路磁通 {loop_index: flux}'''
        self._ensure_built()
        chord_to_index = {loop['chord_key']: i for i, loop in enumerate(self._loops)}
        return {chord_to_index[ck]: flux for ck, flux in self._fundamental_fluxes.items()}

    @property
    def mutuals(self):
        '''当前互感字典 {mutual_id: MutualInductance}'''
        return dict(self._mutuals)

    # 打印信息

    def print_nodes(self):
        '''打印所有活跃节点及其引用计数'''
        if not self._node_refcount:
            print("电路中没有节点")
            return
        print(f"共 {len(self._node_refcount)} 个活跃节点:")
        for nid in sorted(self._node_refcount.keys()):
            print(f"  节点 {nid}: 被 {self._node_refcount[nid]} 个元件引用")

    def print_components(self):
        '''打印所有元件'''
        if not self._components:
            print("电路中没有元件")
            return
        print(f"共 {len(self._components)} 个元件:")
        for cid in sorted(self._components.keys()):
            comp = self._components[cid]
            print(f"  元件 {cid}: {comp.type} ({comp.u}-{comp.v}) 值: {comp.value}")

    def print_edges(self):
        '''打印所有边 (树枝 + 连支) 的信息'''
        self._ensure_built()
        print(f"树枝数: {self._nt}, 连支数: {self._m - self._nt}, 总支路数: {self._m}")
        for k in sorted(self._edge_map.keys()):
            v = self._edge_map[k]
            role = "树枝" if k < self._nt else "连支"
            print(f"  Key {k} [{role}]: {v['type']} ({v['u']}->{v['v']}) "
                  f"值: {v['value']} [元件ID: {v['comp_id']}]")

    def print_loops(self):
        '''打印所有基本回路信息及当前外磁通'''
        self._ensure_built()
        try:
            external_fluxes = self._compute_external_fluxes()
            flux_ready = True
        except ValueError as e:
            external_fluxes = {loop['chord_key']: 0 for loop in self._loops}
            flux_ready = False
            flux_err_msg = str(e)

        # 预计算物理回路的有符号边向量 (用于显示磁通来源)
        phys_vectors = {}
        for fid, pf in self._physical_fluxes.items():
            try:
                phys_vectors[fid] = self._physical_loop_to_signed_vector(pf)
            except (KeyError, ValueError):
                pass

        print(f"电路共有 {self.num_loops} 个基本回路 "
              f"(树枝数: {self.nt}, 总支路数: {self.m})")
        print()

        for i, loop in enumerate(self._loops):
            chord_key = loop['chord_key']
            chord = loop['chord_info']
            flux = external_fluxes.get(chord_key, 0)

            print(f"回路 {i}: (由连支 Key {chord_key} 闭合)")
            # 回路正方向: 连支 u→v (+1), 然后树枝从 v 沿树回到 u
            print(f"  连支: {chord['type']} ({chord['u']}->{chord['v']}), "
                  f"值: {chord['value']} [元件ID: {chord['comp_id']}], "
                  f"遍历: {chord['u']}→{chord['v']} (+1)")
            # 树枝路径: 反转 path_nodes (原 path 是 u→v, 正方向树枝部分是 v→u)
            tree_path = list(reversed(loop['path_nodes']))
            tree_keys = list(reversed(loop['tree_branch_keys']))
            full_path = [chord['u']] + tree_path
            path_str = ' -> '.join(map(str, full_path))
            print(f"  回路路径: {path_str}")
            print(f"  包含树枝:")
            for j, bk in enumerate(tree_keys):
                info = self._edge_map[bk]
                fn = tree_path[j]
                tn = tree_path[j + 1]
                sign = "+1" if fn == info['u'] else "-1"
                print(f"    Key {bk}: {info['type']} ({info['u']}->{info['v']}), "
                      f"值: {info['value']} [元件ID: {info['comp_id']}], "
                      f"遍历: {fn}→{tn} ({sign})")

            # 显示外磁通及其来源
            if flux_ready:
                print(f"  外磁通: {flux}")
                if self._fundamental_fluxes:
                    # 直接指定模式
                    if chord_key in self._fundamental_fluxes:
                        print(f"  磁通来源: 直接指定 (set_external_flux)")
                    else:
                        print(f"  磁通来源: 默认值 0 (未指定)")
                elif phys_vectors:
                    sources = []
                    for fid, pvec in phys_vectors.items():
                        coeff = pvec[chord_key]
                        if coeff != 0:
                            pf = self._physical_fluxes[fid]
                            sign = "+" if coeff > 0 else ""
                            sources.append(f"{sign}{coeff} × {pf.flux} (磁通ID {fid})")
                    if sources:
                        print(f"  磁通来源: {', '.join(sources)}")
            else:
                print(f"  外磁通: 未确定")
            print()

        if not flux_ready:
            print(f"[注意] {flux_err_msg}")

    def _trace_comp_loop_directions(self, comp_loop):
        '''追踪元件回路中每个元件的遍历方向, 返回 [(from_node, to_node), ...]'''
        n = len(comp_loop)
        if n < 2:
            return None

        comps = []
        for cid in comp_loop:
            if cid not in self._components:
                return None
            comps.append(self._components[cid])

        e0_nodes = {comps[0].u, comps[0].v}
        e1_nodes = {comps[1].u, comps[1].v}
        shared = e0_nodes & e1_nodes
        if not shared:
            return None

        if len(shared) == 2:
            # 两元件并联 (共享两个节点): 按第一个元件的正方向 (u→v) 遍历
            start_node = comps[0].u
        else:
            shared_node = shared.pop()
            start_node = (e0_nodes - {shared_node}).pop()

        directions = []
        current = start_node
        for comp in comps:
            if current == comp.u:
                directions.append((comp.u, comp.v))
                current = comp.v
            elif current == comp.v:
                directions.append((comp.v, comp.u))
                current = comp.u
            else:
                return None

        if current != start_node:
            return None
        return directions

    def print_physical_fluxes(self):
        '''打印所有物理磁通的设置'''
        if not self._physical_fluxes:
            print("未设置任何物理磁通")
            return

        print(f"共 {len(self._physical_fluxes)} 个物理磁通:")
        for fid, pf in self._physical_fluxes.items():
            comp_str = ' -> '.join(str(c) for c in pf.comp_loop)
            # 追踪遍历方向
            dirs = self._trace_comp_loop_directions(pf.comp_loop)
            if dirs:
                details = []
                for cid, (fn, tn) in zip(pf.comp_loop, dirs):
                    comp = self._components[cid]
                    details.append(f"{comp.type}({fn}→{tn})")
                detail_str = ' -> '.join(details)
            else:
                details = []
                for cid in pf.comp_loop:
                    if cid in self._components:
                        comp = self._components[cid]
                        details.append(f"{comp.type}({comp.u}-{comp.v})")
                    else:
                        details.append("[已删除]")
                detail_str = ' -> '.join(details)
            print(f"  磁通 {fid}: Φ = {pf.flux}")
            print(f"    回路元件: [{comp_str}] = {detail_str}")

    # 哈密顿量 

    def hamiltonian(self):
        '''
        计算并返回哈密顿量

        返回:
            H_total: 总哈密顿量 sympy 表达式
            info: dict 包含中间结果
                - M:          质量矩阵 (电容相关)
                - phi_t:      树枝磁通变量 (广义坐标)
                - q_t:        树枝电荷变量 (广义动量)
                - ext_fluxes: 实际使用的外磁通符号列表
                - Phi_vec:    全支路磁通表达式
                - F_C:        基本割集矩阵
                - D_C:        电容矩阵
                - L_plus:     电感逆矩阵
                - D_J:        约瑟夫森矩阵
        '''
        self._ensure_built()
        F_C = self.F_C
        D_C, L_plus, D_J = self._build_parameter_matrices()
        external_fluxes = self._compute_external_fluxes()

        nt = self._nt
        m = self._m

        # 1. 定义动态变量
        phi_t = sp.Matrix([sp.symbols(f'Phi_t_{i}') for i in range(nt)])
        q_t = sp.Matrix([sp.symbols(f'Q_t_{i}') for i in range(nt)])

        # 2. 构建外磁通向量 (由物理磁通分解或直接指定)
        phi_ext_list = []
        ext_flux_symbols = set()

        for k in range(m):
            if k < nt:
                phi_ext_list.append(0)
            else:
                flux = external_fluxes.get(k, 0)
                phi_ext_list.append(flux)
                if isinstance(flux, sp.Basic):
                    ext_flux_symbols.update(flux.free_symbols)

        phi_ext_vec = sp.Matrix(phi_ext_list)
        ext_flux_vars = sorted(ext_flux_symbols, key=lambda s: str(s))

        # 3. 全支路磁通: Phi_all = F_C^T * Phi_tree + Phi_ext
        Phi_vec = F_C.T * phi_t + phi_ext_vec

        # 4. 动能: H_kin = ½ q^T M^{-1} q
        M = F_C * D_C * F_C.T
        try:
            M_inv = M.inv()
        except ValueError:
            raise ValueError("质量矩阵 M = F_C·D_C·F_C^T 奇异，无法求逆。请检查电容参数。")
        H_kin = (sp.Rational(1, 2) * q_t.T * M_inv * q_t)[0]

        # 5. 电感势能: H_L = ½ Phi^T L_plus Phi
        H_mag_lin = (sp.Rational(1, 2) * Phi_vec.T * L_plus * Phi_vec)[0]

        # 6. 约瑟夫森势能: H_JJ = -Σ Ej cos(Phi_k)
        H_jj = 0
        for k in range(m):
            coeff = D_J[k, k]
            if coeff != 0:
                phi_k = Phi_vec[k]
                H_jj -= coeff * sp.cos(phi_k)

        H_total = H_kin + H_mag_lin + H_jj

        info = {
            'M': M,
            'phi_t': phi_t,
            'q_t': q_t,
            'ext_fluxes': ext_flux_vars,
            'Phi_vec': Phi_vec,
            'F_C': F_C,
            'D_C': D_C,
            'L_plus': L_plus,
            'D_J': D_J,
        }

        return H_total, info


