import sympy as sp
import networkx as nx


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
        '''生成树权重 (最小权重优先进入树枝)'''
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
    超导电路量子化对象 (动态构建版)

    使用流程:
        1. 创建空电路:      circuit = Circuit()
        2. 添加元件:        cid = circuit.add_component(u, v, type, value)
        3. 添加互感 (可选): mid = circuit.add_mutual(cid1, cid2, value)
        4. 添加物理磁通:    fid = circuit.add_physical_flux(comp_loop, flux)
        5. 查看信息:        circuit.print_edges() / circuit.print_loops()
        6. 计算哈密顿量:    H, info = circuit.hamiltonian()

    元件/互感/磁通可随时增删, 所有计算按需自动重建图结构
    '''

    def __init__(self):
        self._components = {}        # comp_id -> Component
        self._mutuals = {}           # mutual_id -> MutualInductance
        self._physical_fluxes = {}   # flux_id -> PhysicalFlux
        self._next_comp_id = 0
        self._next_mutual_id = 0
        self._next_flux_id = 0
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

    def _ensure_built(self):
        '''确保图结构已构建, 按需重建'''
        if self._dirty:
            if not self._components:
                raise ValueError("电路中没有元件，无法构建")
            self._rebuild()
            self._dirty = False

    # ==================== 元件管理 ====================

    def add_component(self, u: int, v: int, tp: str, value) -> int:
        '''
        添加元件

        参数:
            u, v: 节点编号
            tp: 元件类型 ('C', 'L', 'JJ')
            value: 参数值 (字符串自动转为 sympy 符号)

        返回: 元件 ID (用于后续引用)
        '''
        comp_id = self._next_comp_id
        self._next_comp_id += 1
        self._components[comp_id] = Component(u, v, tp, value)
        self._mark_dirty()
        return comp_id

    def remove_component(self, comp_id: int):
        '''删除元件 (关联的互感会一并删除)'''
        if comp_id not in self._components:
            raise KeyError(f"元件 {comp_id} 不存在")
        del self._components[comp_id]
        # 移除引用了该元件的互感
        to_remove = [k for k, m in self._mutuals.items()
                     if m.comp_id1 == comp_id or m.comp_id2 == comp_id]
        for k in to_remove:
            del self._mutuals[k]
        self._mark_dirty()

    # ==================== 互感管理 ====================

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
        mutual_id = self._next_mutual_id
        self._next_mutual_id += 1
        self._mutuals[mutual_id] = MutualInductance(comp_id1, comp_id2, value)
        self._mark_dirty()
        return mutual_id

    def remove_mutual(self, mutual_id: int):
        '''删除互感'''
        if mutual_id not in self._mutuals:
            raise KeyError(f"互感 {mutual_id} 不存在")
        del self._mutuals[mutual_id]
        self._mark_dirty()

    # ==================== 物理磁通管理 ====================

    def add_physical_flux(self, comp_loop: list[int], flux) -> int:
        '''
        添加物理外磁通

        参数:
            comp_loop: 围成物理孔洞的元件 ID 列表 (按顺序排列)
                       顺序决定磁通正方向 (右手定则)
                       元件首尾相连形成闭合回路
            flux: 磁通值 (字符串自动转为 sympy 符号, 0 表示无磁通)

        返回: 磁通 ID
        '''
        flux_id = self._next_flux_id
        self._next_flux_id += 1
        self._physical_fluxes[flux_id] = PhysicalFlux(comp_loop, flux)
        return flux_id

    def remove_physical_flux(self, flux_id: int):
        '''删除物理磁通'''
        if flux_id not in self._physical_fluxes:
            raise KeyError(f"物理磁通 {flux_id} 不存在")
        del self._physical_fluxes[flux_id]

    # ==================== 内部: 图构建 ====================

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

        T_new = nx.MultiGraph()
        T_new.add_nodes_from(G_orig.nodes)
        current_key = 0

        for item in tree_items:
            comp = item['comp']
            T_new.add_edge(comp.u, comp.v, key=current_key,
                           type=comp.type, value=comp.value)
            edge_map[current_key] = {
                'u': comp.u, 'v': comp.v,
                'type': comp.type, 'value': comp.value,
                'comp_id': item['comp_id'],
            }
            comp_to_edge[item['comp_id']] = current_key
            current_key += 1

        nt = current_key

        G_new = T_new.copy()
        for item in chord_items:
            comp = item['comp']
            G_new.add_edge(comp.u, comp.v, key=current_key,
                           type=comp.type, value=comp.value)
            edge_map[current_key] = {
                'u': comp.u, 'v': comp.v,
                'type': comp.type, 'value': comp.value,
                'comp_id': item['comp_id'],
            }
            comp_to_edge[item['comp_id']] = current_key
            current_key += 1

        m = current_key

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

    # ==================== 内部: 物理磁通分解 ====================

    def _physical_loop_to_signed_vector(self, phys_flux: PhysicalFlux) -> list:
        '''
        将物理回路 (元件 ID 序列) 转换为有符号边向量

        算法:
            1. 根据元件顺序追踪节点路径
            2. 每条边按遍历方向标记 +1 (u→v) 或 -1 (v→u)

        约定: 边 k 的正方向为 u→v (u < v)
        '''
        comp_loop = phys_flux.comp_loop
        m = self._m

        for cid in comp_loop:
            if cid not in self._comp_to_edge:
                raise KeyError(f"物理回路引用的元件 {cid} 不存在或已被删除")

        edge_keys = [self._comp_to_edge[cid] for cid in comp_loop]
        edge_infos = [self._edge_map[ek] for ek in edge_keys]
        n = len(comp_loop)

        if n < 2:
            raise ValueError("物理回路至少需要 2 个元件")

        # 找第一条边和第二条边的共享节点, 确定起始方向
        e0_nodes = {edge_infos[0]['u'], edge_infos[0]['v']}
        e1_nodes = {edge_infos[1]['u'], edge_infos[1]['v']}
        shared = e0_nodes & e1_nodes

        if not shared:
            raise ValueError(
                f"元件 {comp_loop[0]} ({edge_infos[0]['u']}-{edge_infos[0]['v']}) "
                f"与元件 {comp_loop[1]} ({edge_infos[1]['u']}-{edge_infos[1]['v']}) 不相邻"
            )

        # 从第一条边的非共享端点出发
        shared_node = shared.pop()
        start_node = (e0_nodes - {shared_node}).pop()

        # 追踪路径, 构建有符号边向量
        p = [0] * m
        current_node = start_node

        for i in range(n):
            ek = edge_keys[i]
            ei = edge_infos[i]
            u, v = ei['u'], ei['v']

            if current_node == u:
                p[ek] = 1        # 正方向 u→v
                current_node = v
            elif current_node == v:
                p[ek] = -1       # 反方向 v→u
                current_node = u
            else:
                raise ValueError(
                    f"元件 {comp_loop[i]} 与路径不连续 "
                    f"(当前节点 {current_node}, 元件端点 {u}-{v})"
                )

        if current_node != start_node:
            raise ValueError(
                f"物理回路未闭合: 起点 {start_node}, 终点 {current_node}")

        return p

    def _compute_external_fluxes(self) -> dict:
        '''
        从所有物理磁通计算各基本回路的外磁通

        原理:
            任何回路 p 都可以唯一分解为基本回路的线性组合
            分解系数 c_i = p 在连支 k_i 上的分量 (±1 或 0)
            这是因为基本回路矩阵在连支列构成单位阵
        '''
        external_fluxes = {}
        for loop in self._loops:
            external_fluxes[loop['chord_key']] = 0

        for phys_flux in self._physical_fluxes.values():
            p = self._physical_loop_to_signed_vector(phys_flux)

            for loop in self._loops:
                chord_key = loop['chord_key']
                coeff = p[chord_key]
                if coeff != 0:
                    external_fluxes[chord_key] += coeff * phys_flux.flux

        return external_fluxes

    # ==================== 内部: 基本割集矩阵 ====================

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

    # ==================== 内部: 参数矩阵 ====================

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
            except:
                print("警告：电感矩阵奇异，可能包含无电感回路或未定义的互感。")
                Gamma_sub = sp.zeros(n_L, n_L)

            L_plus = sp.zeros(num_edges, num_edges)
            for i in range(n_L):
                for j in range(n_L):
                    L_plus[local_to_global[i], local_to_global[j]] = Gamma_sub[i, j]

        return D_C, L_plus, D_J

    # ==================== 公开属性 ====================

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

    # ==================== 打印信息 ====================

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
        external_fluxes = self._compute_external_fluxes()

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
            print(f"  连支: {chord['type']} ({chord['u']}->{chord['v']}), "
                  f"值: {chord['value']} [元件ID: {chord['comp_id']}]")
            path_str = ' -> '.join(map(str, loop['path_nodes']))
            print(f"  回路路径: {path_str} -> {loop['path_nodes'][0]}")
            print(f"  包含树枝:")
            for bk in loop['tree_branch_keys']:
                info = self._edge_map[bk]
                print(f"    Key {bk}: {info['type']} ({info['u']}->{info['v']}), "
                      f"值: {info['value']} [元件ID: {info['comp_id']}]")

            # 显示外磁通及其来源
            print(f"  外磁通: {flux}")
            if phys_vectors:
                sources = []
                for fid, pvec in phys_vectors.items():
                    coeff = pvec[chord_key]
                    if coeff != 0:
                        pf = self._physical_fluxes[fid]
                        sign = "+" if coeff > 0 else ""
                        sources.append(f"{sign}{coeff} × {pf.flux} (磁通ID {fid})")
                if sources:
                    print(f"  磁通来源: {', '.join(sources)}")
            print()

    def print_physical_fluxes(self):
        '''打印所有物理磁通的设置'''
        if not self._physical_fluxes:
            print("未设置任何物理磁通")
            return

        print(f"共 {len(self._physical_fluxes)} 个物理磁通:")
        for fid, pf in self._physical_fluxes.items():
            comp_str = ' -> '.join(str(c) for c in pf.comp_loop)
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

    # ==================== 哈密顿量 ====================

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

        # 2. 构建外磁通向量 (由物理磁通自动计算)
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
            H_kin = (sp.Rational(1, 2) * q_t.T * M_inv * q_t)[0]
        except:
            M_inv_sym = sp.Symbol("M^{-1}")
            H_kin = (sp.Rational(1, 2) * q_t.T * M_inv_sym * q_t)[0]

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


# ==================== 测试示例 ====================

if __name__ == "__main__":

    # ---- 创建空电路, 逐步添加元件 ----
    circuit = Circuit()

    jj1 = circuit.add_component(0, 2, 'JJ', 'EJ1')
    jj2 = circuit.add_component(1, 2, 'JJ', 'EJ2')
    L   = circuit.add_component(0, 1, 'L',  'L')
    C1  = circuit.add_component(0, 2, 'C',  'C1')
    C2  = circuit.add_component(1, 2, 'C',  'C2')

    print("=" * 50)
    print("元件列表:")
    print("=" * 50)
    circuit.print_components()

    print()
    print("=" * 50)
    print("重排后的边信息:")
    print("=" * 50)
    circuit.print_edges()

    print()
    print("=" * 50)
    print("基本回路 (无外磁通):")
    print("=" * 50)
    circuit.print_loops()

    # ---- 添加物理外磁通 ----
    # 物理孔洞由 JJ1, L, JJ2 围成, 加外磁通 Phi_e
    fid = circuit.add_physical_flux(comp_loop=[jj1, L, jj2], flux='Phi_e')

    print("=" * 50)
    print("设置物理磁通后:")
    print("=" * 50)
    circuit.print_physical_fluxes()
    print()
    circuit.print_loops()

    # ---- 计算哈密顿量 ----
    print("=" * 50)
    print("基本割集矩阵 F_C:")
    print("=" * 50)
    sp.pprint(circuit.F_C)

    H, info = circuit.hamiltonian()

    print(f"\n广义坐标 (树枝磁通):")
    sp.pprint(info['phi_t'])

    print(f"\n外磁通变量:")
    if info['ext_fluxes']:
        sp.pprint(info['ext_fluxes'])
    else:
        print("无")

    print(f"\n全支路磁通 Phi_vec:")
    sp.pprint(info['Phi_vec'])

    print(f"\n哈密顿量 H:")
    sp.pprint(sp.expand(H))
