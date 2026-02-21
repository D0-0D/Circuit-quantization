import sympy as sp
import networkx as nx
import numpy as np

class Component:
    '''定义元件, 用于校验输入'''
    def __init__(self, u: int, v: int, tp: str, value):
        self.u = u
        self.v = v
        self.type = tp.upper()
        valid_types = {'C', 'L', 'JJ'}
        if self.type not in valid_types:
            raise ValueError(f"元件类型 '{self.type}' 不在支持列表 {valid_types} 中")
        if isinstance(value, str):
            self.value = sp.symbols(value)
        else:
            self.value = value
    @property
    def weight(self):
        # 电容优先进入生成树 (权重小)
        return 1 if self.type == 'C' else (2 if self.type == 'JJ' else 3)

class MutualInductance:
    '''定义互感系数'''
    def __init__(self, L1_id: int, L2_id: int, value):
        self.L1_id = L1_id
        self.L2_id = L2_id
        self.value = sp.symbols(value) if isinstance(value, str) else value

class Circuit:
    '''
    超导电路量子化对象

    使用流程:
        1. 创建电路对象:  circuit = Circuit(components, mutuals)
        2. 查看基本回路:  circuit.print_loops()
        3. 编辑外磁通:    circuit.set_external_flux(loop_index, flux)
        4. 计算哈密顿量:  H, info = circuit.hamiltonian()
    '''

    def __init__(self, components: list[Component], mutuals: list[MutualInductance] | None = None):
        self._components = components
        self._mutuals = mutuals

        # 构建图和生成树
        self._build_graph()
        # 计算基本回路信息
        self._compute_loops()
        # 初始化默认外磁通符号
        self._init_default_fluxes()
        # 缓存
        self._F_C_cache = None

    # ==================== 图构建 ====================

    def _build_graph(self):
        '''构建电路图和生成树，重新编号边'''
        components = self._components
        mutuals = self._mutuals

        # 1. 初始构建图
        G_orig = nx.MultiGraph()
        edge_map_orig = []
        for i, comp in enumerate(components):
            G_orig.add_edge(comp.u, comp.v, key=i, type=comp.type,
                            value=comp.value, weight=comp.weight)
            edge_map_orig.append({
                'u': comp.u, 'v': comp.v,
                'type': comp.type, 'value': comp.value, 'orig_id': i
            })

        G_orig.add_nodes_from(range(max(G_orig.nodes) + 1 if G_orig.nodes else 0))

        # 2. 生成树 (最小权重优先：C -> JJ -> L)
        T_orig = nx.minimum_spanning_tree(G_orig, weight='weight')

        # 3. 区分树枝和连支，重新编号
        tree_edges = []
        chord_edges = []
        tree_orig_keys = set(k for _, _, k in T_orig.edges(keys=True))

        for k, info in enumerate(edge_map_orig):
            if k in tree_orig_keys:
                tree_edges.append(info)
            else:
                chord_edges.append(info)

        # 重新构建映射
        new_edge_map = {}
        old_to_new_key = {}
        current_new_key = 0

        # 处理树枝
        T_new = nx.MultiGraph()
        T_new.add_nodes_from(G_orig.nodes)
        for info in tree_edges:
            u, v = info['u'], info['v']
            if u > v: u, v = v, u
            T_new.add_edge(u, v, key=current_new_key,
                           type=info['type'], value=info['value'])
            new_edge_map[current_new_key] = {
                'u': u, 'v': v, 'type': info['type'], 'value': info['value']
            }
            old_to_new_key[info['orig_id']] = current_new_key
            current_new_key += 1

        # 处理连支
        G_new = T_new.copy()
        for info in chord_edges:
            u, v = info['u'], info['v']
            if u > v: u, v = v, u
            G_new.add_edge(u, v, key=current_new_key,
                           type=info['type'], value=info['value'])
            new_edge_map[current_new_key] = {
                'u': u, 'v': v, 'type': info['type'], 'value': info['value']
            }
            old_to_new_key[info['orig_id']] = current_new_key
            current_new_key += 1

        # 处理互感字典
        new_mutual_dict = {}
        if mutuals:
            for m in mutuals:
                nk1 = old_to_new_key[m.L1_id]
                nk2 = old_to_new_key[m.L2_id]
                new_mutual_dict[tuple(sorted((nk1, nk2)))] = m.value

        self._T = T_new
        self._G = G_new
        self._edge_map = new_edge_map
        self._mutual_dict = new_mutual_dict
        self._nt = T_new.number_of_edges()
        self._m = G_new.number_of_edges()

    # ==================== 基本回路 ====================

    def _compute_loops(self):
        '''计算每个连支对应的基本回路 (在树中查找路径)'''
        T = self._T
        edge_map = self._edge_map
        nt = self._nt
        m = self._m

        # 多重图转单重图用于路径查找
        tree_simple = nx.Graph()
        for u, v, k in T.edges(keys=True):
            tree_simple.add_edge(u, v, key=k)

        loops = []
        for chord_key in range(nt, m):
            chord_info = edge_map[chord_key]
            u_chord, v_chord = chord_info['u'], chord_info['v']

            # 在树中找 u_chord 到 v_chord 的唯一路径
            path_nodes = nx.shortest_path(tree_simple, u_chord, v_chord)

            # 提取路径上的树枝 key
            tree_branch_keys = []
            for i in range(len(path_nodes) - 1):
                n1, n2 = path_nodes[i], path_nodes[i + 1]
                tree_branch_keys.append(tree_simple.edges[n1, n2]['key'])

            loops.append({
                'chord_key': chord_key,
                'chord_info': chord_info,
                'tree_branch_keys': tree_branch_keys,
                'path_nodes': path_nodes,
            })

        self._loops = loops

    def _init_default_fluxes(self):
        '''为每个基本回路初始化默认外磁通符号 Phi_ext_<连支编号>'''
        self._external_fluxes = {}
        for loop in self._loops:
            chord_key = loop['chord_key']
            self._external_fluxes[chord_key] = sp.symbols(f'Phi_ext_{chord_key}')

    # ==================== 基本割集矩阵 ====================

    def _compute_fundamental_cut_matrix(self):
        '''计算基本割集矩阵 Q_f, 行对应树枝编号'''
        T, G, edge_map = self._T, self._G, self._edge_map
        m = self._m
        nt = self._nt
        Q_list = [[0] * m for _ in range(nt)]

        # 多重图转单重
        tree_graph = nx.Graph()
        for u, v, k in T.edges(keys=True):
            tree_graph.add_edge(u, v, key=k)

        # 逐行生成矩阵
        for row_key in range(nt):
            Q_list[row_key][row_key] = 1
            u_tree, v_tree = edge_map[row_key]['u'], edge_map[row_key]['v']

            T_temp = tree_graph.copy()
            T_temp.remove_edge(u_tree, v_tree)

            # 找出生成树剪开枝后与小编号节点相连的连通分量
            comp_u = set(nx.node_connected_component(T_temp, u_tree))

            for col_key in range(nt, m):
                u_chord, v_chord = edge_map[col_key]['u'], edge_map[col_key]['v']
                u_in_U = u_chord in comp_u
                v_in_U = v_chord in comp_u

                if u_in_U and not v_in_U:
                    Q_list[row_key][col_key] = 1
                elif not u_in_U and v_in_U:
                    Q_list[row_key][col_key] = -1

        return sp.Matrix(Q_list)

    # ==================== 参数矩阵 ====================

    def _build_parameter_matrices(self):
        '''构建 D_C (电容), L_plus (电感逆), D_J (约瑟夫森) 矩阵'''
        edge_map = self._edge_map
        mutual_dict = self._mutual_dict
        num_edges = len(edge_map)

        # 1. 电容矩阵 D_C
        Dc_list = [0] * num_edges
        for k, info in edge_map.items():
            if info['type'] == 'C':
                Dc_list[k] = info['value']
        D_C = sp.diag(*Dc_list)

        # 2. 约瑟夫森结矩阵 D_J
        Dj_list = [0] * num_edges
        for k, info in edge_map.items():
            if info['type'] == 'JJ':
                Dj_list[k] = info['value']
        D_J = sp.diag(*Dj_list)

        # 3. 电感矩阵 L_plus
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
                    row = local_to_global[i]
                    col = local_to_global[j]
                    L_plus[row, col] = Gamma_sub[i, j]

        return D_C, L_plus, D_J

    # ==================== 公开属性 ====================

    @property
    def nt(self):
        '''树枝数 (自由度)'''
        return self._nt

    @property
    def m(self):
        '''总支路数'''
        return self._m

    @property
    def num_loops(self):
        '''基本回路数 (= 连支数 = m - nt)'''
        return self._m - self._nt

    @property
    def edge_map(self):
        '''边映射表 {key: {u, v, type, value}}'''
        return self._edge_map

    @property
    def loops(self):
        '''基本回路信息列表'''
        return self._loops

    @property
    def F_C(self):
        '''基本割集矩阵 (带缓存)'''
        if self._F_C_cache is None:
            self._F_C_cache = self._compute_fundamental_cut_matrix()
        return self._F_C_cache

    # ==================== 外磁通编辑 ====================

    def set_external_flux(self, loop_index: int, flux):
        '''
        设置某个基本回路的外磁通

        参数:
            loop_index: 回路编号 (0-based, 对应 print_loops 中显示的回路编号)
            flux: 外磁通值, 可以是 sympy 表达式、字符串(自动转符号) 或数值
                  设为 0 表示该回路无外磁通
        '''
        if loop_index < 0 or loop_index >= self.num_loops:
            raise IndexError(f"回路编号 {loop_index} 超出范围 [0, {self.num_loops - 1}]")

        if isinstance(flux, str):
            flux = sp.symbols(flux)

        chord_key = self._loops[loop_index]['chord_key']
        self._external_fluxes[chord_key] = flux

    def get_external_flux(self, loop_index: int):
        '''获取某个基本回路当前设置的外磁通'''
        if loop_index < 0 or loop_index >= self.num_loops:
            raise IndexError(f"回路编号 {loop_index} 超出范围 [0, {self.num_loops - 1}]")
        chord_key = self._loops[loop_index]['chord_key']
        return self._external_fluxes.get(chord_key, 0)

    # ==================== 打印信息 ====================

    def print_edges(self):
        '''打印所有边 (树枝 + 连支) 的信息'''
        print(f"树枝数: {self._nt}, 连支数: {self._m - self._nt}, 总支路数: {self._m}")
        for k, v in self._edge_map.items():
            role = "树枝" if k < self._nt else "连支"
            print(f"  Key {k} [{role}]: {v['type']} ({v['u']}->{v['v']}) 值: {v['value']}")

    def print_loops(self):
        '''打印所有基本回路信息及当前外磁通设置'''
        print(f"电路共有 {self.num_loops} 个基本回路 (树枝数: {self.nt}, 总支路数: {self.m})")
        print()

        for i, loop in enumerate(self._loops):
            chord_key = loop['chord_key']
            chord = loop['chord_info']
            flux = self._external_fluxes.get(chord_key, 0)

            print(f"回路 {i}: (由连支 Key {chord_key} 闭合)")
            print(f"  连支: {chord['type']} ({chord['u']}->{chord['v']}), 值: {chord['value']}")
            # 回路路径: 节点序列 + 回到起点
            path_str = ' -> '.join(map(str, loop['path_nodes']))
            print(f"  回路路径: {path_str} -> {loop['path_nodes'][0]}")
            print(f"  包含树枝:")
            for bk in loop['tree_branch_keys']:
                info = self._edge_map[bk]
                print(f"    Key {bk}: {info['type']} ({info['u']}->{info['v']}), 值: {info['value']}")
            print(f"  外磁通: {flux}")
            print()

    # ==================== 哈密顿量计算 ====================

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
        F_C = self.F_C
        D_C, L_plus, D_J = self._build_parameter_matrices()

        nt = self._nt
        m = self._m

        # 1. 定义动态变量
        phi_t = sp.Matrix([sp.symbols(f'Phi_t_{i}') for i in range(nt)])
        q_t = sp.Matrix([sp.symbols(f'Q_t_{i}') for i in range(nt)])

        # 2. 构建外磁通向量 (使用用户设置的外磁通)
        phi_ext_list = []
        ext_flux_symbols = set()

        for k in range(m):
            if k < nt:
                phi_ext_list.append(0)
            else:
                flux = self._external_fluxes.get(k, 0)
                phi_ext_list.append(flux)
                # 收集外磁通中的自由符号
                if isinstance(flux, sp.Basic):
                    ext_flux_symbols.update(flux.free_symbols)

        phi_ext_vec = sp.Matrix(phi_ext_list)
        ext_flux_vars = sorted(ext_flux_symbols, key=lambda s: str(s))

        # 3. 全支路磁通: Phi_all = F_C^T * Phi_tree + Phi_ext
        Phi_vec = F_C.T * phi_t + phi_ext_vec

        # 4. 动能 (电场能量): H_kin = ½ q^T M^{-1} q
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



'''
if __name__ == "__main__":

    # ---- 第 1 步: 定义元件并创建电路 ----
    my_components = [
        Component(0, 2, 'JJ', 'EJ1'),
        Component(1, 2, 'JJ', 'EJ2'),
        Component(0, 1, 'L', 'L'),
        Component(0, 2, 'C', 'C1'),
        Component(1, 2, 'C', 'C2'),
    ]

    circuit = Circuit(my_components)

    # ---- 第 2 步: 查看电路信息和基本回路 ----
    print("=" * 50)
    print("电路边信息:")
    print("=" * 50)
    circuit.print_edges()

    print()
    print("=" * 50)
    print("基本回路信息:")
    print("=" * 50)
    circuit.print_loops()

    # ---- 第 3 步: 编辑外磁通 ----
    # 例如: 将回路 0 的外磁通设为符号 Phi_e, 回路 1 设为 0 (无外磁通)
    circuit.set_external_flux(0, 'Phi_e')
    circuit.set_external_flux(1, 0)

    print("修改外磁通后:")
    circuit.print_loops()

    # ---- 第 4 步: 计算哈密顿量 ----
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
'''
