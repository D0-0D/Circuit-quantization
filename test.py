import sympy as sp
import networkx as nx
import numpy as np

# ==================== 基础类定义 (保持不变) ====================
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

# ==================== 图构建逻辑 (保持不变) ====================
def build_circuit_graph(components: list[Component], mutuals: None | list[MutualInductance] = None):
    # 1. 初始构建图
    G_orig = nx.MultiGraph()
    edge_map_orig = {}
    for i, comp in enumerate(components):
        G_orig.add_edge(comp.u, comp.v, key=i, type=comp.type, value=comp.value, weight=comp.weight)
        edge_map_orig[i] = {'u': comp.u, 'v': comp.v, 'type': comp.type, 'value': comp.value, 'orig_id': i}

    G_orig.add_nodes_from(range(max(G_orig.nodes) + 1 if G_orig.nodes else 0))

    # 2. 生成树 (最小权重优先：C -> JJ -> L)
    T_orig = nx.minimum_spanning_tree(G_orig, weight='weight')

    # 3. 区分树枝和连支，重新编号
    tree_edges = []
    chord_edges = []
    tree_orig_keys = set(k for _, _, k in T_orig.edges(keys=True))

    sorted_orig_keys = sorted(edge_map_orig.keys())
    
    for k in sorted_orig_keys:
        info = edge_map_orig[k]
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
        T_new.add_edge(u, v, key=current_new_key, type=info['type'], value=info['value'])
        new_edge_map[current_new_key] = {'u': u, 'v': v, 'type': info['type'], 'value': info['value']}
        old_to_new_key[info['orig_id']] = current_new_key
        current_new_key += 1
        
    # 处理连支
    G_new = T_new.copy() 
    for info in chord_edges:
        u, v = info['u'], info['v']
        if u > v: u, v = v, u
        G_new.add_edge(u, v, key=current_new_key, type=info['type'], value=info['value'])
        new_edge_map[current_new_key] = {'u': u, 'v': v, 'type': info['type'], 'value': info['value']}
        old_to_new_key[info['orig_id']] = current_new_key
        current_new_key += 1

    # 处理互感字典
    new_mutual_dict = {}
    if mutuals:
        for m in mutuals:
            nk1 = old_to_new_key[m.L1_id]
            nk2 = old_to_new_key[m.L2_id]
            new_mutual_dict[tuple(sorted((nk1, nk2)))] = m.value

    return T_new, G_new, new_edge_map, new_mutual_dict

# ==================== 矩阵生成核心函数 (保持不变) ====================

def fundamental_cut_matrix(T, G, edge_map):
    m = G.number_of_edges()
    nt = T.number_of_edges()
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
            u_chord, v_chord = edge_map[col_key]['u'], edge_map[col_key]['v']
            u_in_U = u_chord in comp_u
            v_in_U = v_chord in comp_u
            
            if u_in_U and not v_in_U:
                Q_list[row_key][col_key] = 1
            elif not u_in_U and v_in_U:
                Q_list[row_key][col_key] = -1

    return sp.Matrix(Q_list)

def build_parameter_matrices(edge_map, mutual_dict):
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

# ==================== 哈密顿量计算 (已修改：支持外磁通) ====================
def calculate_hamiltonian(F_C, D_C, L_plus, D_J):
    """
    计算最终的哈密顿量表达式，包含外磁通
    """
    # 获取维度信息
    nt = F_C.shape[0]  # 树枝数 (自由度)
    m = F_C.shape[1]   # 总支路数
    
    # 1. 定义动态变量
    # Phi_t: 树枝磁通 (广义坐标)
    # Q_t: 树枝电荷 (广义动量)
    phi_t = sp.Matrix([sp.symbols(f'Phi_t_{i}') for i in range(nt)])
    q_t = sp.Matrix([sp.symbols(f'Q_t_{i}') for i in range(nt)])
    
    # 2. 定义外磁通 (External Flux)
    # 对于每一个连支 (index >= nt)，它闭合了一个回路，分配一个外磁通 Phi_ext_k
    # 对于树枝 (index < nt)，外磁通偏移量为 0
    phi_ext_list = []
    ext_flux_vars = []
    
    for k in range(m):
        if k < nt:
            phi_ext_list.append(0)
        else:
            # 连支 k 对应的回路外磁通
            sym = sp.symbols(f'Phi_ext_{k}')
            phi_ext_list.append(sym)
            ext_flux_vars.append(sym)
            
    phi_ext_vec = sp.Matrix(phi_ext_list)
    
    # 3. 构建全支路磁通向量 Phi_vec
    # 关系式: Phi_all = F^T * Phi_tree + Phi_ext_offset
    # 原理：连支磁通 = (通过树枝闭合的路径磁通) + (回路穿过的外磁通)
    # 注意：这里的 F_C 是基本割集矩阵，其转置 F_C.T 将树枝磁通映射到全支路磁通
    Phi_vec = F_C.T * phi_t + phi_ext_vec
    
    # 4. 动能部分 (电场能量) - 仅与电荷有关，不受外磁通直接影响
    # M = F_C * D_C * F_C.T
    M = F_C * D_C * F_C.T
    try:
        M_inv = M.inv()
        H_kin = (sp.Rational(1, 2) * q_t.T * M_inv * q_t)[0]
    except:
        # 如果不可逆 (例如没有电容)，使用符号表示
        M_inv_sym = sp.Symbol("M^{-1}") 
        H_kin = (sp.Rational(1, 2) * q_t.T * M_inv_sym * q_t)[0]
        
    # 5. 势能部分 (磁场能量 - 线性部分)
    # H_L = 0.5 * Phi_all.T * L_plus * Phi_all
    # 这里直接使用全支路磁通向量计算，自然包含了外磁通产生的线性项和常数项
    H_mag_lin = (sp.Rational(1, 2) * Phi_vec.T * L_plus * Phi_vec)[0]
    
    # 6. 势能部分 (约瑟夫森结 - 非线性部分)
    # H_JJ = -Sum(Ej * cos(Phi_k))
    H_jj = 0
    # 遍历所有支路，寻找 JJ
    for k in range(m):
        coeff = D_J[k, k]
        if coeff != 0:
            phi_k = Phi_vec[k] # 这里取出的标量包含 Phi_t 和 Phi_ext
            H_jj -= coeff * sp.cos(phi_k)
            
    # 总哈密顿量
    H_total = H_kin + H_mag_lin + H_jj
    
    return H_total, M, phi_t, q_t, ext_flux_vars, Phi_vec

# ==================== 主程序 ====================
if __name__ == "__main__":
    # 示例电路：Flux Qubit 或者是 rf-SQUID 变种
    # 0, 1, 2 三个节点
    # 3条支路：
    # 0-1: JJ1
    # 1-2: JJ2
    # 2-0: L (电感)
    # 再加一个电容以便有质量矩阵
    
    my_components = [
        Component(0, 1, 'JJ', 'EJ1'),
        Component(1, 2, 'JJ', 'EJ2'),
        Component(0, 2, 'L', 'L_loop'), 
        Component(0, 1, 'C', 'C1'), # 并联在 JJ1 上的电容
        Component(1, 2, 'C', 'C2'), # 并联在 JJ2 上的电容
        Component(0, 2, 'C', 'C3')  # 并联在 L 上的电容
    ]
    my_mutuals = None

    print("--- 1. 构建图结构 ---")
    T, G, edge_map, mutual_dict = build_circuit_graph(my_components, my_mutuals)
    
    nt = T.number_of_edges()
    m = G.number_of_edges()
    print(f"树枝数 (自由度): {nt}, 总支路数: {m}")
    for k, v in edge_map.items():
        role = "树枝 (Tree)" if k < nt else "连支 (Chord)"
        print(f"  Key {k} [{role}]: {v['type']} ({v['u']}->{v['v']}) Val: {v['value']}")

    print("\n--- 2. 基本割集矩阵 F^(C) ---")
    F_C = fundamental_cut_matrix(T, G, edge_map)
    sp.pprint(F_C)

    print("\n--- 3. 参数矩阵 ---")
    D_C, L_plus, D_J = build_parameter_matrices(edge_map, mutual_dict)
    
    print("\n--- 4. 计算哈密顿量 (含外磁通) ---")
    H, M_mat, phi_vars, q_vars, ext_fluxes, Phi_all_vec = calculate_hamiltonian(F_C, D_C, L_plus, D_J)
    
    print(f"\n生成的广义坐标 (树枝):")
    sp.pprint(phi_vars)
    
    print(f"\n生成的外部磁通变量 (对应连支):")
    if ext_fluxes:
        sp.pprint(ext_fluxes)
    else:
        print("无 (无闭合回路)")

    print("\n>>> 全支路磁通表达式 (Phi_branch = F^T * Phi_t + Phi_ext) <<<")
    print("注意观察连支部分包含了 Phi_ext：")
    sp.pprint(Phi_all_vec)
    
    print("\n>>> 最终哈密顿量 H <<<")
    # 展开并简化一下以便查看
    H_expanded = sp.expand(H)
    sp.pprint(H_expanded)