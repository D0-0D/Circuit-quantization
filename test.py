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

# ==================== 图构建逻辑 (微调以适应计算) ====================
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

    # 保持原始定义的顺序稳定性，先按原始ID排序
    sorted_orig_keys = sorted(edge_map_orig.keys())
    
    for k in sorted_orig_keys:
        info = edge_map_orig[k]
        if k in tree_orig_keys:
            tree_edges.append(info)
        else:
            chord_edges.append(info)

    # 重新构建映射
    # 新编号：0 ~ nt-1 为树枝，nt ~ m-1 为连支
    new_edge_map = {}
    old_to_new_key = {}
    
    current_new_key = 0
    
    # 处理树枝
    T_new = nx.MultiGraph()
    T_new.add_nodes_from(G_orig.nodes)
    for info in tree_edges:
        u, v = info['u'], info['v']
        # 强制方向 u < v 方便矩阵处理，但记录原始方向并不影响无向图的生成树
        if u > v: u, v = v, u
        T_new.add_edge(u, v, key=current_new_key, type=info['type'], value=info['value'])
        new_edge_map[current_new_key] = {'u': u, 'v': v, 'type': info['type'], 'value': info['value']}
        old_to_new_key[info['orig_id']] = current_new_key
        current_new_key += 1
        
    # 处理连支
    G_new = T_new.copy() # 包含树枝
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

# ==================== 矩阵生成核心函数 ====================

def fundamental_cut_matrix(T, G, edge_map):
    """
    构建基本割集矩阵 F^(C) (SymPy Matrix).
    行：树枝 (0 ~ nt-1)
    列：所有支路 (0 ~ m-1)
    """
    m = G.number_of_edges()
    nt = T.number_of_edges()
    
    # 初始化矩阵 (nt x m)
    # Q_f[i][j] 表示第i个割集是否包含第j条支路
    Q_list = [[0] * m for _ in range(nt)]

    # 构建树的辅助图用于连通性分析
    tree_graph = nx.Graph()
    for u, v, k in T.edges(keys=True):
        tree_graph.add_edge(u, v, key=k)

    for row_key in range(nt):
        # 1. 树枝本身在割集中系数为 +1
        Q_list[row_key][row_key] = 1
        
        # 获取当前树枝的端点
        u_tree, v_tree = edge_map[row_key]['u'], edge_map[row_key]['v']
        
        # 2. 移除该树枝，产生两个连通分量
        T_temp = tree_graph.copy()
        T_temp.remove_edge(u_tree, v_tree)
        
        # 确定包含 u_tree 的分量节点集
        comp_u = set(nx.node_connected_component(T_temp, u_tree))
        # comp_v = set(G.nodes) - comp_u
        
        # 3. 检查所有连支 (key >= nt)
        for col_key in range(nt, m):
            u_chord, v_chord = edge_map[col_key]['u'], edge_map[col_key]['v']
            
            # 判断连支是否跨越两个分量
            u_in_U = u_chord in comp_u
            v_in_U = v_chord in comp_u
            
            if u_in_U and not v_in_U:
                # 连支 u 端在 U 集，v 端在 V 集。
                # 树枝方向定义为 U -> V (因为 u_tree < v_tree 且 u_tree 在 U 中)
                # 连支方向定义为 u_chord -> v_chord
                # 方向一致
                Q_list[row_key][col_key] = 1
            elif not u_in_U and v_in_U:
                # 连支 u 端在 V 集，v 端在 U 集。
                # 方向相反
                Q_list[row_key][col_key] = -1
            # 否则连支在同一侧，不穿过割集，为 0

    return sp.Matrix(Q_list)

def build_parameter_matrices(edge_map, mutual_dict):
    """
    构建 D_C, L_plus, D_J 矩阵
    """
    num_edges = len(edge_map)
    
    # 1. 电容矩阵 D_C (对角)
    Dc_list = [0] * num_edges
    for k, info in edge_map.items():
        if info['type'] == 'C':
            Dc_list[k] = info['value']
    D_C = sp.diag(*Dc_list)
    
    # 2. 约瑟夫森结矩阵 D_J (对角)
    Dj_list = [0] * num_edges
    for k, info in edge_map.items():
        if info['type'] == 'JJ':
            Dj_list[k] = info['value'] # 这里存的是 Ic 或 Ej
    D_J = sp.diag(*Dj_list)
    
    # 3. 电感矩阵处理 (L -> L_inverse -> L_plus)
    # 首先识别所有 L 类型的支路索引
    L_indices = [k for k, info in edge_map.items() if info['type'] == 'L']
    
    if not L_indices:
        L_plus = sp.zeros(num_edges, num_edges)
    else:
        # 构建电感子矩阵
        n_L = len(L_indices)
        L_sub = sp.zeros(n_L, n_L)
        
        # 建立 局部索引 <-> 全局索引 的映射
        local_to_global = {i: k for i, k in enumerate(L_indices)}
        
        for i in range(n_L):
            k_i = local_to_global[i]
            # 自感
            L_sub[i, i] = edge_map[k_i]['value']
            # 互感
            for j in range(i + 1, n_L):
                k_j = local_to_global[j]
                key_tuple = tuple(sorted((k_i, k_j)))
                if key_tuple in mutual_dict:
                    val = mutual_dict[key_tuple]
                    L_sub[i, j] = val
                    L_sub[j, i] = val
        
        # 求逆得到电感系数矩阵 (Gamma)
        try:
            Gamma_sub = L_sub.inv()
        except:
            # 如果符号计算无法求逆，或者奇异，抛出错误
            print("警告：电感矩阵奇异，可能包含无电感回路或未定义的互感。")
            Gamma_sub = sp.zeros(n_L, n_L)

        # 映射回全尺寸矩阵 L^+
        L_plus = sp.zeros(num_edges, num_edges)
        for i in range(n_L):
            for j in range(n_L):
                row = local_to_global[i]
                col = local_to_global[j]
                L_plus[row, col] = Gamma_sub[i, j]
                
    return D_C, L_plus, D_J

# ==================== 哈密顿量计算 ====================
def calculate_hamiltonian(F_C, D_C, L_plus, D_J):
    """
    计算最终的哈密顿量表达式
    """
    # 1. 计算质量矩阵 M 和 势能矩阵 K
    M = F_C * D_C * F_C.T
    K = F_C * L_plus * F_C.T
    
    nt = M.shape[0] # 自由度数量 (树枝数)
    
    # 2. 定义动态变量
    # Phi_t: 树枝磁通 (广义坐标)
    # Q_t: 树枝电荷 (广义动量)
    phi_t = sp.Matrix([sp.symbols(f'Phi_t_{i}') for i in range(nt)])
    q_t = sp.Matrix([sp.symbols(f'Q_t_{i}') for i in range(nt)])
    
    # 3. 动能部分 (电场能量)
    # H_C = 0.5 * Q_t.T * M^(-1) * Q_t
    try:
        M_inv = M.inv()
        # [0] 用于提取 1x1 矩阵中的标量元素
        H_kin = (sp.Rational(1, 2) * q_t.T * M_inv * q_t)[0]
    except:
        print("警告：电容矩阵 M 不可逆，使用符号 M_inv 代替。")
        # 创建一个符号代表逆矩阵，或者提示用户需要降维
        # 注意：实际上如果是奇异矩阵，意味着存在约束，不能简单求逆。
        # 这里为了让代码跑通，我们把 M_inv 视为一个标量系数或者保持矩阵形式但取其元素
        # 修正点：即使 M_inv 是符号，q_t.T * Scalar * q_t 依然会返回 1x1 矩阵
        M_inv_sym = sp.Symbol("M^{-1}") 
        # 强制取 [0] 转换为标量
        H_kin = (sp.Rational(1, 2) * q_t.T * M_inv_sym * q_t)[0]
        
    # 4. 势能部分 (磁场能量 - 线性部分)
    # H_L = 0.5 * Phi_t.T * K * Phi_t
    # 同样确保取出标量
    H_mag_lin = (sp.Rational(1, 2) * phi_t.T * K * phi_t)[0]
    
    # 5. 势能部分 (约瑟夫森结 - 非线性部分)
    # 变换回全支路磁通向量: Phi = F_C.T * Phi_t
    Phi_vec = F_C.T * phi_t 
    
    H_jj = 0
    m = D_J.shape[0]
    for k in range(m):
        coeff = D_J[k, k]
        if coeff != 0:
            phi_k = Phi_vec[k] # 这里取出的已经是标量
            # 约瑟夫森势能: -E_J * cos(phi)
            H_jj -= coeff * sp.cos(phi_k)
            
    # 总哈密顿量 (现在所有项都是标量了)
    H_total = H_kin + H_mag_lin + H_jj
    
    return H_total, M, K, phi_t, q_t


# ==================== 主程序 ====================
if __name__ == "__main__":
    # 示例电路：rf-SQUID 或者 简单的 LCJ 回路
    # 节点 0, 1, 2
    # 0-1: JJ
    # 1-2: L1
    # 2-0: L2 (与 L1 互感)
    # 1-0: C (分流电容)
    
    my_components = [
        Component(0, 1, 'JJ', 'E_J'),  # 支路 0
        Component(1, 2, 'L', 'L1'),    # 支路 1
        Component(2, 0, 'L', 'L2'),    # 支路 2
        Component(1, 0, 'C', 'C1'),    # 支路 3
        Component(1, 4, 'C', 'C2')
    ]

    my_mutuals = [
        MutualInductance(L1_id=1, L2_id=2, value='M_12')
    ]

    print("--- 1. 构建图结构 ---")
    T, G, edge_map, mutual_dict = build_circuit_graph(my_components, my_mutuals)
    
    nt = T.number_of_edges()
    m = G.number_of_edges()
    print(f"树枝数: {nt}, 总支路数: {m}")
    for k, v in edge_map.items():
        role = "树枝" if k < nt else "连支"
        print(f"  Key {k} ({role}): {v['type']} ({v['u']}->{v['v']}), Val: {v['value']}")

    print("\n--- 2. 构建基本割集矩阵 F^(C) ---")
    F_C = fundamental_cut_matrix(T, G, edge_map)
    sp.pprint(F_C)

    print("\n--- 3. 构建参数矩阵 ---")
    D_C, L_plus, D_J = build_parameter_matrices(edge_map, mutual_dict)
    
    print("电感逆矩阵 L^+ (非零部分示意):")
    sp.pprint(L_plus)

    print("\n--- 4. 计算哈密顿量 ---")
    H, M_mat, K_mat, phi_vars, q_vars = calculate_hamiltonian(F_C, D_C, L_plus, D_J)
    
    print("\n广义坐标 (树枝磁通):")
    sp.pprint(phi_vars)
    
    print("\n电容矩阵 M (M = F C F^T):")
    sp.pprint(M_mat)
    
    print("\n电感矩阵 K (K = F L^+ F^T):")
    sp.pprint(K_mat)
    
    print("\n>>> 最终哈密顿量 H(Q_t, Phi_t) <<<")
    # 简化输出以便阅读
    sp.pprint(H)