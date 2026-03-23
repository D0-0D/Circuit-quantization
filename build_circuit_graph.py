'''边的正方向定义为从小号节点指向大号节点,无视输入的顺序'''




import sympy as sp
import networkx as nx


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
        return 2 if self.type == 'C' else 1

class MutualInductance:
    '''定义互感系数'''
    def __init__(self, L1_id: int, L2_id: int, value):
        self.L1_id = L1_id
        self.L2_id = L2_id
        self.value = sp.symbols(value) if isinstance(value, str) else value

def build_circuit_graph(components: list[Component], mutuals: None | list[MutualInductance] = None):
    """
    根据元件列表构建电路多重图(MultiGraph)，并对边的key重新编号，
    使得在最小生成树T中的边编号在前，不在其中的边编号在后。
    返回的 edge_map 和 mutual_dict 均使用新编号作为支路标识。

    参数:
        components: 包含 Component 对象的列表
        mutuals: 包含 MutualInductance 对象的列表

    返回:
        T (MultiGraph): 最小生成树（key已重新编号）
        G (MultiGraph): 完整电路图（key已重新编号）
        edge_map (dict): 新key -> 边信息字典，包含 'u', 'v', 'type', 'value'
        mutual_dict (dict): (新key1, 新key2) -> 互感值
    """
    # 1. 初始构建图（使用原始key）
    G_orig = nx.MultiGraph()
    edge_map_orig = {}          # 原始key -> (u, v, old_key, data)
    for i, comp in enumerate(components):
        G_orig.add_edge(
            comp.u,
            comp.v,
            key=i,
            type=comp.type,
            value=comp.value,
            weight=comp.weight
        )
        edge_map_orig[i] = (comp.u, comp.v, i, {
            'type': comp.type,
            'value': comp.value,
            'weight': comp.weight
        })

    # 补足节点（确保所有编号节点存在）
    max_node = max(G_orig.nodes) if G_orig.nodes else 0
    G_orig.add_nodes_from(range(max_node + 1))

    # 2. 生成含尽量少C的生成树（基于原始key）
    T_orig = nx.minimum_spanning_tree(G_orig, weight='weight')

    # 3. 处理互感（基于原始支路ID）
    mutual_dict_orig = {}       # (原始key1, 原始key2) -> 互感值
    if mutuals:
        for m in mutuals:
            if m.L1_id not in edge_map_orig or m.L2_id not in edge_map_orig:
                raise ValueError(f"互感指定的支路 ID ({m.L1_id} 或 {m.L2_id}) 不存在。")
            u1, v1, k1, d1 = edge_map_orig[m.L1_id]
            u2, v2, k2, d2 = edge_map_orig[m.L2_id]
            type1 = d1['type']
            type2 = d2['type']
            if type1 != 'L' or type2 != 'L':
                raise ValueError(f"互感只能定义在电感(L)之间！")
            # 统一存储为(小key, 大key)的形式，但原始key仍可能大小不一，暂时保持原样
            # 后续重编号后会重新排序
            mutual_dict_orig[(m.L1_id, m.L2_id)] = m.value

    # ------------------ 重新编号边的key ------------------
    # 收集所有原始边的信息（按原始key排序，保证确定性）
    edges_info = []  # 每个元素: (u, v, old_key, data)
    for old_key, (u, v, _, data) in edge_map_orig.items():
        edges_info.append((u, v, old_key, data))
    edges_info.sort(key=lambda x: x[2])   # 按原始key升序

    # 获取树边的原始key集合
    tree_keys = set()
    for u, v, key in T_orig.edges(keys=True):
        tree_keys.add(key)

    # 分离树边和非树边
    tree_edges = []
    non_tree_edges = []
    for u, v, old_key, data in edges_info:
        if old_key in tree_keys:
            tree_edges.append((u, v, old_key, data))
        else:
            non_tree_edges.append((u, v, old_key, data))

    # 创建新图，重新分配key
    G_new = nx.MultiGraph()
    G_new.add_nodes_from(G_orig.nodes)

    old_to_new_key = {}          # 原始key -> 新key
    new_edge_map = {}             # 新key -> 边信息字典

    new_key = 0

    # 添加树边
    for u, v, old_key, data in tree_edges:
        G_new.add_edge(u, v, key=new_key, **data)
        old_to_new_key[old_key] = new_key
        mi, ma = (u, v) if u < v else (v, u)
        new_edge_map[new_key] = {
            'u': mi,
            'v': ma,
            'type': data['type'],
            'value': data['value']
        }
        new_key += 1

    # 添加非树边
    for u, v, old_key, data in non_tree_edges:
        G_new.add_edge(u, v, key=new_key, **data)
        old_to_new_key[old_key] = new_key
        mi, ma = (u, v) if u < v else (v, u)
        new_edge_map[new_key] = {
            'u': mi,
            'v': ma,
            'type': data['type'],
            'value': data['value']
        }
        new_key += 1

    # 从树边信息重新构建最小生成树T_new（使用新key）
    T_new = nx.MultiGraph()
    T_new.add_nodes_from(G_new.nodes)
    for u, v, old_key, data in tree_edges:
        new_key = old_to_new_key[old_key]
        T_new.add_edge(u, v, key=new_key, **data)

    # 转换互感字典为使用新key
    new_mutual_dict = {}
    for (k1_orig, k2_orig), value in mutual_dict_orig.items():
        new_k1 = old_to_new_key[k1_orig]
        new_k2 = old_to_new_key[k2_orig]
        # 统一存储为(小key, 大key)形式
        coupling_key = tuple(sorted((new_k1, new_k2)))
        new_mutual_dict[coupling_key] = value

    return T_new, G_new, new_edge_map, new_mutual_dict



def fundamental_loop_matrix(T, G, edge_map):
    """
    构建基本回路矩阵，所有边的参考方向均为小节点→大节点。
    返回:
        B_f: 列表的列表，大小为 (连支数) × (总支路数)
        loop_fluxes: 每个回路的磁通量符号（sympy Symbol）
    """
    m = G.number_of_edges()
    nt = T.number_of_edges()
    # 构建树的无向图（用于路径查找）
    tree_graph = nx.Graph()
    # 树边信息字典：键为排序后的节点对，值为 (key, u, v) 其中 u<v
    tree_edge_info = {}
    for u, v, key in T.edges(keys=True):
        # 确保 u<v（已在构建 T 时保证）
        tree_graph.add_edge(u, v)
        a, b = (u, v) if u < v else (v, u)   # a<v
        tree_edge_info[(a, b)] = (key, u, v)

    n_loops = m - nt
    B_f = [[0] * m for _ in range(n_loops)]
    loop_fluxes = []

    loop_idx = 0
    for key in range(nt, m):  # 连支
        info = edge_map[key]
        u, v = info['u'], info['v']   # 已保证 u<v
        B_f[loop_idx][key] = 1        # 连支方向为正（小→大）

        # 在树中找 u 到 v 的路径
        path_nodes = nx.shortest_path(tree_graph, source=v, target=u)
        for i in range(len(path_nodes) - 1):
            a, b = path_nodes[i], path_nodes[i+1]
            # 获取树边key'
            key_a, key_b = (a, b) if a < b else (b, a)
            t_key, t_u, t_v = tree_edge_info[(key_a, key_b)]
            # 确定符号：如果路径方向 a->b 与边的参考方向 t_u->t_v 一致，则 +1，否则 -1
            if a < b:
                sign = 1
            else:
                sign = -1
            B_f[loop_idx][t_key] += sign

        # 用户输入磁通量符号
        flux_str = input(f"请输入回路 {loop_idx + nt} (对应连支 key={key}) 的磁通量符号: ").strip()
        if flux_str == "":
            flux_sym = sp.symbols(f'Phi_{loop_idx + nt}')
        else:
            flux_sym = sp.symbols(flux_str)
        loop_fluxes.append(flux_sym)

        loop_idx += 1

    return B_f, loop_fluxes


def fundamental_cut_matrix(T, G, edge_map):
    """
    构建基本割集矩阵，所有边的参考方向均为小节点→大节点。
    每个割集对应于一条树支，方向与该树支的方向一致（小→大）。

    返回:
        Q_f: 列表的列表，大小为 (树支数) × (总支路数)
    """
    m = G.number_of_edges()
    nt = T.number_of_edges()
    # 构建树的无向图（用于分割）
    tree_graph = nx.Graph()
    # 记录每条树支的信息：键为排序后的节点对，值为 key
    tree_edge_dict = {}
    for u, v, key in T.edges(keys=True):
        tree_graph.add_edge(u, v)
        a, b = (u, v) if u < v else (v, u)
        tree_edge_dict[(a, b)] = key

    # 初始化割集矩阵，行：树支，列：所有支路
    Q_f = [[0] * m for _ in range(nt)]

    # 对每条树支构造基本割集
    for row_key in range(nt):
        # 找到该树支对应的边信息
        # 由于树支 key 从 0 到 nt-1，可以直接从 T 中获取？
        # 但需要知道它的两个端点。由于 T 是 MultiGraph，我们可以遍历找到 key 对应的边。
        # 更可靠：根据 key 从 T 的边数据中获取
        u = v = None
        for uu, vv, kk in T.edges(keys=True):
            if kk == row_key:
                u, v = uu, vv
                break
        if u is None:
            raise ValueError(f"树中找不到 key={row_key} 的边")
        # 确保 u < v（T 中应该已保证）
        if u > v:
            u, v = v, u

        # 在树中移除该边，得到两个连通分量
        # 复制树图，移除该边
        tree_graph_copy = tree_graph.copy()
        tree_graph_copy.remove_edge(u, v)
        # 获取包含 u 的节点集合（通过 BFS/DFS）
        comp_u = set(nx.node_connected_component(tree_graph_copy, u))
        comp_v = set(tree_graph.nodes) - comp_u  # 另一个分量

        # 树支自身
        Q_f[row_key][row_key] = 1

        # 遍历所有连支（key >= nt）
        for col_key in range(nt, m):
            info = edge_map[col_key]
            a, b = info['u'], info['v']  # a < b
            # 判断连支的两个端点是否分别位于两个不同分量
            if (a in comp_u and b in comp_v):
                # 方向一致：连支的小节点在 u 侧，大节点在 v 侧
                Q_f[row_key][col_key] = 1
            elif (a in comp_v and b in comp_u):
                # 方向相反：连支的小节点在 v 侧，大节点在 u 侧
                Q_f[row_key][col_key] = -1
            # 否则不在割集中，保持 0

    return Q_f



class CircuitAnalyzer:
    def __init__(self, T, G, edge_map, mutual_dict, FLM, FCM, loop_fluxes):
        """
        初始化分析器。

        参数:
            T: 最小生成树 (MultiGraph)，边已重新编号，树支 key 0..nt-1
            G: 完整电路图 (MultiGraph)，边已重新编号，总支路数 m
            edge_map: dict, key -> {'u','v','type','value'}
            mutual_dict: dict, (key1,key2) -> 互感值 (符号或数字)
            FLM: 基本回路矩阵 (列表的列表) 大小为 (连支数) × (总支路数)
            FCM: 基本割集矩阵 (列表的列表) 大小为 (树支数) × (总支路数)
            loop_fluxes: 每个回路的磁通量符号列表，长度 = 连支数
        """
        self.T = T
        self.G = G
        self.edge_map = edge_map
        self.mutual_dict = mutual_dict
        self.FLM = sp.Matrix(FLM)   # 转换为 sympy 矩阵
        self.FCM = sp.Matrix(FCM)
        self.loop_fluxes = loop_fluxes

        self.m = G.number_of_edges()           # 总支路数
        self.nt = T.number_of_edges()          # 树支数
        self.nl = self.m - self.nt              # 连支数

        # 收集所有节点（用于后续节点通量，但这里我们以树支通量为坐标）
        self.nodes = list(G.nodes)

        # 符号定义：树支通量 Φ_t (长度为 nt)
        self.Phi_t = sp.Matrix(sp.symbols(f'Phi_t0:{self.nt}', real=True))
        # 树支通量的时间导数
        self.Phi_t_dot = sp.Matrix(sp.symbols(f'dPhi_t0:{self.nt}', real=True))
        # 广义动量 p (与树支通量共轭)
        self.p = sp.Matrix(sp.symbols(f'p0:{self.nt}', real=True))

        # 构建各种矩阵
        self._build_component_matrices()
        self._build_M_K_Q0_I0()
        self._build_lagrangian()
        self._build_hamiltonian()

    def _build_component_matrices(self):
        """构建电容矩阵 D_C、电感矩阵 L_full、临界电流向量 J_C"""
        # 1. 电容矩阵 D_C (对角)
        D_C_diag = []
        for key in range(self.m):
            tp = self.edge_map[key]['type']
            val = self.edge_map[key]['value']
            if tp == 'C':
                # 确保 val 是 sympy 表达式
                if isinstance(val, str):
                    val = sp.symbols(val)
                D_C_diag.append(val)
            else:
                D_C_diag.append(0)
        self.D_C = sp.diag(*D_C_diag)   # 大小为 m×m

        # 2. 临界电流向量 J_C (长度为 m)
        J_C_list = []
        for key in range(self.m):
            tp = self.edge_map[key]['type']
            val = self.edge_map[key]['value']
            if tp == 'JJ':
                # 假设 value 是约瑟夫森能量 E_J，临界电流 I_c = E_J (在 ℏ=2e=1 单位下)
                if isinstance(val, str):
                    val = sp.symbols(val)
                J_C_list.append(val)
            else:
                J_C_list.append(0)
        self.J_C = sp.Matrix(J_C_list)   # 列向量

        # 3. 电感矩阵 L_full (考虑互感)
        # 初始化 m×m 零矩阵
        L_full = sp.zeros(self.m, self.m)
        # 填充自感
        for key in range(self.m):
            tp = self.edge_map[key]['type']
            val = self.edge_map[key]['value']
            if tp == 'L':
                if isinstance(val, str):
                    val = sp.symbols(val)
                L_full[key, key] = val
        # 填充互感
        for (k1, k2), M_val in self.mutual_dict.items():
            # 互感值可能是字符串
            if isinstance(M_val, str):
                M_val = sp.symbols(M_val)
            # 对称填入
            L_full[k1, k2] = M_val
            L_full[k2, k1] = M_val
        self.L_full = L_full
        # 电感矩阵的伪逆 L⁺
        self.L_pinv = L_full.pinv()   # Moore-Penrose 伪逆

    def _build_M_K_Q0_I0(self):
        """构建质量矩阵 M、刚度矩阵 K、偏移电荷 Q0、偏移电流 I0 (假设外部磁通静态)"""
        # M = FCM * D_C * FCM^T
        self.M = self.FCM * self.D_C * self.FCM.T
        # K = FCM * L⁺ * FCM^T
        self.K = self.FCM * self.L_pinv * self.FCM.T

        # 外部磁通向量 ext_flux (长度为 m)
        # 对于树支 (key < nt): 0
        # 对于连支 (key >= nt): 对应回路的磁通量符号
        ext_flux_list = [0] * self.m
        for i in range(self.nl):
            key = self.nt + i
            ext_flux_list[key] = self.loop_fluxes[i]
        self.ext_flux = sp.Matrix(ext_flux_list)

        # Q0 = FCM * D_C * [0; d(ext_flux)/dt]
        # 假设静态 => d(ext_flux)/dt = 0，因此 Q0 = 0
        self.Q0 = sp.zeros(self.nt, 1)

        # I0 = FCM * L⁺ * [0; ext_flux]
        # 构建一个长度为 m 的向量 [0; ext_flux]，其中前 nt 个为0，后 nl 个为 ext_flux 对应分量
        ext_flux_vector = sp.Matrix.vstack(
            sp.zeros(self.nt, 1),
            sp.Matrix([self.ext_flux[self.nt + i] for i in range(self.nl)])
        )
        self.I0 = self.FCM * self.L_pinv * ext_flux_vector

    def _build_lagrangian(self):
        """构建拉格朗日量 L = T - U (公式 B29)"""
        # 动能项 T = 1/2 * Φ_t_dot^T * M * Φ_t_dot + Q0^T * Φ_t_dot
        T_kin = sp.Rational(1,2) * self.Phi_t_dot.T * self.M * self.Phi_t_dot + self.Q0.T * self.Phi_t_dot
        # 势能项 U = 1/2 * Φ_t^T * K * Φ_t + I0^T * Φ_t - J_C^T * cos( (FCM^T * Φ_t) + ext_flux )
        # 注意 cos 是逐元素余弦，作用于向量
        arg = self.FCM.T * self.Phi_t + self.ext_flux   # 长度为 m 的向量
        cos_term = sp.Matrix([sp.cos(arg[i]) for i in range(self.m)])
        U_pot = (sp.Rational(1,2) * self.Phi_t.T * self.K * self.Phi_t +
                 self.I0.T * self.Phi_t -
                 self.J_C.T * cos_term)
        # 拉格朗日量是标量，但表达式中有矩阵乘积，需要提取 [0,0] 元素（因为结果是 1×1 矩阵）
        self.L = (T_kin - U_pot)[0,0]   # 转换为标量

    def _build_hamiltonian(self):
        """勒让德变换得到哈密顿量"""
        # 广义动量 p = ∂L/∂Φ_t_dot = M * Φ_t_dot + Q0
        # 因此 Φ_t_dot = M^{-1} (p - Q0)
        # 需要 M 可逆
        M_inv = self.M.inv()   # 假设 M 可逆
        # 构建表达式
        p_minus_Q0 = self.p - self.Q0
        # 动能部分 (p - Q0)^T * M^{-1} * (p - Q0) / 2
        kinetic = sp.Rational(1,2) * p_minus_Q0.T * M_inv * p_minus_Q0
        # 势能部分与拉格朗日量中的 U 相同，但 Φ_t_dot 替换为用 p 表示，而 p 是独立变量
        # 哈密顿量 H = p^T * Φ_t_dot - L
        # 但更简单直接代入公式 H = kinetic + potential
        potential = (sp.Rational(1,2) * self.Phi_t.T * self.K * self.Phi_t +
                     self.I0.T * self.Phi_t -
                     self.J_C.T * sp.Matrix([sp.cos((self.FCM.T * self.Phi_t + self.ext_flux)[i]) for i in range(self.m)]))
        self.H_classical = (kinetic + potential)[0,0]

    def quantize(self):
        """量子化：将 Φ_t 和 p 替换为算符，满足对易关系 [Φ_ti, pj] = i ℏ δ_ij
           在自然单位 ℏ=1 下，[Φ_ti, pj] = i δ_ij
           返回量子哈密顿量符号表达式（仍用符号表示算符，但引入对易关系）
        """
        # 此处我们仅返回符号表达式，不进行算符代数运算
        # 可以用符号替换，但不对易关系不自动处理
        # 实际模拟时通常将坐标和动量用产生湮灭算符展开
        # 这里我们仅做形式上的替换
        # 定义算符符号（也可以用 sympy 的 NonCommutative 符号）
        Phi_hat = sp.Matrix([sp.Symbol(f'\\hat{{\\Phi}}_{i}', commutative=False) for i in range(self.nt)])
        p_hat = sp.Matrix([sp.Symbol(f'\\hat{{p}}_{i}', commutative=False) for i in range(self.nt)])
        # 将经典表达式中的 Φ_t 和 p 替换为算符符号
        H_quantum = self.H_classical.subs({self.Phi_t[i]: Phi_hat[i] for i in range(self.nt)})
        H_quantum = H_quantum.subs({self.p[i]: p_hat[i] for i in range(self.nt)})
        return H_quantum

    def linearize(self):
        """线性化：在势能最小值附近展开，得到简正模式
           返回：简正频率列表，以及变换矩阵
        """
        # 先找到势能的极小点
        # 势能部分为 U = 0.5 * Φ_t^T * K * Φ_t + I0^T * Φ_t - J_C^T * cos( (FCM^T * Φ_t) + ext_flux )
        # 对于线性化，我们通常将余弦展开到二阶，然后找到平衡点
        # 但更系统的方法是：将 Φ_t 用平衡点加上小量，展开到二阶
        # 这里简化处理：假设平衡点为 Φ_t = 0，然后对余弦展开到二阶
        # 即 cos(x) ≈ 1 - x^2/2
        # 那么 U ≈ 0.5 Φ_t^T K Φ_t + I0^T Φ_t - J_C^T [1 - 0.5*(FCM^T Φ_t + ext_flux)^2]
        # 常数项忽略，线性项为 I0^T Φ_t + J_C^T * (FCM^T Φ_t + ext_flux) 的线性部分？需要小心处理
        # 这里不展开详细计算，仅示意
        # 实际应用中，通常先做变量平移消除线性项
        pass

    def truncate_to_qubit(self):
        """截断到二能级（假设已经通过简正模式分解为多个模，每个模具有小非线性）
           此处仅示意
        """
        pass


# ------------------ 测试 ------------------
if __name__ == "__main__":
    # 沿用之前的元件和互感
    my_components = [
        Component(0, 1, 'JJ', 'J'),
        Component(1, 2, 'L', 'L1'),
        Component(2, 0, 'L', 'L2'),
        Component(1, 0, 'C', 'C_shunt'),   # 电容
    ]

    my_mutuals = [
        MutualInductance(L1_id=1, L2_id=2, value='M_12')
    ]

    T, G, edge_map, mutual_dict = build_circuit_graph(my_components, my_mutuals)

    print("树T中的边（新key）：", list(T.edges(keys=True)))
    print("全图G中的边（新key）：", list(G.edges(keys=True)))
    print("新的edge_map（新key -> 边信息）：")
    for k, info in edge_map.items():
        print(f"  {k}: {info}")
    print("新的mutual_dict（新key对 -> 互感值）：", mutual_dict)

    # 构建基本回路矩阵
    print("\n开始构建基本回路矩阵...")
    FLM, loop_fluxes = fundamental_loop_matrix(T, G, edge_map)

    print("\n基本回路矩阵 FLM：")
    for row in FLM:
        print(row)
    print("回路磁通量符号：", loop_fluxes)

    # 构建基本割集矩阵
    print("\n开始构建基本割集矩阵...")
    FCM = fundamental_cut_matrix(T, G, edge_map)

    print("\n基本割集矩阵 FCM：")
    for row in FCM:
        print(row)

    print(f"\n树T的边数: {T.number_of_edges()}")
    print(f"全图G的边数: {G.number_of_edges()}")

    # 创建分析器
    analyzer = CircuitAnalyzer(T, G, edge_map, mutual_dict, FLM, FCM, loop_fluxes)

    print("\n质量矩阵 M:")
    sp.pprint(analyzer.M)

    print("\n刚度矩阵 K:")
    sp.pprint(analyzer.K)

    print("\n偏移电流 I0:")
    sp.pprint(analyzer.I0)

    print("\n拉格朗日量 L:")
    sp.pprint(analyzer.L)

    print("\n经典哈密顿量 H:")
    sp.pprint(analyzer.H_classical)

    # 量子化
    H_quantum = analyzer.quantize()
    print("\n量子哈密顿量 (形式):")
    sp.pprint(H_quantum)


# ------------------ 测试 ------------------
if __name__ == "__main__":
    my_components = [
        Component(0, 1, 'JJ', 'J'),
        Component(1, 2, 'L', 'L1'),
        Component(2, 0, 'L', 'L2'),
        Component(1, 0, 'C', 2e-6),
    ]

    my_mutuals = [
        MutualInductance(L1_id=1, L2_id=2, value='M_12')
    ]

    T, G, edge_map, mutual_dict = build_circuit_graph(my_components, my_mutuals)

    print("树T中的边（新key）：", list(T.edges(keys=True)))
    print("全图G中的边（新key）：", list(G.edges(keys=True)))
    print("新的edge_map（新key -> 边信息）：")
    for k, info in edge_map.items():
        print(f"  {k}: {info}")
    print("新的mutual_dict（新key对 -> 互感值）：", mutual_dict)

    # 输出边数统计
    et = T.number_of_edges()
    eg = G.number_of_edges()
    print(f"\n树T的边数: {et}")
    print(f"全图G的边数: {eg}")
    
    
    # 构建基本回路矩阵
    print("\n开始构建基本回路矩阵...")
    FLM, loop_fluxes = fundamental_loop_matrix(T, G, edge_map)

    print("\n基本回路矩阵 FLM：")
    for row in FLM:
        print(row)
    print("回路磁通量符号：", loop_fluxes)
    
    
    print("\n开始构建基本割集矩阵...")
    FCM = fundamental_cut_matrix(T, G, edge_map)

    print("\n基本割集矩阵 FCM：")
    for row in FCM:
        print(row)
    

    print(f"\n树T的边数: {T.number_of_edges()}")
    print(f"全图G的边数: {G.number_of_edges()}")
    
    
analyzer = CircuitAnalyzer(T, G, edge_map, mutual_dict, FLM, FCM, loop_fluxes)


print("质量矩阵 M:")
analyzer.M.pprint()
print("\n刚度矩阵 K:")
analyzer.K.pprint()
print("\n拉格朗日量 L:")
analyzer.L.pprint()
print("\n经典哈密顿量 H:")
analyzer.H_classical.pprint()
print("\n量子哈密顿量 (形式):")
analyzer.quantize().pprint()