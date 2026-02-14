'''边的正方向定义为从小号节点指向大号节点'''




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
        flux_str = input(f"请输入回路 {loop_idx} (对应连支 key={key}) 的磁通量符号: ").strip()
        if flux_str == "":
            flux_sym = sp.symbols(f'Phi_{loop_idx}')
        else:
            flux_sym = sp.symbols(flux_str)
        loop_fluxes.append(flux_sym)

        loop_idx += 1

    return B_f, loop_fluxes



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
    print(f"\n树T的边数: {T.number_of_edges()}")
    print(f"全图G的边数: {G.number_of_edges()}")
    
    
    # 构建基本回路矩阵
    print("\n开始构建基本回路矩阵...")
    B_f, loop_fluxes = fundamental_loop_matrix(T, G, edge_map)

    print("\n基本回路矩阵 B_f：")
    for row in B_f:
        print(row)
    print("回路磁通量符号：", loop_fluxes)

    print(f"\n树T的边数: {T.number_of_edges()}")
    print(f"全图G的边数: {G.number_of_edges()}")