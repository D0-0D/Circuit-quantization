import sympy as sp
import networkx as nx

class Component:
    '''
    定义元件, 用于校验输入
    '''
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
    '''
    定义互感系数
    输入: 节点1, 节点2, 互感系数
    '''
    def __init__(self, L1_id: int, L2_id: int, value):
        self.L1_id = L1_id
        self.L2_id = L2_id
        self.value = sp.symbols(value) if isinstance(value, str) else value



def build_circuit_graph(components:list[Component], mutuals: None|list[MutualInductance] = None):
    """
    根据元件列表构建电路多重图(MultiGraph)。
    参数:
        components: 包含 Component 对象的列表
        mutuals: 包含MutualInductance对象的列表
    返回: 
        T (MultiGraph): 最小生成树
        G (MultiGraph): 完整电路图
        edge_map (dict): 边索引到 (u, v, key) 的映射
        mutual_dict: 互感映射字典 {(支路A, 支路B): 互感值}
    """
    G = nx.MultiGraph()
    edge_map = {}
    # 1. 遍历列表并添加支路
    for i, comp in enumerate(components):
        G.add_edge(
            comp.u, 
            comp.v, 
            key=i, 
            type=comp.type, 
            value=comp.value, 
            weight=comp.weight
        )       
        edge_map[i] = (comp.u, comp.v, i)

    # 补足缺少编号的节点, 接地节点是最大编号节点, 所以接地节点一定在电路中
    # TODO: 添加校验电路是否是完整的一部分的代码(如果电路是孤立的两部分,接地也会出问题)
    max_node = max(G.nodes)
    G.add_nodes_from(range(max_node + 1))

    # 2. 生成含尽量少C的生成树
    T = nx.minimum_spanning_tree(G, weight='weight')

    # 3. 处理并校验互感
    mutual_dict = {}
    if mutuals:
        for m in mutuals:
            if m.L1_id not in edge_map or m.L2_id not in edge_map:
                raise ValueError(f"互感指定的支路 ID ({m.L1_id} 或 {m.L2_id}) 不存在。")

            u1, v1, k1 = edge_map[m.L1_id]
            u2, v2, k2 = edge_map[m.L2_id]

            type1 = G.edges[u1, v1, k1]['type']
            type2 = G.edges[u2, v2, k2]['type']

            # 校验：互感只能定义在电感(L)之间
            if type1 != 'L' or type2 != 'L':
                raise ValueError(
                    f"互感只能定义在电感(L)之间！"
                )

            # 字典里只存第一个节点编号小于第二个节点编号的互感值
            coupling_key = tuple(sorted((m.L1_id, m.L2_id)))
            mutual_dict[coupling_key] = m.value

    return T,G, edge_map, mutual_dict

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


print(T.edges)
print(G.edges)
print(edge_map)
print(mutual_dict[(1,2)])#读取时第一个参数小于第二个
