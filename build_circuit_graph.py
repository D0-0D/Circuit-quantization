import sympy as sp
import networkx as nx

class Component:
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
    

def build_circuit_graph(components:list[Component]):
    """
    根据元件列表构建电路多重图(MultiGraph)。
    参数:
        components: 包含 Component 对象的列表
    返回: 
        T (MultiGraph): 最小生成树
        G (MultiGraph): 完整电路图
        edge_map (dict): 边索引到 (u, v, key) 的映射
    """
    G = nx.MultiGraph()
    edge_map = {}

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
    T = nx.minimum_spanning_tree(G, weight='weight')

    return T,G, edge_map

my_components = [
    Component(0, 1, 'JJ', 'J'),   
    Component(1, 2, 'L', 10e-3),  
    Component(2, 0, 'C', 0.5),   
    Component(1, 0, 'C', 2e-6),   
]

T, G, edge_map = build_circuit_graph(my_components)


print(f"生成树支路 ID: {[key for u, v, key in T.edges(keys=True)]}")

