import sympy as sp
import networkx as nx


def build_circuit_graph(components):
    """
    根据元件列表构建电路图。
    参数 components 格式: [(node1, node2, type, value), ...]
    返回: G (MultiGraph), edge_map (dict)
    """
    G = nx.MultiGraph()
    edge_map = {}
    
    # 遍历列表并添加边
    for i, (u, v, tp, val) in enumerate(components):
        tp = tp.upper()#兼容小写字母
        wt = 1 if tp == 'C' else 2
        #连支中尽量多C就是生成树尽量少C, 所以C的wt应该是1, 其它的是2
        G.add_edge(u, v, key=i, type=tp, value=val, weight=wt)
        edge_map[i] = (u, v, i)
    
    return G, edge_map

my_components = [
    (0, 1, 'JJ', 0.5),   
    (1, 2, 'L', 10e-3),  
    (2, 0, 'C', 0.5),   
    (1, 0, 'C', 2e-6),   
]

G, edge_map = build_circuit_graph(my_components)

T = nx.minimum_spanning_tree(G, weight='weight')

print(f"生成树支路 ID: {[key for u, v, key in T.edges(keys=True)]}")

