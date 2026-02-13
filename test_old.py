import sympy as sp
import networkx as nx

G = nx.MultiGraph() #Q:应该不要用有向图吧
edge_map = {}


# automatically restore the id for future reference
def add_edge_with_id(G, u, v, edge_id, **attr):
    G.add_edge(u, v, key=edge_id, **attr)
    edge_map[edge_id] = (u, v, edge_id)
    return edge_id


n = int(input()) # number of nodes
G.add_nodes_from(range(0, n)) # node "n - 1" is by default the ground node

m = int(input()) # number of edges
for i in range(0, m):
    n1 = int(input())
    n2 = int(input())
    tp = str(input()) # L, C, JJ
    val = float(input())
    if(tp == C):
        wt = 2
    else:
        wt = 1
    add_edge_with_id(G, n1, n2, i, type = tp, value = val, weight = wt) # directed from n1 to n2


# construct the spanning tree
T = nx.minimum_spanning_tree(G, weight = 'weight')