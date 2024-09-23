import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D


G = nx.DiGraph()

# paper table
# # operation lo . 
# 1 S3 = A11 - A21 C11 
# 2 S1 = A21 + A22 A21 
# 3 T1 = B12 - B11 C22 
# 4 T3 = B22 - B12 B12 
# 5 P7 = $\times$(S3T3) C21 
# 6 S2 = S1 - A11 C12
# 7 P1 = $\times$(A11B11) C11 
# 8 T2 = B22 - T1 B11 
# 9 P5 = $\times$(S1T1) A11 
# 10 T4 = T2 - B21 C22 
# 11 P4 = $\times$(A22T4) A21 
# 12 S4 = A12 - S2 A22
# 13 P6 = $\times$(S2T2) C22
# 14 U2 = P1 + P6 C22
# 15 P2 = $\times$(A12B21) C12
# 16 U1 = P1 + P2 C11
# 17 U4 = U2 + P5 C12
# 18 U3 = U2 + P7 C22
# 19 U6 = U3 - P4 C21
# 20 U7 = U3 + P5 C22
# 21 P3 = $\times$(S4B22) A12
# 22 U5 = U4 + P3 C12

edge2label = {
    ('A11', 'S3'): r'$-$',
    ('A21', 'S3'): r'$-$',
    ('A21', 'S1'): r'$+$',
    ('A22', 'S1'): r'$+$',
    ('B12', 'T1'): r'$-$',
    ('B11', 'T1'): r'$-$',
    ('B22', 'T3'): r'$-$',
    ('B12', 'T3'): r'$-$',
    ('S3', 'P7'): r'$\times$',
    ('T3', 'P7'): r'$\times$',
    ('S1', 'S2'): r'$-$',
    ('A11', 'S2'): r'$-$',
    ('A11', 'P1'): r'$\times$',
    ('B11', 'P1'): r'$\times$',
    ('B22', 'T2'): r'$-$',
    ('T1', 'T2'): r'$-$',
    ('S1', 'P5'): r'$\times$',
    ('T1', 'P5'): r'$\times$',
    ('T2', 'T4'): r'$-$',
    ('B21', 'T4'): r'$-$',
    ('A22', 'P4'): r'$\times$',
    ('T4', 'P4'): r'$\times$',
    ('A12', 'S4'): r'$-$',
    ('S2', 'S4'): r'$-$',
    ('S2', 'P6'): r'$\times$',
    ('T2', 'P6'): r'$\times$',
    ('P1', 'U2'): r'$+$',
    ('P6', 'U2'): r'$+$',
    ('A12', 'P2'): r'$\times$',
    ('B21', 'P2'): r'$\times$',
    ('P1', 'U1'): r'$+$',
    ('P2', 'U1'): r'$+$',
    ('U2', 'U4'): r'$+$',
    ('P5', 'U4'): r'$+$',
    ('U2', 'U3'): r'$+$',
    ('P7', 'U3'): r'$+$',
    ('U3', 'U6'): r'$-$',
    ('P4', 'U6'): r'$-$',
    ('U3', 'U7'): r'$+$',
    ('P5', 'U7'): r'$+$',
    ('S4', 'P3'): r'$\times$',
    ('B22', 'P3'): r'$\times$',
    ('U4', 'U5'): r'$+$',
    ('P3', 'U5'): r'$+$'
}


label2color = {
    r'$-$': '#FF6347',
    r'$+$': '#32CD32',
    r'$\times$': '#1E90FF'
}
edge2color = {key: label2color[label] for key, label in edge2label.items()}

node2color = {}
for edge, label in edge2label.items():
    node2color[edge[1]] = label2color[label]


G.add_edges_from(edge2label.keys())

node2layer = dict.fromkeys(list(G.nodes), 0)
def dfs_max_depth(node, layer):
    if node2layer[node] <= layer:
        node2layer[node] = layer
    for child in G.neighbors(node):
        dfs_max_depth(child, layer + 1)

for root_node in ['A11', 'A12', 'A21', 'A22', 'B11', 'B12', 'B21', 'B22']:
    dfs_max_depth(root_node, 0)
    
layer2nodes = {}
for node, layer in node2layer.items():
    if layer not in layer2nodes:
        layer2nodes[layer] = []
    layer2nodes[layer].append(node)
layer2nodes = {key: list(sorted(value)) for key, value in layer2nodes.items()}

positions = {
    'A11': (0, 9), 'A12': (0, 6), 'A21': (0, 7), 'A22': (0, 8),
    'B11': (0, 5), 'B12': (0, 4), 'B21': (0, 3), 'B22': (0, 2),
    'P1': (1, 8), 'P2': (1, 7), 
    'S1': (1, 6), 'S3': (1, 9),
    'T1': (1, 4), 'T3': (1, 3),
    'P5': (2, 6), 'P7': (2, 5), 'S2': (2, 8), 'T2': (2, 3), 'U1': (2, 7),
    'P6': (3, 5), 'S4': (3, 4), 'T4': (3, 3),
    'P3': (4, 4), 'P4': (4, 3), 'U2': (4, 6),
    'U3': (5, 3), 'U4': (5, 6),
    'U5': (6, 6), 'U6': (6, 5), 'U7': (6, 4)
}


labels = {
    'A11': '$a$',
    'A12': '$b$',
    'A21': '$c$',
    'A22': '$d$',
    'B11': '$A$',
    'B12': '$C$',
    'B21': '$B$',
    'B22': '$D$',
    'S3': '$t$',
    'S1': '$c$',
    'T1': '$w$',
    'T3': '$C$',
    'P7': '$v$',
    'S2': '$u$',
    'P1': '$t$',
    'T2': '$A$',
    'P5': '$a$',
    'T4': '$w$',
    'P4': '$c$',
    'S4': '$d$',
    'P6': '$w$',
    'U2': '$w$',
    'P2': '$u$',
    'U1': '$t$',
    'U4': '$u$',
    'U3': '$w$',
    'U6': '$v$',
    'U7': '$w$',
    'P3': '$b$',
    'U5': '$u$'
}


positions = {key: [val[1], 6 - val[0]] for key, val in positions.items()}
node_colors = [node2color.get(node, 'lightblue') for node in G.nodes()]

nx.draw(
    G, 
    positions, 
    labels=labels, 
    with_labels=True, 
    node_size=600, 
    node_color=node_colors, 
    font_size=10, 
    font_weight='bold',
    arrows=True
)
# print(edge2color)
# nx.draw_networkx_edges(
#     G, 
#     positions,
#     edge_color=edge2color.values(),
#     arrows=True
# )

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Initial Node', markersize=10, markerfacecolor='lightblue'),
    Line2D([0], [0], marker='o', color='w', label=r'$-$ (Subtraction)', markersize=10, markerfacecolor='#FF6347'),
    Line2D([0], [0], marker='o', color='w', label=r'$+$ (Addition)', markersize=10, markerfacecolor='#32CD32'),
    Line2D([0], [0], marker='o', color='w', label=r'$\times$ (Multiplication)', markersize=10, markerfacecolor='#1E90FF')
]

plt.legend(handles=legend_elements, loc='lower right')
plt.title('Strassen-Winograd Algorithm Operations Graph')
# plt.show()
plt.savefig('./visualization_scripts/strassen_opeartions_graph.png', dpi=300)