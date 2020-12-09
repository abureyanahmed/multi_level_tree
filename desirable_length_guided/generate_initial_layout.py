import networkx as nx
import math

#G = nx.full_rary_tree(9, 50)
#G = nx.full_rary_tree(3, 4)

def numberOfNodes(G, root, parent, count):  
  
    count[root] = 1
    for u in G.neighbors(root):  
          
        # Condition to omit reverse path  
        # path from children to parent  
        if u == parent:  
            continue
          
        # recursive call for DFS  
        numberOfNodes(G, u, root, count)  
          
        # update count[] value of parent  
        # using its children  
        count[root] += count[u]

def get_drawing_coordinates(G, cur_vertex, parent, start_angle, end_angle, cur_vertex_crd_x, cur_vertex_crd_y, crd_x, crd_y, count, label_to_id, edges_to_index, edge_distance):
    crd_x[label_to_id[cur_vertex]] = cur_vertex_crd_x
    crd_y[label_to_id[cur_vertex]] = cur_vertex_crd_y
    if parent==-1:
        n_chld = len(list(G.neighbors(cur_vertex)))
    else:
        n_chld = len(list(G.neighbors(cur_vertex)))-1
    chld_count = 0
    chld_frac = 0
    for u in G.neighbors(cur_vertex):
        if u == parent:
            continue
        chld_start_angle = chld_frac*(end_angle-start_angle) + start_angle
        chld_count = chld_count + count[u]
        chld_frac = chld_count/(count[cur_vertex]-1)
        chld_end_angle = chld_frac*(end_angle-start_angle) + start_angle
        chld_mid_angle = (chld_start_angle+chld_end_angle)/2
        edge_len = 50
        if (cur_vertex, u) in edges_to_index.keys():
          edge_len = edge_distance[edges_to_index[(cur_vertex, u)]]
        elif (u, cur_vertex) in edges_to_index.keys():
          edge_len = edge_distance[edges_to_index[(u, cur_vertex)]]
        else:
          print('edge not found in get_drawing_coordinates!')
          quit()
        chld_vertex_x = cur_vertex_crd_x + edge_len*math.cos(chld_mid_angle)
        chld_vertex_y = cur_vertex_crd_y + edge_len*math.sin(chld_mid_angle)
        get_drawing_coordinates(G, u, cur_vertex, chld_start_angle, chld_end_angle, chld_vertex_x, chld_vertex_y, crd_x, crd_y, count, label_to_id, edges_to_index, edge_distance)

'''
cnt = {}
src = list(G.nodes())[0]
numberOfNodes(G, src, -1, cnt)
crd_x = {}
crd_y = {}
crd_x[src] = 500
crd_y[src] = 500
get_drawing_coordinates(G, src, -1, 0, 2*math.pi, crd_x[src], crd_y[src], crd_x, crd_y, cnt)
print(crd_x)
print(crd_y)
'''



