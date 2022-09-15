import sys
import numpy as np
import networkx as nx

root = int(sys.argv[1])
#edge_file_name = "datasets/raw/Graph_8.txt.weighted.mtx"
#label_file_name = "datasets/raw/Graph_8.txt.labels"
#out_dir = "datasets/output/"
edge_file_name = sys.argv[2]
label_file_name = sys.argv[3]
out_dir = sys.argv[4]

my_edges = []
f = open(edge_file_name)
ln = f.readline()
ln = f.readline()
ln = ln.split()
n, m = int(ln[1]), int(ln[2])
for i in range(m):
 ln = f.readline()
 ln = ln.split()
 if len(ln)>2:
  u, v, w = int(ln[0])-1, int(ln[1])-1, float(ln[2])
 else:
  u, v, w = int(ln[0])-1, int(ln[1])-1, 100
 my_edges.append([u, v, w])
f.close()

label_to_index = {}
index_to_label = {}
f = open(label_file_name)
ln = f.readline()
cnt = 0
for i in range(n):
 label = ln[:len(ln)-1]
 label_to_index[label] = cnt
 index_to_label[cnt] = label
 cnt += 1
 ln = f.readline()
f.close()

def edges2graph(my_edges, i2k=None, label2i=None):
    nodes = set()
    edges = set()
    for i, edge in enumerate(my_edges):
            source, target, edge_len = edge
            nodes.update([i2k[source], i2k[target]])
            edges.add( (i2k[source], i2k[target], edge_len) )

    if label2i is None:
        label2i = {k:i for i,k in enumerate(nodes)}
        #i2k = list(range(len(nodes)))
        i2k = {label2i[k]:k for k in label2i.keys()}
    g = nx.Graph()

    nodes = [dict(id=label2i[k], label=k) for i,k in enumerate(nodes)]
    ids = [n['id'] for n in nodes]
    g.add_nodes_from( zip(ids, nodes) )

    #edges = [(i2k[label2i[e[0]]],i2k[label2i[e[1]]], e[2]) for e in edges]
    edges = [(label2i[e[0]],label2i[e[1]], e[2]) for e in edges]
    g.add_weighted_edges_from(edges)
    return g, i2k, label2i

def fan(nodes, origin=[0,0], radii=[], phaseCenter=0, phaseRange=np.pi, weights=[1,1], mode='random'):
  pos = {}
  phases = {}
  ranges = {}
  n = len(nodes)
  cos, sin = np.cos, np.sin
  weightTotal = sum(weights)
  weights = [w/weightTotal for w in weights]
  nr = sorted(zip(nodes, weights, radii), key=lambda x:x[1])

  nr2 = []
  for i in list(range(len(nr)))[::-1]:
    if i%2 == 0:
      nr2.append(nr[i])
    else:
      nr2.insert(0, nr[i])

  nodes, weights, radii = zip(*nr2)

  weightCumSum = [sum(weights[:i]) for i in range(len(weights)+1)]
  for i in range(n):
    angle_offset = (weightCumSum[i]+weightCumSum[i+1])/2 * phaseRange
    angle_i = phaseCenter - phaseRange/2 + angle_offset
    ri = radii[i]
    pos[nodes[i]] = [origin[0] + ri*cos(angle_i), origin[1] + ri*sin(angle_i)]
    phases[nodes[i]] = angle_i
    ranges[nodes[i]] = weights[i] * phaseRange * 0.9
  return pos, phases, ranges

def radial_layout(g, root=None, mode='center', origin=[0,0], phase0=0, range0=np.pi*2):
  g0 = g
  g = nx.bfs_tree(g, source=root)
  pos = {}
  phases = {}
  ranges = {}
  depth_from_root = nx.shortest_path_length(g, root)
  pos[root] = origin
  phases[root] = phase0
  ranges[root] = range0
  roots = [root, ]
  depth = 1
  while len(pos) < len(g.nodes):
    newRoots = []
    for root in roots:
      neighbors = [n for n in g.neighbors(root) if n not in pos]
      if len(neighbors) > 0:
        edge_lengths = [g0.edges[(root, n)]['weight'] for n in neighbors]
        subTreeSizes = [len(nx.bfs_tree(g, i).nodes) for i in neighbors]
        degrees = [g.degree[i] for i in neighbors]
        depths = [depth_from_root[i] for i in neighbors]
        weights = [z for x, y, z in zip(degrees, depths, subTreeSizes)]
        newRoots += neighbors
        newPos, newPhases, newRanges = fan(neighbors, mode=mode, origin=origin, radii=[depth for e in edge_lengths], phaseCenter=phases[root], phaseRange=ranges[root], weights=weights)
        pos.update(newPos)
        phases.update(newPhases)
        ranges.update(newRanges)
    roots = newRoots
    depth+=1
  return pos

subgraph, _, _ = edges2graph(my_edges, index_to_label, label_to_index)
g = subgraph
pos0 = radial_layout(g, root, mode='center')

filename = out_dir + "/radial.txt"
f = open(filename, 'w')
for i in range(n):
  f.write(str(10*pos0[i][0]) + "\t" + str(10*pos0[i][1]) + "\t" + str(i+1) + "\n")
f.close()

