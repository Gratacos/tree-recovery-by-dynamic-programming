import queue

from networkx.algorithms.distance_measures import center
from numpy.lib.utils import source
import partitiongraph as pg
import networkx as nx
import numpy as np
from queue import PriorityQueue, Queue
import scipy.linalg as la
import scipy.sparse.linalg as sp_la
import scipy
import time
mode = 1
SPECTRAL_TIME = 0
FIEDLER_TIME = 0
CUT_TIME = 0
STORE_ADJ = False
ADJ_TIMES = {'lap':0, 'shuffle':0, 'shuffle2':0, 'loops':0, 'adj':0}

# Example use:
#   Compute merge sequence bottom up:
#       merge_seq = priority_merge(G.copy(), ctype="outdegree")
#       T = iter_trigraph(G, merge_seq)
#       Mt = merge_tree(G, T, T.graph['root edge'])
#   Compute merge sequence top-bottom:
#       T = balanced_cut_trigraph(G, normalized=True)
#       Mt = merge_tree(G, T, T.graph['root edge'], rename_nodes=True)
#       Mt, total_cost = compute_costs_l(G, Mt)

meta = {
    'unspecified root':False,
    'bell':True,
    'cost function':'idempotent'
}


def compute_merge_sequence(G):
    # Obtain merge sequence
    FIEDLER_TIME = 0
    CUT_TIME = 0
    mode = 2 # eigsh with sigma 0.0 ... 01
    # mseq.mode = 4 # Fielder
    method = 2   # Graph Cut

    new_g, pre_mseq, pre_merge_map = pre_merge(G)

    if new_g.number_of_nodes() != 1:
        Mt = topdown_merge_tree(new_g, normalized=False, method=method)

        # Save balanced progress
        r = [n for n in Mt if nx.degree(Mt,n) == 2][0]
        ms = merge_sequence(Mt, new_g.number_of_nodes(), r)
        if pre_mseq != []:
            seq = merge_sequences(pre_mseq, ms.copy(), pre_merge_map)

        Mt, cost_change = obtain_root(Mt, r, new_g)

        # save root progress
        r = [n for n in Mt if nx.degree(Mt,n) == 2][0]
        ms = merge_sequence(Mt, new_g.number_of_nodes(), r)

        if pre_mseq != []:
            seq = merge_sequences(pre_mseq, ms, pre_merge_map)

        r = [n for n in Mt if len(list(Mt.predecessors(n)))==0][0]
        Mt, total_cost = compute_costs_l(new_g, Mt, r=r, N=None)
        ccl = compute_cost_change_list(Mt)
        pivot_global(Mt, ccl)

        seq = merge_sequence(Mt, new_g.number_of_nodes(), r)
        if pre_mseq != []:
            seq = merge_sequences(pre_mseq, seq, pre_merge_map)
    else:
        seq = pre_mseq

    return seq



def merge(G, u, v, next_name=None, distances=False, names=False):
    """ Physical merge of 2 nodes """
    if next_name is None:
        next_name = max(list(G.nodes())) + 1
    udeg = nx.degree(G, u)
    vdeg = nx.degree(G, v)

    #print('Merging:', [G.edges[e]['name'] for e in G.edges(u, keys=True)], '&', [G.edges[e]['name'] for e in G.edges(v, keys=True)])

    # Merge (u, v) into u, remove all in_edges, and put all out edges into next_name
    in_edges = [e for e in G.edges(u, keys=True) if (e[0] == u and e[1] == v) or (e[0] == v and e[1] == u)]

    u_out_edges = [(e[0], e[1], e[2], G.edges[e]['name']) for e in G.edges(u, keys=True) if e[0] != v and e[1] != v]
    v_out_edges = [(e[0], e[1], e[2], G.edges[e]['name']) for e in G.edges(v, keys=True) if e[0] != u and e[1] != u]
    m = {u:next_name, v:next_name}
    #new_edges = [(m[e[0]], m[e[1]], e[2], e[3]) for e in list(u_out_edges) + list(v_out_edges)]
    
    # Remove nodes and all edges

    if distances:
        min_in_length = min([G.edges[e]['len'] for e in in_edges])
        # new_edges_u = [(m.get(e[0], e[0]), m.get(e[1], e[1]), e[2], e[3], G.edges[(e[0], e[1], e[2])]['len']) for e in list(u_out_edges)]
        # new_edges_v = [(m.get(e[0], e[0]), m.get(e[1], e[1]), e[2], e[3], G.edges[(e[0], e[1], e[2])]['len']) for e in list(v_out_edges)]

        new_edges_u = [(m.get(e[0], e[0]), m.get(e[1], e[1]), e[2], e[3], G.edges[(e[0], e[1], e[2])]) for e in list(u_out_edges)]
        new_edges_v = [(m.get(e[0], e[0]), m.get(e[1], e[1]), e[2], e[3], G.edges[(e[0], e[1], e[2])]) for e in list(v_out_edges)]
        G.add_node(next_name)
        if names:
            if nx.degree(G, u) >= 3 and len([neigh for neigh in nx.neighbors(G, v)]) == 2:#nx.degree(G, v) == 2:
                #edges = G.nodes[u]['edges']
                to_edge = [neigh for neigh in nx.neighbors(G, v) if neigh != u][0]
                new_edge = {next_name:G.edges[(u, v, 0)]['name'], to_edge:G.edges[(to_edge, v, 0)]['name']}
                # edges[(to_edge, 0)] = edges[(v, 0)]
                # for (n1, n2, k) in in_edges:
                #     del edges[(v, k)]
                # G.nodes[next_name]['edges'] = edge
                

                # from_edges = G.nodes[to_edge]['edges']
                # from_edges[(next_name,0)] = from_edges[(v, 0)]
                # #del from_edges[(v, 0)]
                # G.nodes[to_edge]['edges'] = from_edges
                # for (n1, k) in from_edges.copy():
                #     if n1 == u or n1 == v:
                #         print('\t', (next_name, k) , '->', (n1, k))
                #         next_k = 0
                #         while (next_name, next_k) in 
                #         from_edges[(next_name, k)] = from_edges[(n1, k)]
                #         del from_edges[(n1, k)]


                # # Check
                # for (n1, k) in edges:
                #     assert n1 in G, 'Line 68, n1 NOT in G'
                # print('(u*, v, to_edge)', (u, v, to_edge))
                # for (n1, k) in from_edges:
                #     assert n1 in G, 'Line 71, n1 (' + str(n1) + ') NOT in G'
                
            elif nx.degree(G, v) >= 3 and len([neigh for neigh in nx.neighbors(G, u)]) == 2:#nx.degree(G, u) == 2:
                #edges = G.nodes[v]['edges']
                to_edge = [neigh for neigh in nx.neighbors(G, u) if neigh != v][0]
                new_edge = {next_name:G.edges[(u, v, 0)]['name'], to_edge:G.edges[(to_edge, u, 0)]['name']} # For edge: to_edge, next_name
                # edges[(to_edge, 0)] = edges[(u, 0)]
                # for (n1, n2, k) in in_edges:
                #     del edges[(u, k)]
                # G.nodes[next_name]['edges'] = edges
                

                # from_edges = G.nodes[to_edge]['edges']
                # from_edges[(next_name,0)] = from_edges[(u, 0)]
                # #del from_edges[(u, 0)]
                # G.nodes[to_edge]['edges'] = from_edges
                # for (n1, k) in from_edges.copy():
                #     if n1 == v:
                #         from_edges[(next_name, k)] = from_edges[(n1, k)]
                #         del from_edges[(n1, k)]

                # # Check
                # for (n1, k) in edges:
                #     assert n1 in G, 'Line 87, n1 NOT in G' 
                # print('(u, v*, to_edge)', (u, v, to_edge))
                # for (n1, k) in from_edges:
                #     assert n1 in G, 'Line 93, n1 NOT in G'
            # for e in in_edges:
            #     G.nodes[next_name]['edges'][e] = G.edges[e]['name']
            else: 
                to_edge = False
        
        if udeg <= 2:
            min_length_u = min_in_length
            min_length_v = 0
            if vdeg >= 3: G.nodes[next_name]['junction'] = G.nodes[v]['junction']
        else:
            min_length_u = 0
            min_length_v = min_in_length
            if udeg >= 3: G.nodes[next_name]['junction'] = G.nodes[u]['junction']

        # re-add node and add new edges
        G.remove_nodes_from([u, v])

        for (n1, n2, k, name, props) in new_edges_u:
            length = props['len']
            k2 = G.add_edge(n1, n2, name=name, len=length+min_length_u)
            if to_edge and (to_edge in (n1, n2) and next_name in (n1, n2)):
                G.edges[(n1, n2, k2)]['names'] = new_edge
            elif 'names' in props:
                G.edges[(n1, n2, k2)]['names'] = {m.get(n, n):props['names'][n] for n in props['names']}

        for (n1, n2, k, name, props) in new_edges_v:
            length = props['len']
            k2 = G.add_edge(n1, n2, name=name, len=length+min_length_v)
            if to_edge and (to_edge in (n1, n2) and next_name in (n1, n2)):
                G.edges[(n1, n2, k2)]['names'] = new_edge
            elif 'names' in props:
                G.edges[(n1, n2, k2)]['names'] = {m.get(n, n):props['names'][n] for n in props['names']}

    else:
        new_edges = [(m.get(e[0], e[0]), m.get(e[1], e[1]), e[2], e[3]) for e in list(u_out_edges) + list(v_out_edges)]
        # re-add node and add new edges
        G.add_node(next_name)
        if names:
            if nx.degree(G, u) >= 3:
                edges = G.nodes[u]['edges']
            else:
                edges = G.nodes[v]['edges']
            print(next_name, '->', edges)
            G.nodes[next_name]['edges'] = edges
            for e in in_edges:
                G.nodes[next_name]['edges'][e] = G.edges[e]['name']
        G.remove_nodes_from([u, v])
        
        for (n1, n2, k, name) in new_edges:
            G.add_edge(n1, n2, name=name)

    

    #print('\t\t', [G.edges[e]['name'] for e in G.edges(next_name, keys=True)])
    if 'root' in G.graph and (G.graph['root'] == u or G.graph['root'] == v):
        G.graph['root'] = next_name
    
    return next_name + 1


def priority_cost(G, e, ctype):
    """ Given a graph G and an edge e, compute the cost """
    if ctype == "outdegree":
        overcomplete_degree = nx.degree(G, e[0]) + nx.degree(G, e[1]) 
        in_degree = G.number_of_edges(e[0], e[1])
        cost = overcomplete_degree - 2*in_degree
    elif ctype == "indegree":
        in_degree = G.number_of_edges(e[0], e[1])
    elif ctype == "product":
        in_degree = G.number_of_edges(e[0], e[1])
        cost = (nx.degree(G, e[0]) - in_degree) * (nx.degree(G, e[1]) - in_degree)
    return cost


def priority_merge(G, ctype="outdegree"):
    """ Bottom-up merge sequence computation """
    # Inititalize 
    next_name = None 
    merge_seq = [[-1, -1, None]]*(G.number_of_nodes()-1)
    parents = [-1]*G.number_of_nodes() + [-1]*(G.number_of_nodes()-1)
    outnodes = [-1]*G.number_of_nodes() + [-1]*(G.number_of_nodes()-1)
    next_pair = PriorityQueue()
    exists = [1]*G.number_of_nodes() + [-1]*(G.number_of_nodes()-1)
    for e in G.edges(keys=True):
        w = priority_cost(G, e, ctype)
        next_pair.put((w, e))
    
    # Compute merge sequence
    for i in range(len(merge_seq)):
        # Pick u, v:
        (w, (u, v, k)) = next_pair.get()
        while not exists[u-1] or not exists[v-1]:
            (w, (u, v, k)) = next_pair.get()

        # Merge
        next_name = merge(G, u, v, next_name)
        last_name = next_name - 1

        exists[u-1] = 0
        exists[v-1] = 0
        exists[last_name-1] = 1

        merge_seq[i] = [u, v, last_name]
        parents[u-1] = last_name
        parents[v-1] = last_name
        outnodes[last_name-1] = [u, v]

        # Add to merge sequence
        for e in G.edges(last_name, keys=True):
            w = priority_cost(G, e, ctype)
            next_pair.put((w, e))
    
    return merge_seq


# Pivoting:

import scipy.linalg as la
import scipy.sparse.linalg as sp_la

def iter_trigraph(G, merge_seq):
    """ Undirected trigraph computation """
    merge_seq = np.array(merge_seq)
    (R, C) = merge_seq.shape
    if C == 2: 
        N = G.number_of_nodes()
        merge_seq = np.hstack((merge_seq, np.arange(N+1, 2*N).reshape((-1, 1)) ))
    edges1 = merge_seq[:, [0, 2]]
    edges2 = merge_seq[:, [1, 2]]
    edges = np.vstack((edges1, edges2))
    T = nx.Graph()
    T.add_edges_from(edges)
    
    # Fix the final edge, as trigraphs should be rootless
    T.add_edge(merge_seq[-1][0], merge_seq[-1][1])
    re = (merge_seq[-1][0], merge_seq[-1][1])
    T.remove_node(merge_seq[-1][2])
    
    return T, re

def merge_tree(G, T, re, N=None, rename_nodes=False):
    """ Given a tree and a root edge, return merge tree """
    if N is None:
        N = max([n for n in T.nodes() if nx.degree(T, n) == 1])
    
    # Insert a root node:
    r = max(T.nodes())+1
    T.add_edge(r, re[0])
    T.add_edge(r, re[1])
    T.remove_edge(re[0], re[1])
    Mt = nx.dfs_tree(T, r)
    Mt.graph['root'] = r
    
    new_names = {n:n for n in Mt.nodes()}
    if rename_nodes:
        bfs_ord = reversed(bfs_nodes(Mt, r))
        it = 1
        for n in bfs_ord:
            if Mt.out_degree(n) == 2:
                new_names[n] = N + it
                it += 1
    nx.relabel_nodes(Mt, new_names, copy=False)
    
    return Mt, new_names

def merge_sequence(Mt, N, r):
    """ Given a merge tree, compute merge sequence """
    S = reversed(bfs_nodes(Mt, r))
    merge_seq = [0]*(N-1)
    nmap = {n : n for n in range(1, N+1)}
    it = 0
    for n in S:
        if Mt.out_degree(n) == 2:
            [u, v] = Mt.successors(n)
            merge_seq[it] = [nmap[u], nmap[v], N+it+1]
            nmap[n] = N+it+1
            it += 1
    return merge_seq


def bfs_nodes(Mt, r):
    Q = Queue()
    Q.put(r)
    S = [r]
    while not Q.empty():
        n = Q.get()
        succ = list(Mt.successors(n))
        if len(succ) == 2:
            [u, v] = succ
            S += [u, v]
            Q.put(u)
            Q.put(v)
    
    return S


def cost_function(n):
    if 'cost function' not in meta or meta['cost function'] == 'outdegree':
        return n
    elif meta['cost function'] == 'bell':
        return float(cost_function(n))
    else:
        return float(pg.number_directed_partitions(n))


def cost_product(outd_u, outd_w, outd_n=0):
    """
    Given two nodes u, w merged into n, and their degrees,
    compute the cost of merging them
    """
    w1 = 1
    w2 = 0
    w3 = 0
    w4 = 0
    #return cost_function(outd_u)*cost_function(outd_w)

    logn = max(np.log(cost_function(outd_n)), 0)
    return cost_function(outd_u)*cost_function(outd_w)*\
        (w1*(outd_u + outd_w) + w2*logn+w3)+w4


def obtain_root(Mt_d, r, g=None):
    """
    Compute best root, given merge tree Mt
    """
    Mt = nx.Graph(Mt_d)
    # Step 1: compute total cost of Mt (to do)

    total_cost_r = 0

    # Step 2 iterate over edges - propagating cost
    # 2.1 label root edges
    for e in Mt.edges(r):
        Mt.edges[e]['root_cost'] = total_cost_r
        Mt.nodes[e[0]]['base_cost'] = total_cost_r
        Mt.nodes[e[1]]['base_cost'] = total_cost_r
        Mt.nodes[e[1]]['base_neighbor'] = e[0]
        Mt.nodes[e[0]]['base_neighbor'] = e[1]


    Q = queue.Queue()
    diffs = []
    for neigh in nx.neighbors(Mt, r):
        if nx.degree(Mt, neigh) == 3:
            Q.put((0, (r, neigh)))
        diffs.append((0, (r, neigh)))
    # Q = list(nx.bfs_successors(Mt), source=r)
    # Q = [q for q in Q if len(q[1]) == 2 and q[0] != r]
    while not Q.empty():
        (diff, (a, x)) = Q.get()
        [b, c] = [neigh for neigh in nx.neighbors(Mt, x) if neigh != a] # If 1 neighbor?
        Ca = Mt.edges[(x, a)]['cutsize']
        Cb = Mt.edges[(x, b)]['cutsize']
        Cc = Mt.edges[(x, c)]['cutsize']
        diffs1 = diff - cost_product(Cc, Cb, Ca) - cost_product(Ca, Ca, 0) + \
            cost_product(Cc, Ca, Cb) + cost_product(Cb, Cb, 0)
        diffs2 = diff - cost_product(Cb, Cc, Ca) - cost_product(Ca, Ca, 0) + \
            cost_product(Cb, Ca, Cc) + cost_product(Cc, Cc, 0)
        diffs += [(diffs1, (x, b)), (diffs2, (x, c))]
        if nx.degree(Mt, b) == 3:
            Q.put((diffs1, (x, b)))
        if nx.degree(Mt, c) == 3:
            Q.put((diffs2, (x, c)))
        #tri_map[39] = 9999
        # print("{" + str(tri_map[x]) + ', ' + str(tri_map[b]) + ', ' + str(diff) + ")-> a:(" + str(tri_map[a]) + \
        #     " " + str(Ca) + ")-> b:(" + str(tri_map[b]) + " " + str(Cb) + ")-> c:(" + str(tri_map[c]) + " " + str(Cc) + ") -> " + str(diffs1))
        # print("{" + str(tri_map[x]) + ', ' + str(tri_map[c]) + ', ' + str(diff) + ")-> a:(" + str(tri_map[a]) + \
        #     " " + str(Ca) + ")-> b:(" + str(tri_map[c]) + " " + str(Cc) + ")-> c:(" + str(tri_map[b]) + " " + str(Cb) + ") -> " + str(diffs2))

    (min_rc, min_edge) = min(diffs)



    # CUT SIZE CHECK 2
    
    # for e in Mt.edges():
    #     cutsize1 = Mt.edges[e]['cutsize']
    #     mt = Mt.copy()
    #     mt.remove_edge(e[0], e[1])
    #     conn = list(nx.connected_components(mt))
    #     mt1 = nx.subgraph(mt, conn[0])
    #     mt2 = nx.subgraph(mt, conn[1])
    #     leaves1 = [n for n in mt1 if nx.degree(mt1, n) <= 1]
    #     leaves2 = [n for n in mt2 if nx.degree(mt2, n) <= 1]
    #     cutsize2 = nx.cut_size(g, leaves1, leaves2)
    #     assert cutsize2 == cutsize1, str(cutsize1) + ', '+ str(cutsize2)

    if min_rc == 0:
        return Mt_d, 0
    
        
    # # 2.2 propagate:
    # edges = nx.bfs_edges(Mt, r)
    # edges = [e for e in edges if r not in e]
    # for e in edges:
    #     # base_cost is F(Ta)
    #     if 'base_cost' in Mt.nodes[e[0]]:
    #         x = e[0]
    #         neigh_b = e[1]
    #     else:
    #         x = e[1]
    #         neigh_b = e[0]
    #     base_cost = Mt.nodes[x]['base_cost']
    #     neigh_a = Mt.nodes[x]['base_neighbor']
    #     Ca  = Mt.edges[(x, neigh_a)]['cutsize'] # TODO - where is this value? - it's on the edge, not node
    #     neigh_c = [n for n in nx.neighbors(Mt, x) if n != neigh_a and n != neigh_b][0]
    #     Cc = Mt.edges[(x, neigh_c)]['cutsize']
    #     Cb = Mt.edges[(x, neigh_b)]['cutsize']

    #     tri_map[39] = 9999
    
    #     root_cost = base_cost - cost_product(Cc, Cb, Ca) - cost_product(Ca, Ca, 0) + \
    #         cost_product(Cc, Ca, Cb) + cost_product(Cb, Cb, 0)
        
    #     print("{" + str(tri_map[x]) + ', ' + str(tri_map[neigh_b]) + ', ' + str(base_cost) + ")-> a:(" + str(tri_map[neigh_a]) + \
    #         " " + str(Ca) + ")-> b:(" + str(tri_map[neigh_b]) + " " + str(Cb) + ")-> c:(" + str(tri_map[neigh_c]) + " " + str(Cc) + ") -> " + str(root_cost))
    #     Mt.edges[e]['root_cost'] = root_cost
    #     Mt.nodes[neigh_b]['base_cost'] = root_cost
    #     Mt.nodes[neigh_b]['base_neighbor'] = x

    # Obtain minimum cost edge
    #(min_rc, min_edge) = min([(Mt.edges[e]['root_cost'], e) for e in Mt.edges()])
    
    # Remove root:
    (u, v) = nx.neighbors(Mt, r)
    cut = Mt.edges[(r, u)]['cutsize']
    Mt.remove_node(r)
    Mt.add_edge(u,v, cutsize=cut)

    # Add root to minimum cost edge: (where does e come from?)
    cut = Mt.edges[(min_edge[0], min_edge[1])]['cutsize']
    Mt.remove_edge(min_edge[0], min_edge[1])
    Mt.add_edge(min_edge[0], r, cutsize=cut)
    Mt.add_edge(min_edge[1], r, cutsize=cut)

    #print((tri_map[min_edge[0]], tri_map[min_edge[1]]))
    Mt = nx.dfs_tree(Mt, r)
    return Mt, min_rc



def compute_costs_l(G, Mt, r=None, N=None):
    """ Region outdegree computation using numpy arrays (slow) """
    if N is None:
        N = G.number_of_nodes()
    if r is None:
        r = [n for n in Mt.nodes() if Mt.in_degree(n) == 0][0]
    
    # Get cost
    # Obtain a reverse BFS ordering in S
    S = bfs_nodes(Mt, r)
    
    # From the reverse BFS ordering (S), compute cost
    total_cost = 0
    while S != []:
        p = S.pop()
        
        succ = list(Mt.successors(p))
        if len(succ) == 2:
            [c1, c2] = succ
            # The cost between c1 and c2 
            b1 = Mt.nodes[c1]['boundary']
            b2 = Mt.nodes[c2]['boundary']
            cl1 = Mt.nodes[c1]['cluster']
            cl2 = Mt.nodes[c2]['cluster']

            Mt.nodes[p]['boundary'] = (b1 + b2)*np.logical_not(cl1 + cl2)#np.logical_and(b1, cl2)
            Mt.nodes[p]['cluster'] = cl1+cl2
            Mt.nodes[p]['outdegree'] = Mt.nodes[p]['boundary'].sum()
            Mt.nodes[p]['cost'] = float(cost_function(Mt.nodes[c1]['outdegree']))*float(cost_function(Mt.nodes[c2]['outdegree']))#Mt.nodes[p]['boundary'].sum()
        else:
            Mt.nodes[p]['boundary'] = np.zeros((N,))
            Mt.nodes[p]['boundary'][[neigh-1 for neigh in nx.neighbors(G, p)]] = 1
            Mt.nodes[p]['cluster'] = np.zeros((N,))
            Mt.nodes[p]['cluster'][p-1] = 1
            Mt.nodes[p]['outdegree'] = Mt.nodes[p]['boundary'].sum()
            Mt.nodes[p]['cost'] = float(cost_function(nx.degree(G, p)))

        total_cost += Mt.nodes[p]['cost']
    return Mt, total_cost


def compute_cost_change_l(Mt, u, v, a, b, c):
    xclust = Mt.nodes[a]['cluster'] + Mt.nodes[b]['cluster']
    xbound = (Mt.nodes[a]['boundary']+Mt.nodes[b]['boundary'])*np.logical_not(Mt.nodes[a]['cluster']+Mt.nodes[b]['cluster'])

    xout = xbound.sum()
    costx = xout*Mt.nodes[c]['outdegree']

    bx = float(cost_function(xout))
    ba = float(cost_function(Mt.nodes[a]['outdegree']))
    bb = float(cost_function(Mt.nodes[b]['outdegree']))
    bc = float(cost_function(Mt.nodes[c]['outdegree']))
    bv = float(cost_function(Mt.nodes[v]['outdegree']))
    
    # new_cost = (bx*bc + ba*bb)
    # old_cost = (ba*bv + bb*bc)

    new_cost = cost_product(xout,Mt.nodes[c]['outdegree']) + cost_product(Mt.nodes[a]['outdegree'],Mt.nodes[b]['outdegree'])
    old_cost = cost_product(Mt.nodes[a]['outdegree'],Mt.nodes[v]['outdegree']) + cost_product(Mt.nodes[b]['outdegree'],Mt.nodes[c]['outdegree'])

    cost_change = new_cost - old_cost
    
    return (cost_change, xbound, xclust, xout, costx, (u,v,a,b,c))


def node_cost_candidates(Mt, ccl, i):
    
    u = i+1
    ccl[i, :] = np.zeros((4,))
    if Mt.out_degree(i+1) == 2:
        ccc = compute_cost_change_l
        [u1, u2] = sorted(Mt.successors(i+1))

        # Between b and c, which should we give to u?
        cchange = [(0,), (0,), (0,), (0,)]
        if Mt.out_degree(u1) == 2:
            v = u1
            a = u2
            [b, c] = sorted(Mt.successors(u1))
            ccl[i, 0] = ccc(Mt, u, v, a, b, c)[0]
            ccl[i, 1] = ccc(Mt, u, v, a, c, b)[0]
        if Mt.out_degree(u2) == 2:
            v = u2
            a = u1
            [b, c] = sorted(Mt.successors(u2))
            ccl[i, 2] = ccc(Mt, u, v, a, b, c)[0]
            ccl[i, 3] = ccc(Mt, u, v, a, c, b)[0]


def compute_cost_change_list(Mt):
    """ Compute a Nx4 matrix with all possible cost changes """
    N = Mt.number_of_nodes()
    ccl = np.zeros((N,4))
    for i in range(N):
        node_cost_candidates(Mt, ccl, i)

    return ccl


def pivot_global(Mt, ccl):
    """ Globally determine the best way to pivot a merge tree """
    [n, i] = divmod(ccl.argmin(), ccl.shape[1])
    n += 1
    while ccl[n-1,i] < 0:
        #print(n, '\t->', ccl[n-1, i])
        pivot_local(Mt, ccl, n, i)
        [n, i] = divmod(ccl.argmin(), ccl.shape[1])
        n += 1


def pivot_local(Mt, ccl, u, i):
    """ Given a Cost Change List ccl, a node u and index i, do the pivot based on the index i """
    # Arange the nodes based on i
    [u1, u2] = sorted(Mt.successors(u))
    if i == 0:
        v = u1
        a = u2
        [b, c] = sorted(Mt.successors(v))
    elif i == 1:
        v = u1
        a = u2
        [c, b] = sorted(Mt.successors(v))
    elif i == 2:
        v = u2
        a = u1
        [b, c] = sorted(Mt.successors(v))
    else:
        v = u2
        a = u1
        [c, b] = sorted(Mt.successors(v))

    # deal with parent of subtree edges, n0
    parent = False
    if Mt.in_degree(u) == 1:  
        n0 = next(Mt.predecessors(u))
        parent = True

    Mt.add_edges_from([(u, c), (v, a)])
    Mt.remove_edges_from([(u, a), (v, c)])
    
    (cost_change, xbound, xclust, xout, costx, (u, v, a, b, c)) = compute_cost_change_l(Mt, u, v, a, b, c)  # REDUNDANT

    Mt.nodes[v]['cluster'] = xclust
    Mt.nodes[v]['boundary'] = xbound
    Mt.nodes[v]['outdegree'] = xout
        
    # Update costs
    if parent: node_cost_candidates(Mt, ccl, n0-1)
    node_cost_candidates(Mt, ccl, u-1)
    node_cost_candidates(Mt, ccl, v-1)


def pivot(Mt, u, l=True):
    """ 
    Given a parent node u with successor v, make v the parent, u the successor
    b will be handed from v to u\
    Possibilities with these variables when pivoting u and v are:
     - u, v, a, b, c
     - u, v, a, c, b
    """
    ccc = compute_cost_change_l
    
    if Mt.out_degree(u) == 2:
        [u1, u2] = Mt.successors(u)
        
        # Between b and c, which should we give to u?
        cchange = [(0,), (0,), (0,), (0,)]
        if Mt.out_degree(u1) == 2:
            v = u1
            a = u2
            [b, c] = Mt.successors(u1)
            cchange[0] = ccc(Mt, u, v, a, b, c)
            cchange[1] = ccc(Mt, u, v, a, c, b)
        if Mt.out_degree(u2) == 2:
            v = u2
            a = u1
            [b, c] = Mt.successors(u2)
            cchange[2] = ccc(Mt, u, v, a, b, c)
            cchange[3] = ccc(Mt, u, v, a, c, b)

        best_pivot = np.argmin([cc[0] for cc in cchange])

        if cchange[best_pivot][0] >= 0:
            cost_change = 0
        else:
            (cost_change, xbound, xclust, xout, costx, (u, v, a, b, c)) = cchange[best_pivot]
            # deal with parent of subtree edges, n0
            if Mt.in_degree(u) == 1:  
                n0 = next(Mt.predecessors(u))
                Mt.remove_edge(n0, u)
                Mt.add_edge(n0, v)

            # Change root if necessary
            if u == Mt.graph['root']:
                Mt.graph['root'] = v

            # Deal with u edges
            Mt.remove_edge(u, v)
            Mt.add_edge(u, b)

            # Deal with v edges
            Mt.remove_edge(v, b)
            Mt.add_edge(v, u)

            # Change both the outdegree and clusters of u and v:
            #prevcost = Mt.nodes[u]['cost']
            Mt.nodes[v]['cluster'] = Mt.nodes[u]['cluster']
            Mt.nodes[v]['boundary'] = Mt.nodes[u]['boundary']
            Mt.nodes[v]['outdegree'] = Mt.nodes[u]['outdegree']
            Mt.nodes[u]['cluster'] = xclust
            Mt.nodes[u]['boundary'] = xbound
            Mt.nodes[u]['outdegree'] = xout
    else:
        cost_change = 0

    return cost_change


def to_mathematica(G, Mt, root):
    # Convert G
    coords = []
    for n in sorted(G):
        if 'x' in G.nodes[n]:
            coords.append([G.nodes[n]['x'], G.nodes[n]['y']])
        else:
            coords.append([0.0, 0.0])
    adj_list = [[u, v] for (u, v) in G.edges()]
    merge_tree = []
    rindex = -1
    for n in range(1, 2*G.number_of_nodes()):
        if n <= G.number_of_nodes():
            merge_tree.append([3, n, int(Mt.nodes[n]["outdegree"])])
        elif n == root:
            succ = [int(m) for m in Mt.successors(n)]
            rindex = n
            merge_tree.append([1, succ, 0])
        else:
            succ = [int(m) for m in Mt.successors(n)]
            #boundary = (np.where(Mt.nodes[n]['boundary'])[0]+1).tolist()
            merge_tree.append([2, succ, int(Mt.nodes[n]['outdegree']) ])
    return [rindex, merge_tree], [coords, adj_list]






def outd_n(G):
    """ Naive """
    g = G.copy()
    def f(e): return (nx.degree(g, e[0])+nx.degree(g, e[1]) - 2*g.number_of_edges(e[0], e[1]))
    N = g.number_of_nodes()
    

    seq = [[0, 0, i] for i in range(N+1, 2*N)]
    for i in range(len(seq)):
        #print(i, ' /', G.number_of_nodes(), ' :: ', len(seq))
        (u,v,k) = min(g.edges, key=f)
        seq[i] = [u, v, seq[i][2]]
        merge(g, u, v, next_name=seq[i][2])
    
    return seq

def outd(G):
    """ Faster """
    g = G.copy()
    def f(e): return c[tuple(sorted([e[0], e[1]]))]
    c = {
        tuple(sorted(e)) : (
            cost_function(nx.degree(g, e[0])+nx.degree(g, e[1]) - 2*g.number_of_edges(e[0], e[1]))
            ) for e in G.edges()
    }
    N = g.number_of_nodes()
    
    seq = [[0, 0, i] for i in range(N+1, 2*N)]
    for i in range(len(seq)):
        (u,v,k) = min(g.edges, key=f)
        seq[i] = [u, v, seq[i][2]]
        merge(g, u, v, next_name=seq[i][2])
        del c[tuple(sorted([u, v]))]
        for e in g.edges(seq[i][2]):
            ep = tuple(sorted(e))
            c[ep] = cost_function(nx.degree(g, e[0])+nx.degree(g, e[1]) - 2*g.number_of_edges(e[0], e[1]))

    return seq


def prod(G):
    """ Faster """
    g = G.copy()
    def f(e): return c[tuple(sorted([e[0], e[1]]))]
    # c = {
    #     tuple(sorted(e)) : (nx.degree(g, e[0])*nx.degree(g, e[1])) for e in G.edges()
    # }
    c = {
        tuple(sorted(e)) : cost_product(nx.degree(g, e[0]),nx.degree(g, e[1]),0) for e in G.edges()
    }
    N = g.number_of_nodes()
    
    seq = [[0, 0, i] for i in range(N+1, 2*N)]
    for i in range(len(seq)):
        (u,v,k) = min(g.edges, key=f)
        seq[i] = [u, v, seq[i][2]]
        merge(g, u, v, next_name=seq[i][2])
        del c[tuple(sorted([u, v]))]
        for e in g.edges(seq[i][2]):
            ep = tuple(sorted(e))
            # c[ep] = (nx.degree(g, e[0])*nx.degree(g, e[1]))
            c[ep] = cost_product(nx.degree(g, e[0]),nx.degree(g, e[1]),0)

    return seq



# TESTING ALTERNATIVES TO NETWORKX
def get_fied(G, nodelist, nmap, normalized=False):

    #rev_nmap = {nmap[k]:k for k in nmap}
    #nodelist = [rev_nmap[n] for n in sorted(rev_nmap)]

    # print("G stats", G.number_of_nodes(), G.number_of_edges(), nx.number_connected_components(G))
    if 'laplacian' not in G or not STORE_ADJ:
        if normalized:
            L = nx.normalized_laplacian_matrix(G, nodelist)
        else:
            L = nx.laplacian_matrix(G, nodelist).astype(np.float)
        G.graph['laplacian'] = L
    else:
        t0 = time.time()
        L = G.graph['laplacian']
        ADJ_TIMES['lap'] += time.time() - t0
    # # Alt 1: Using eigsh: Non singular matrices fail completely because of sigma=0 parameter...
    # # ERROR FACTOR IS EXACTLY SINGULAR
    if mode == 2:
        # print(L.shape, np.linalg.det(L.todense()), end=' ')
        # if L.shape == (13,13):
        #     print(L.todense())
        try:
            eig_values, v= scipy.sparse.linalg.eigsh(L, 2, sigma=0.000001, which='LM')
        except TypeError:
            eig_values, v= scipy.linalg.eigh(L.toarray())
        ev = v[:,1]
        sd = ev.argsort().tolist()
    # # Alt 2: Problems with sparsity in small matrices... problems with normalized Laplacian
    elif mode ==3:
        if G.number_of_nodes() <= 12: L = L.toarray()
        eig_values, eig_vectors = scipy.sparse.linalg.eigs(L)
        fiedler_pos = np.where(eig_values.real == np.sort(eig_values.real)[1])[0][0]
        fiedler_vector = np.transpose(eig_vectors)[fiedler_pos]
        sd = fiedler_vector.argsort().tolist()
    elif mode == 4:
        nodelist2 = [n for n in G.nodes()]
        nmap2 = {nodelist2[i]:i for i in range(len(nodelist2))}
        t0 = time.time()
        fied = nx.fiedler_vector(G, method='tracemin_lu')
        global FIEDLER_TIME
        FIEDLER_TIME += time.time() - t0
        sd = fied.argsort().tolist()
        return [nodelist2[i] for i in sd]

    # # Alt 2: Using nx spectral ordering - Completes well... but slow. Also, freezes with small matrices without tracemin_lu
    # sd = [nidx[n] for n in nx.spectral_ordering(G, normalized=normalized, method='tracemin_lu', seed=1)]
    
    # Return NODES
    return [nodelist[i] for i in sd]

def topdown_split(g, a, nodes, nmap, normalized=False, method=1):
    # print(g.number_of_nodes())
    if mode == 1:
        fied = list(nx.spectral_ordering(g, normalized=normalized, method='tracemin_lu', seed=1))
    else:
        t0 = time.time()
        fied = get_fied(g, nodes, nmap, normalized)
        global SPECTRAL_TIME
        SPECTRAL_TIME += time.time() - t0
    #fied = list(nx.spectral_ordering(g, normalized=normalized, seed=1))
    fm = [nmap[f] for f in fied]
    N = len(fm)
    cut_sizes = np.ones((N,))*np.inf
    if method == 1:
        # MOD
        low = int(np.ceil(N/3))
        high = N-low
        for i in range(low, high+1):
            cut_sizes[i] = a[np.ix_(fm[:i], fm[i:])].sum()/2
        # /MOD
        # for i in range(int(N/3)+1, N-int(N/3)):
        #     cut_sizes[i] = a[np.ix_(fm[:i], fm[i:])].sum()/2

        c = cut_sizes.argmin()
        # print("Index:", c)
        nodes1 = fied[:c]
        nodes2 = fied[c:]
    elif method == 2:
        global CUT_TIME
        low = int(np.ceil(N/3))
        high = N-low
        source_set = fied[:low]
        target_set = fied[high:]
        center_set = fied[low:high]

        t0 = time.time()
        nodes1, nodes2, c = set_minimum_edge_cut(g, a, nmap, source_set, target_set, center_set)
        CUT_TIME += time.time() - t0

    return nodes1, nodes2, c


def set_minimum_edge_cut(g, a, nmap, source_set, target_set, center_set):
    """ Min cut from source set to target set. """
    # Base case - size 4 graph is just cut between source and target
    if center_set == []:
        #c = np.count_nonzero(a[np.ix_([nmap[n] for n in source_set], [nmap[n] for n in target_set])])
        c = np.sum(a[np.ix_([nmap[n] for n in source_set], [nmap[n] for n in target_set])])
        return source_set, target_set, c

    M = max(g.nodes())
    v1, v2 = M+1, M+2

    # Create a virtual vertex representing the source set, v1
    sgm = create_virtual_vertex(g, a, nmap, center_set+target_set, source_set, v1)
    nmap[v1] = g.number_of_nodes()
    
    it = 0
    nmap2 = {}
    for k in sorted(nmap):
        if k in sgm:
            nmap2[k] = it
            it += 1 
    # a = nx.adjacency_matrix(sgm, nmap)
    a = nx.adjacency_matrix(sgm, nmap2) # MOD
    #  = nx.adjacency_matrix(sgm, {k:v for (k,v) in nmap.items() if k in sgm})

    # Create a virtual vertex representing the target set, v2
    create_virtual_vertex(sgm, a, nmap2, center_set + [v1], target_set, v2, copy=False)

    # Convert multigraph to normal graph
    sg = nx.Graph()
    sg.add_nodes_from(sgm.nodes())
    for e in sgm.edges():
        sg.add_edge(e[0], e[1], capacity=sgm.number_of_edges(e[0], e[1]))#weight=sgm.number_of_edges(e[0], e[1]))

    # Perform cut and retrieve original nodes from the virtual vertices
    val, partition = nx.minimum_cut(sg, v1, v2)#, capacity='weight')
    if v1 in partition[0]:
        nodes1 = list(partition[0] - {v1}) + source_set
        nodes2 = list(partition[1] - {v2}) + target_set
    else:
        nodes1 = list(partition[1] - {v1}) + source_set
        nodes2 = list(partition[0] - {v2}) + target_set

    return nodes1, nodes2, val

from scipy.sparse import csr_matrix
def create_virtual_vertex(g, a, nmap, set1, set2, v, copy=True):
    """ Replace set2 with virtual vertex """
    n1m = [nmap[n] for n in set1]
    n2m = [nmap[n] for n in set2]

    # Find edges set1 -> set2, replace with set1 -> v
    v_w = a[np.ix_(n1m, n2m)].sum(axis=1).flatten()

    if len(v_w.shape) > 1:
        v_n = np.where(v_w)[1]
    else:
        v_n = np.where(v_w)[0]
        v_w = v_w.reshape((1, v_w.size))
    edges = []
    for vn in v_n:
        #edges += [(set1map[n1m[vn]], v)]*v_w[0, vn]
        edges += [(set1[vn], v)]*v_w[0, vn]

    # Graph metadata
    if STORE_ADJ:
        t0 = time.time()
        a2 = a[np.ix_(n1m, n1m)]
        t1 = time.time()
        #a2 = csr_matrix((a2.data, a2.indices, np.pad(a2.indptr, (0,1), "edge")),  shape=(len(n1m)+1, len(n1m)+1)) # Padding for sparse matrices
        a2 = np.pad(a2, ((0,1),(0,1)))
        #a2 = csr_matrix((a2.data, a2.indices, a2.indptr), shape=(len(n1m)+1, len(n1m)+1))
        a2[-1,:-1] = v_w
        a2[:-1,-1] = v_w

        # laplacian
        l = g.graph['laplacian']
        l2 = a[np.ix_(n1m, n1m)]
        t1 = time.time()
        #l2 = csr_matrix((l2.data, l2.indices, np.pad(l2.indptr, (0,1), "edge")),  shape=(len(n1m)+1, len(n1m)+1)) # Padding for sparse matrices
        l2 = np.pad(l2, ((0,1),(0,1)))
        ADJ_TIMES['shuffle2'] += time.time() - t1
        #a2 = csr_matrix((a2.data, a2.indices, a2.indptr), shape=(len(n1m)+1, len(n1m)+1))
        l2[-1,:-1] = -v_w
        l2[:-1,-1] = -v_w
        l2[-1,-1] = v_w.sum()

        # n1m: a2 ix -> a1 ix
        # nmap: node -> a1 ix
        # want: sg node -> a2 ix
        nmap2 = {n1m[i]:i for i in range(len(n1m))} # a1 ix -> a2 ix
        nmap2 = {n:nmap2[nmap[n]] for n in set1} # sg node -> a2 ix
        nmap2[v] = a2.shape[0]-1
        set1 = set1 + [v]
        ADJ_TIMES['shuffle'] += time.time() - t0
    # Remove set2 and add set1 -> v edges
    if copy:
        sg = g.subgraph(set1).copy()
        sg.add_node(v)
        sg.add_edges_from(edges)
        if STORE_ADJ:
            t0 = time.time()
            sg.graph['nodes'] = set1
            sg.graph['adjacency'] = a2
            sg.graph['laplacian'] = l2
            sg.graph['nmap'] = nmap2
            ADJ_TIMES['shuffle'] += time.time() - t0
        # HERE: MODIFY a AND nmap TO CONTAIN v
        return sg
    else:
        g.remove_nodes_from(set2)
        g.add_node(v)
        g.add_edges_from(edges)
        if STORE_ADJ:
            t0 = time.time()
            g.graph['nodes'] = set1
            g.graph['adjacency'] = a2
            g.graph['laplacian'] = l2
            g.graph['nmap'] = nmap2
            ADJ_TIMES['shuffle'] += time.time() - t0
        return None



def topdown_merge_tree(g, normalized=False, method=1):
    N = g.number_of_nodes()

    # Create graph with all edges connected to N+1
    tr = nx.Graph()
    tr.add_edges_from([(n, N+1) for n in g])
    for n in g:
        tr.edges[(n, N+1)]['cutsize'] = nx.degree(g, n)
    tr.nodes[N+1]['subgraph'] = g.copy()

    q = Queue()
    if nx.degree(tr, N+1) == 3:
        u = min([(nx.degree(g,n), n) for n in g])[1]
        [v, w] = [n for n in g if n != u]
        tr.add_edges_from([(N+2, v), (N+2, w), (N+1, N+2)])
        tr.edges[(v, N+2)]['cutsize'] = tr.edges[(v, N+1)]['cutsize']
        tr.edges[(w, N+2)]['cutsize'] = tr.edges[(w, N+1)]['cutsize']
        tr.edges[(N+1, N+2)]['cutsize'] = nx.degree(g, u)
        tr.remove_edges_from([(v, N+1), (w, N+1)])
    else:
        q.put(N+1)

    vnode = N+2
    it = 1
    while not q.empty():
        # Extract graph and information
        v = q.get()
        sg = tr.nodes[v]['subgraph']
        if 'adjacency' not in sg.graph or not STORE_ADJ:
            nodes = sorted(list(sg))
            nmap = {nodes[i]:i for i in range(len(nodes))}
            a = nx.adjacency_matrix(sg, nodes).todense()#nx.adjacency_matrix(sg, nodes)
            sg.graph['laplacian'] = nx.laplacian_matrix(sg, nodes).astype(np.float).todense()#nx.laplacian_matrix(sg, nodes).astype(np.float)
        else:
            t0 = time.time()
            nodes = sg.graph['nodes']
            nmap = sg.graph['nmap']
            a = sg.graph['adjacency']
            ADJ_TIMES['adj'] += time.time() - t0

        # Balanced split
        nodes1, nodes2, cutsize = topdown_split(sg, a, nodes, nmap, normalized=normalized, method=method)


        # Modify tree
        cuts1 = {n1:tr.edges[(n1, v)]['cutsize'] for n1 in nodes1}
        cuts2 = {n2:tr.edges[(n2, v)]['cutsize'] for n2 in nodes2}
        tr.remove_node(v)
        if it != 1:
            tr.add_node(vnode,   subgraph=create_virtual_vertex(sg, a, nmap, nodes1, nodes2, vnode+1))  # Replace nodes2 with virtual vertex vnode+1
            tr.add_node(vnode+1, subgraph=create_virtual_vertex(sg, a, nmap, nodes2, nodes1, vnode))
            tr.add_edge(vnode, vnode+1, cutsize=cutsize)
            tr.add_edges_from([(vnode, n) for n in nodes1] + [(vnode+1, n) for n in nodes2])#, cutsize=cutsize) # Doesn't consider multigraph
            for v1 in cuts1:
                tr.edges[(vnode, v1)]['cutsize'] = cuts1[v1]
            for v2 in cuts2:
                tr.edges[(vnode+1, v2)]['cutsize'] = cuts2[v2]
        else:
            tr.add_node(vnode,   subgraph=create_virtual_vertex(sg, a, nmap, nodes1, nodes2, v))  
            tr.add_node(vnode+1, subgraph=create_virtual_vertex(sg, a, nmap, nodes2, nodes1, v))
            tr.add_edges_from([(vnode, v), (vnode+1, v)], cutsize=cutsize)
            tr.add_edges_from([(vnode, n) for n in nodes1] + [(vnode+1, n) for n in nodes2])#, cutsize=cutsize) # Doesn't consider multigraph
            for v1 in cuts1:
                tr.edges[(vnode, v1)]['cutsize'] = cuts1[v1]
            for v2 in cuts2:
                tr.edges[(vnode+1, v2)]['cutsize'] = cuts2[v2]


        
        # ##########
        # # Cut size check
        # ##########
        # cutsize1 = nx.cut_size(sg, nodes1, nodes2)
        # trSplit = tr.copy()
        # if it != 1: trSplit.remove_edge(vnode, vnode+1)
        # else: trSplit.remove_edges_from([(vnode, v), (vnode+1, v)])
        # comps = list(nx.connected_components(trSplit))
        # tr1 = nx.subgraph(trSplit, comps[0])
        # tr2 = nx.subgraph(trSplit, comps[1])
        # leaves1 = [n for n in tr1 if nx.degree(tr1, n) == 1]
        # leaves2 = [n for n in tr2 if nx.degree(tr2, n) == 1]
        # cutsize2 = nx.cut_size(g, leaves1, leaves2)
        # assert cutsize==cutsize1 and cutsize==cutsize2, "CUTSIZES: " + str(cutsize) + ", " + str(cutsize1) + ", " + str(cutsize2)



        for neigh in nx.neighbors(tr, vnode+1):
            if neigh > N and neigh != N+1:
                if v in tr.nodes[neigh]['subgraph']:
                    nx.relabel_nodes(tr.nodes[neigh]['subgraph'], {v:vnode+1}, copy=False)
                    if STORE_ADJ:
                        t0 = time.time()
                        sg_nodes = tr.nodes[neigh]['subgraph'].graph['nodes']
                        for i in range(len(sg_nodes)):
                            if sg_nodes[i] == v:
                                sg_nodes[i] = vnode+1
                        tr.nodes[neigh]['subgraph'].graph['nmap'][vnode+1] = tr.nodes[neigh]['subgraph'].graph['nmap'][v]
                        del tr.nodes[neigh]['subgraph'].graph['nmap'][v]
                        ADJ_TIMES['loops'] += time.time() - t0
        for neigh in nx.neighbors(tr, vnode):
            if neigh > N and neigh != N+1:
                if v in tr.nodes[neigh]['subgraph']:
                    nx.relabel_nodes(tr.nodes[neigh]['subgraph'], {v:vnode}, copy=False)
                    if STORE_ADJ:
                        t0 = time.time()
                        sg_nodes = tr.nodes[neigh]['subgraph'].graph['nodes']
                        for i in range(len(sg_nodes)):
                            if sg_nodes[i] == v:
                                sg_nodes[i] = vnode
                        tr.nodes[neigh]['subgraph'].graph['nmap'][vnode] = tr.nodes[neigh]['subgraph'].graph['nmap'][v]
                        del tr.nodes[neigh]['subgraph'].graph['nmap'][v]
                        ADJ_TIMES['loops'] += time.time() - t0

        if nx.degree(tr, vnode) > 3:
            q.put(vnode)
        if nx.degree(tr, vnode+1) > 3:
            q.put(vnode+1)
        vnode += 2
        it += 1

    # Create a directed merge tree
    mt = nx.dfs_tree(tr, N+1)
    nx.set_edge_attributes(mt, {e:tr.edges[e]['cutsize'] for e in mt.edges()}, 'cutsize')
    nodes = sorted(list(mt))
    relabel = {
        nodes[i]:(i+1) for i in range(len(nodes))
    }
    nx.relabel_nodes(mt, relabel, copy=False)

    return mt


def pre_merge(G, merge_deg_1=True, merge_one_2=True, names=False, distances=False):
    g = G.copy()
    mseq = []
    next_node = g.number_of_nodes()+1

    # Degree 1
    if merge_deg_1:
        q = [(u,v) for (u,v) in g.edges() if nx.degree(g, u) == 1 or nx.degree(g, v) == 1]
        while q != []:
            (u, v) = q.pop()

            if u in g and v in g and (nx.degree(g, u) == 1 or nx.degree(g, v) == 1):
                merge(g, u, v, next_name=next_node, distances=distances, names=names)
                mseq.append([u, v, next_node])
                next_node += 1
                q += [(next_node-1, neigh) for neigh in nx.neighbors(g, next_node-1) if nx.degree(g, neigh) == 1 or nx.degree(g, next_node-1) == 1]
    # print('Degree 1', len([n for n in g if nx.degree(g, n) == 1]))
    # print('Degree 2', len([n for n in g if nx.degree(g, n) == 2]))
    # print()
    # Both degree 2
    q = [(u,v) for (u,v) in g.edges() if nx.degree(g, u) == 2 and nx.degree(g, v) == 2]
    while q != []:
        (u, v) = q.pop()

        if u in g and v in g and (nx.degree(g, u) == 2 and nx.degree(g, v) == 2):

            if names:
                u_node = [n for n in nx.neighbors(g,u) if n!=v][0]
                v_node = [n for n in nx.neighbors(g,v) if n!=u][0]
                u_ename = [g.edges[e]['name'] for e in g.edges(u,keys=True) if e[0]!=v and e[1]!=v][0]
                v_ename = [g.edges[e]['name'] for e in g.edges(v,keys=True) if e[0]!=u and e[1]!=u][0]

            merge(g, u, v, next_name=next_node, distances=distances, names=names)

            if names:
                g.edges[(next_node, u_node, 0)]['name'] = u_ename
                g.edges[(next_node, v_node, 0)]['name'] = v_ename

            mseq.append([u, v, next_node])
            next_node += 1
            if nx.degree(g, next_node-1) == 2:
                q += [(next_node-1, neigh) for neigh in nx.neighbors(g, next_node-1) if nx.degree(g, neigh) == 2]
    # print('Degree 1', len([n for n in g if nx.degree(g, n) == 1]))
    # print('Degree 2', len([n for n in g if nx.degree(g, n) == 2]))
    # print()
    if names:
        for n in g:
            g.nodes[n]['edges'] = {}
            for e in g.edges(n, keys=True):
                g.edges[e]['names'] = {e[0]:g.edges[e]['name'], e[1]:g.edges[e]['name']}
                # if e[0] == n:
                #     if e[0] == 2485:
                #         print('n:', n, ', e0', e[0])
                #         print('    e1 e2:', e[1], e[2])
                #     g.nodes[n]['edges'][(e[1], e[2])] = g.edges[e]['name']
                # else:
                #     if e[1] == 2485:
                #         print('n:', n, ', e1', e[1])
                #         print('    e0 e2:', e[0], e[2])
                #     g.nodes[n]['edges'][(e[0], e[2])] = g.edges[e]['name']

    # One degree 2:
    if merge_one_2:
        q = [(u,v) for (u,v) in g.edges() if nx.degree(g, u) == 2 or nx.degree(g, v) == 2]
        while q != []:
            (u, v) = q.pop()

            #if u in g and v in g and ((nx.degree(g, u) == 2 and len(list(g.neighbors(u))) != 1) or (nx.degree(g, v) == 2 and len(list(g.neighbors(v))) != 1)):
            if u in g and v in g and (nx.degree(g, u) == 2 or nx.degree(g, v) == 2):
                # if names and nx.degree(g, u) == 2:
                #     node = [n for n in nx.neighbors(g,u) if n!=v][0]
                #     ename = [g.edges[e]['name'] for e in g.edges(u,keys=True) if e[0]!=v and e[1]!=v][0]
                #     edges = g.nodes[v]['edges']

                    
                #     edges[node] = ename
                #     g.nodes[next_node]['edges'] = edges
                # elif names and nx.degree(g, v) == 2:
                #     nodes = [n for n in nx.neighbors(g,v) if n!=u]
                #     enames = [g.edges[e]['name'] for e in g.edges(v,keys=True) if e[0]!=u and e[1]!=u]
                #     node = [n for n in nx.neighbors(g,v) if n!=u][0]
                #     ename = [g.edges[e]['name'] for e in g.edges(v,keys=True) if e[0]!=u and e[1]!=u][0]
                #     edges = g.nodes[u]['edges']

                #     merge(g, u, v, next_name=next_node, distances=distances)
                #     edges[node] = ename
                #     g.nodes[next_node]['edges'] = edges

                merge(g, u, v, next_name=next_node, distances=distances, names=names)

                mseq.append([u, v, next_node])
                next_node += 1
                q += [(next_node-1, neigh) for neigh in nx.neighbors(g, next_node-1) if nx.degree(g, neigh) == 2 or nx.degree(g, next_node-1) == 2]
    # print('Degree 1', len([n for n in g if nx.degree(g, n) == 1]))
    # print('Degree 2', len([n for n in g if nx.degree(g, n) == 2]))
    # print()

    # Final pass:
    if g.number_of_nodes() == 2:
        print('HERE')
        [u,v] = list(g.nodes())
        merge(g, u, v, next_name = next_node, distances=distances, names=names)
        mseq.append([u, v, next_node])
        next_node += 1

    sorted_nodes = list(sorted(g))
    pre_merge_map = {sorted_nodes[i]:(i+1) for i in range(len(sorted_nodes))}
    nx.relabel_nodes(g, pre_merge_map, copy=False)
    g.graph['root'] = pre_merge_map[g.graph['root']]

    return g, mseq, pre_merge_map



def merge_sequences(pre_mseq, mseq, pre_merge_map):
    N = pre_mseq[0][2]-1
    Nf = pre_mseq[-1][2]
    M = mseq[0][2]-1
    reverse_map = {pre_merge_map[k]:k for k in pre_merge_map}
    for i in range(len(mseq)):
        for j in range(3):
            if mseq[i][j] <= M:
                mseq[i][j] = reverse_map[mseq[i][j]]
            else:
                mseq[i][j] = mseq[i][j] - M + Nf
    return pre_mseq + mseq


def pseudotop(G):
    g = G.copy()
    def valid_merge(n):
        # True if degree 2 with degree 2 neighbors
        if nx.degree(g, n) == 2 and all([nx.degree(g, neigh) == 2 for neigh in nx.neighbors(g, n)]):
            return True
    q = [e for e in g.edges() if valid_merge(e[0]) or valid_merge(e[1])]
    while q != []:
        (u, v) = q.pop()
        if not valid_merge(v):
            next_name = v
            props = g.nodes[v]
        else:
            next_name = u
            props = g.nodes[u]
        merge(g, u, v, next_name=next_name)
        for k in props:
            g.nodes[next_name][k] = props[k]
    # Rename
    sorted_nodes = sorted(g)
    nmap = {sorted_nodes[i]:(i+1) for i in range(len(sorted_nodes))}
    nx.relabel_nodes(g, nmap, copy=False)
    g.graph['nodemap'] = nmap
    return g


import random

def random_merge_seq(G, adjacent=False):
    """ Generate a random merge sequence - adjacent indicates whether two merged nodes should have an edge between them. """
    mseq = []
    next_node = len(G)+1
    g = G.copy()
    #choices = list(G.edges())
    while len(mseq) != len(G)-1:
        (u,v) = list(g.edges())[random.randint(0, g.number_of_edges()-1)]#choices.pop(random.randint(0, len(choices)))
        mseq.append([u,v,next_node])
        g.remove_edge(u,v)

        nedges = []
        for neigh in nx.neighbors(g, u):
            if neigh != v:
                num = g.number_of_edges(u, neigh)
                nedges += [(next_node, neigh)]*num
        for neigh in nx.neighbors(g, v):
            if neigh != u:
                num = g.number_of_edges(v, neigh)
                nedges += [(next_node, neigh)]*num
        g.remove_node(u)
        g.remove_node(v)
        g.add_edges_from(nedges)

        next_node += 1

    return mseq

