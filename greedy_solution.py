import networkx as nx
import partitiongraph as pg
import partitions as prt 
import dpartitions as dp
import numpy as np
from queue import PriorityQueue, Queue 
from itertools import combinations
import os


def fix_costs(g, cname='cost'):
    c = g.graph[cname]
    #c[g.graph['root']] = {g.graph['root']:0.0 for k in c[g.graph['root']]}
    c[g.graph['root']] = {k:0.0 for k in c[g.graph['root']]}
    deg12 = [n for n in g if nx.degree(g, n) <= 2]
    for n in deg12:
        c[n] = {k:0.0 for k in c[n]}
    # junctions = [n for n in g if nx.degree(g, n) > 2]
    # for n in junctions:
    #     for p in g.nodes[n]['partition_candidates']:
    #         if len(p.partition) == 1:
    #             enc = dp.get_encoding(p)#g.nodes[n]['parent_lists'][dp.get_encoding(p)]
    #             matches = [k for k in c[n] if k[1:].startswith(enc)]
    #             for k in matches:
    #                 #print(k)
    #                 c[n][k] = 0.0


def node_dist(g, u, v):
    u_vec = np.array([g.nodes[u]['x'], g.nodes[u]['y'], g.nodes[u]['z']])
    v_vec = np.array([g.nodes[v]['x'], g.nodes[v]['y'], g.nodes[v]['z']])
    return np.linalg.norm(u_vec - v_vec)


def topological_sort(g, r, dist=node_dist):
    # Initialize
    for n in g:
        g.nodes[n]['greedy_parents'] = np.zeros((nx.degree(g, n),))-1

    # Top sort
    distances = np.zeros((g.number_of_nodes(),))
    for (u, v) in nx.bfs_edges(g, r):
        ename = g.edges[(u,v,0)]['name']
        vdict = g.nodes[v]['parent_dict']
        g.nodes[v]['greedy_parents'][vdict[v]] = 1
        if v != r and distances[v-1]==0:
            distances[v-1] = distances[u-1] + dist(g, u, v)
        elif v != r:
            distances[v-1] = min(distances[u-1] + dist(g, u, v), distances[v-1]) #??? TODO: BFS with queue based on eucledian distance
    return distances.argsort() + 1


def dist_spanning_tree(g, r, dist=node_dist):
    q = PriorityQueue()
    q.put((0, r))
    #parents = [-np.ones((nx.degree(g, n),)) for n in g]
    for n in g: g.nodes[n]['greedy_parents'] = np.ones((nx.degree(g, n),)) 
    seen = [False]*g.number_of_nodes()
    #dists = [0]*g.number_of_nodes()
    for n in g: g.nodes[n]['dist'] = 0 
    while not q.empty():
        (c, n) = q.get()
        seen[n-1] = True
        pdict = g.nodes[n]['parent_dict']
        for u in nx.neighbors(g, n):
            if not seen[u-1]:
                #q.put((dist(g, n, u), u))
                #dists[u-1] = min([g.edges[e]['len'] for e in g.edges(n, keys=True) if u in e[:2]]) + dists[n-1]
                #q.put((dists[u-1], u))
                g.nodes[u]['dist'] = min([g.edges[e]['len'] for e in g.edges(n, keys=True) if u in e[:2]]) + g.nodes[n]['dist']
                q.put((g.nodes[u]['dist'], u))
                g.nodes[n]['greedy_parents'][pdict[g.edges[(n,u,0)]['name']]] = -1
                #parents[n-1][pdict[g.edges[(n,u,0)]['name']]] = +1
    #return parents


def dist_spanning_tree_3(g, r):
    targets = nx.shortest_path_length(g, g.graph['root'], weight='len')
    #paths = nx.shortest_path(g, g.graph['root'], weight='dist')
    nx.set_node_attributes(g, targets, 'dist')
    for n in g: g.nodes[n]['greedy_parents'] = np.ones((nx.degree(g, n),)) 
    for (u, v, k) in g.edges(keys=True):
        n = min((g.nodes[u]['dist'], u), (g.nodes[v]['dist'], v))[1]
        pdict = g.nodes[n]['parent_dict']
        g.nodes[n]['greedy_parents'][pdict[g.edges[(u,v,k)]['name']]] = -1
    



def project_edge(g, e):
    """ Given an edge e, project it towards the edge direction until a junction is found """
    (u, v, k) = e
    edges = [(u, v, k)]
    while nx.degree(g, v) == 2:
        #(u, v, k) = [e for e in g.edges(v, keys=True) if e != (u,v,k) and e != (v,u,k)][0]
        v, u = [n for n in nx.neighbors(g, v) if n != u][0], v
        edges.append((u, v, 0))
    return edges


def fix_dirs(g):
    """ Given parents for each node, make it so that degree 2 vertices have no parents """
    q = Queue()
    #map(q.put, [n for n in g if nx.degree(g,n)==2 and np.all(g.nodes[n]['parents']==1)])
    deg_2_junctions = [q.put(n) for n in g if nx.degree(g,n)==2 and np.all(g.nodes[n]['greedy_parents']==1)]
    #print([n for n in g if nx.degree(g,n)==2 and np.all(g.nodes[n]['parents']==1)])
    while not q.empty():
        n = q.get()
        [n1, n2] = nx.neighbors(g, n)
        edges_1 = project_edge(g, (n, n1, 0))
        edges_2 = project_edge(g, (n, n2, 0))
        
        dist_1 = g.nodes[edges_1[-1][1]]['dist']#sum([g.edges[e]['len'] for e in edges_1]) + g.nodes[n2]['dist']
        dist_2 = g.nodes[edges_2[-1][1]]['dist']#sum([g.edges[e]['len'] for e in edges_2]) + g.nodes[n1]['dist']
        
        option = max((dist_1, edges_1), (dist_2, edges_2)) # The max distance one loses priority
        
        for (u, v, k) in option[1]:
            ename = g.edges[(u,v,k)]['name']
            g.nodes[u]['greedy_parents'][g.nodes[u]['parent_dict'][ename]] = -1
            g.nodes[v]['greedy_parents'][g.nodes[v]['parent_dict'][ename]] = +1


def dist_spanning_tree_2(g, r): # Without node metadata

    q = PriorityQueue()
    seen = [False]*g.number_of_nodes()
    q.put((0, r))
    seen[r-1] = True
    dirs = []
    dists = [0]*g.number_of_nodes()

    while not q.empty():
        (c, n) = q.get()
        for u in nx.neighbors(g, n):
            if not seen[u-1]:
                dists[u-1] = min([g.edges[e]['len'] for e in g.edges(n, keys=True) if u in e[:2]]) + dists[n-1]
                q.put((dists[u-1], u))
                seen[u-1] = True # Incorrect
                dirs.append((n, u)) # n -> u
    return dirs



def all_partitions_fixed_direction(edges, parents):
    """ Given a parentlist, generate all possible partitions consistent with those parentlists """
    edges = np.array(edges)
    (parent_inds,) = np.where(parents==1)
    (child_inds,) = np.where(parents!=1)
    child_edges = [edges[i] for i in child_inds.tolist()]
    parent_edges = [edges[p_ind] for p_ind in parent_inds]
    if child_edges == []:
        return [[[e] for e in parent_edges]]
    all_undirected = [p.partition for p in prt.Partition.all_partitions([edges[i] for i in child_inds.tolist()]) if len(p.partition) <= len(parent_edges)]
    #print('Child partitions', all_undirected)
    #print('Parents', parents, np.where(parents==1))

    has_parent = dict()
    
    for p_ind in parent_inds:
        p_edge = edges[p_ind]
        all_undirected_pind = [part + [[p_edge]] for part in all_undirected if len(part + [[p_edge]]) <= len(parent_edges)]
        #print('\tAdding', p_edge)
        for i in range(len(all_undirected)):
            pcand = all_undirected[i]
            if len(pcand) > len(parent_edges): continue
            for j in range(len(all_undirected[i])):
                pset = pcand[j]
                if (set(parent_edges) & set(pcand[j])) == set():
                    new_cand = pcand[:j] + [list(pcand[j]) + [p_edge]] + pcand[(j+1):]
                    all_undirected_pind.append(new_cand)
        all_undirected = all_undirected_pind
    return all_undirected


def directed_optimal(G):
    """ Given a set of parents for each vertex, set the optimal partition """
    #gs.dist_spanning_tree(G, G.graph['root'])
    costs = G.graph['cost']
    for n in G.nodes():
        # for each node, compute optimal partition, given directions
        if nx.degree(G, n) <= 2:
            continue 
        pdict = G.nodes[n]['parent_dict']
        edges = [0]*G.nodes[n]['greedy_parents'].size
        parents = G.nodes[n]['greedy_parents']
        for (ename, i) in pdict.items():
            edges[i] = int(ename)

        # Compute all partition candidates consistent with parents
        partitions = gs.all_partitions_fixed_direction(edges, parents)
        partitions = [prt.Partition(part) for part in partitions]
        
        # Get candidate with minimal cost
        mincost = (float('inf'), -1)
        for partition in partitions:
            enc = dp.get_encoding(partition, (pdict, parents), np_array=True)
            if costs[n][enc] < mincost[0]:
                mincost = (costs[n][enc], partition)
        g.nodes[n]['greedy_partition'] = mincost[1]



def greedy_solution(G):
    for e in G.edges(keys=True):
        if 'len' not in G.edges[e]:
            vec1 = np.array([G.nodes[e[0]]['x'], G.nodes[e[0]]['y'], G.nodes[e[0]]['z']])
            vec2 = np.array([G.nodes[e[1]]['x'], G.nodes[e[1]]['y'], G.nodes[e[1]]['z']])
            G.edges[e]['len'] = np.linalg.norm(vec1 - vec2)

    dist_spanning_tree_3(G, G.graph['root'])
    fix_dirs(G)
    
    for n in G:
        pdict = G.nodes[n]['parent_dict']
        pairs = [(G.edges[e]['name'], G.nodes[n]['greedy_parents'][pdict[G.edges[e]['name']]]) for e in G.edges(n, keys=True)]
        edges = [pair[0] for pair in pairs]
        parents = [pair[1] for pair in pairs]
        if G.graph['root'] != n:
            parts = all_partitions_fixed_direction(edges, np.array(parents))
            cost_parts = []
            #print(n)
            for part in parts:
                enc = dp.get_encoding(prt.Partition(part), (pdict, G.nodes[n]['greedy_parents']), np_array=True)
                cost_parts.append((G.graph['cost'][n][enc], part))
            G.nodes[n]['greedy_partition'] = prt.Partition(min(cost_parts)[1])
        else:
            child_part = [edges[i] for i in range(len(edges)) if parents[i] == -1]
            part = [[edges[i]] for i in range(len(edges)) if parents[i] == +1] + [child_part]
            G.nodes[n]['greedy_partition'] = prt.Partition(part)


def fix_triangle_partitions(g, n, cands):
    if nx.degree(g, n) > 2:
        deg_2_neighbors = [neigh for neigh in nx.neighbors(g, n) if nx.degree(g, neigh) == 2]
        triangle_edges = [(g.edges[(n,u,0)]['name'], g.edges[(n,v,0)]['name']) \
            for (u, v) in combinations(deg_2_neighbors, 2) if u in g[v]]
        good_candidates = []
        for i in range(len(cands)):
            cand = cands[i]
            # For each candidate, we want to remove those that have
            # both v and u in one pset
            good_candidate = True
            for (u, v) in triangle_edges:
                for pset in cand:
                    if u in pset and v in pset:
                        good_candidate = False
                        break
                    elif u in pset or v in pset:
                        break
                if not good_candidate:
                    break
            if good_candidate:
                good_candidates.append(cand)
        return good_candidates
    else:
        return cands

def get_time():
    return getattr(os.times(), 'user')

def get_greedy(G):

    timings = []

    t0 = get_time()  # Time step 1: obtain shortest paths tree
    for e in G.edges(keys=True):
        if 'len' not in G.edges[e]:
            vec1 = np.array([G.nodes[e[0]]['x'], G.nodes[e[0]]['y'], G.nodes[e[0]]['z']])
            vec2 = np.array([G.nodes[e[1]]['x'], G.nodes[e[1]]['y'], G.nodes[e[1]]['z']])
            G.edges[e]['len'] = np.linalg.norm(vec1 - vec2)

    timings.append(get_time() - t0)
    t0 = get_time()  # Time step 2: Orient directions
    dist_spanning_tree_3(G, G.graph['root'])
    fix_dirs(G)
    timings.append(get_time() - t0)
    t0 = get_time()  # Time step 3: find partitions
    
    ocost = G.graph['cost']
    qcost = {n:{enc:(True, ocost[n][enc]) for enc in ocost[n]} for n in ocost}
    
    for n in G:
        pdict = G.nodes[n]['parent_dict']
        pairs = [(G.edges[e]['name'], G.nodes[n]['greedy_parents'][pdict[G.edges[e]['name']]]) for e in G.edges(n, keys=True)]
        edges = [pair[0] for pair in pairs]
        parents = [pair[1] for pair in pairs]
        if G.graph['root'] != n:
            parts = all_partitions_fixed_direction(edges, np.array(parents))
            parts = fix_triangle_partitions(G, n, parts)
            cost_parts = []
            #print(n)
            for part in parts:
                enc = dp.get_encoding(prt.Partition(part), (pdict, G.nodes[n]['greedy_parents']), np_array=True)
                cost_parts.append((ocost[n][enc], part, enc))
            (chosen_cost, chosen_part, chosen_enc) = min(cost_parts)
            G.nodes[n]['greedy_partition'] = prt.Partition(chosen_part)
            qcost[n][chosen_enc] = (False, chosen_cost)
        else:
            child_part = [edges[i] for i in range(len(edges)) if parents[i] == -1]
            part = [[edges[i]] for i in range(len(edges)) if parents[i] == +1] + [child_part]
            G.nodes[n]['greedy_partition'] = prt.Partition(part) 
            chosen_enc = dp.get_encoding(prt.Partition(part), (pdict, G.nodes[n]['greedy_parents']), np_array=True)
            if chosen_enc in ocost[n]:
                qcost[n][chosen_enc] = (False, ocost[n][chosen_enc])
            else:
                qcost[n][chosen_enc] = (False, 0)
    timings.append(get_time() - t0)
    G.graph['priority cost'] = qcost
    return timings