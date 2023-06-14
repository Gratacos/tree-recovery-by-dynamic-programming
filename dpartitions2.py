import numpy as np
from typing import List, Dict
import partitions as prt

# The format of the directions may change, as such, we make a type for it

Direction = np.array

def generate_intersection_map(P1: List[prt.Partition], D1map, \
                              P2: List[prt.Partition], D2map):
    """
    Input:  P1  - Undirected partitions
            D1  - Direction definitions for P1
                    (Dicitonary of correspondences, Bitmap)
            P2  - Undirected partitions
            D2  - Direction definitions for P2
                    (Dictionary of correspondences, Bitmap)
    Output: Dp1 - Direction definitions for non-border-edges in P1
            Dp2 - Direction definitions for non-border-edges in P2
            B   - Direction definitions for border edges in both with
                    mapping
    """

    #(D1map, D1bits) = D1
    #(D2map, D2bits) = D2

    # Get edges T(# total edges) from a sample partition
    # e1 = list(P1.fullset)
    # e2 = list(P2.fullset)
    e1 = list(D1map.keys())
    e2 = list(D2map.keys())
    edges = e1 + e2
    # print("e1:", e1 )
    # print("e2:", e2)

    # Iterate once over both edges to know which are border edges
    # Result:   is_border[edge] = is edge a border edge?
    #           B_map[edge] = new index  // is this necessary? 
    # T(# edges + # in edges)
    B_map = dict()      # edge -> index
    is_border = []      # index -> bool
    index = 0
    split = float("inf")
    corr = dict()
    for i in range(len(edges)):
        e = edges[i]
        if e in B_map:  # Not a set -> should be O(1)
            is_border[B_map[e]] = False
        else:
            B_map[e] = index
            is_border.append(True)
            #print(i, ">=", len(e1), "for", e)
            if i >= len(e1):
                corr[e] = 1
            else:
                corr[e] = 0
            index += 1
    #print(corr)

    # Extract border/non-border edges 
    # Result:   border_edges = edge (if edge is border)
    #           in_edges = edge (if edge is in-edge)
    # T(# edges)
    in_edges = []
    border_map_1 = []
    border_map_2 = []
    border_edges_1 = []
    border_edges_2 = []
    ind = 0
    
    for e in B_map:
        #print("-")
        if not is_border[B_map[e]]:
            in_edges.append(e)  # in_edges: not border
            ind += 1
        elif not corr[e]:
            border_map_1.append(D1map[e])
            border_edges_1.append(e)
        else:
            border_map_2.append(D2map[e])
            border_edges_2.append(e)

    # Create bitmaps that correspond to each other
    # Result:   Dp1 = bit matrix of in-edges 1
    #           Dp2 = bit matrix of in-edges 2
    # T(# nonborder-edges)
    D1ord = [-1] * len(in_edges)
    D2ord = [-1] * len(in_edges)
    for i in range(len(in_edges)):
        D1ord[i] = D1map[in_edges[i]]
        D2ord[i] = D2map[in_edges[i]]
    
    
    # Dp1 = D1bits[D1ord, :]
    # Dp2 = D2bits[D2ord, :]

    obj = {
        "ords" : (D1ord, D2ord), # Multiply: Dp1[D1ord, :].T @ Dp2[D2ord, :]
        "border_edges" :  (border_map_1, border_edges_1, border_map_2, border_edges_2),
        "num_edges" : len(in_edges),
        "border" : border_edges_1 + border_edges_2,
        "non-border" : in_edges
    }

    return obj
    #return Dp1, Dp2, in_edges, (border_map_1, border_edges_1, border_map_2, border_edges_2)   #(Bbits, Bmap)



def direction_join( join_obj, direction_obj ):
    """
    Input:  List[(encoding, bit1, bit2)]
            direction_obj { "ords": orders for bitmap,
                            "border_edges" : information for border edges}
    Output: Dictionary { "encoding" : bitmap }
            Dictionary { edge : index }
    """
    # Extract information from input
    (D1ord, D2ord) = direction_obj["ords"]
    (border_map_1, border_edges_1, border_map_2, border_edges_2) = direction_obj["border_edges"]
    num_in = direction_obj["num_edges"]

    outbits_dict = {}

    for (partition_encoding, bitmap1, bitmap2) in join_obj:
        R = bitmap1[D1ord, :].T @ bitmap2[D2ord, :]
        (p1chosen, p2chosen) = np.where(R == -num_in)

        chosenD1 = bitmap1[:, p1chosen]
        chosenD2 = bitmap2[:, p2chosen]

        outbits1 = chosenD1[border_map_1, :]
        outbits2 = chosenD2[border_map_2, :]
        outbits = np.vstack((outbits1, outbits2))

        outbits_dict[partition_encoding] = outbits

    bmap = {border_edges_1[i] : i for i in range(len(border_edges_1))}
    bmap = {**bmap, **{border_edges_2[i] : i + len(border_edges_1) for i in range(len(border_edges_2))}}

    return outbits_dict, bmap

import dpartitions as dp
import time
def construct_matrices(G, u, v, new_node, back_track, cost, new_part_candidates):
    pcands1 = G.nodes[u]['partition_candidates']
    pcands2 = G.nodes[v]['partition_candidates']
    dcands1 = G.nodes[u]['parent_lists']
    dcands2 = G.nodes[v]['parent_lists']
    bdict1 = G.nodes[u]['parent_dict']
    bdict2 = G.nodes[v]['parent_dict']

    all_edges = pcands1[0].fullset | pcands2[0].fullset
    intern_edges = pcands1[0].fullset & pcands2[0].fullset
    num_in = len(intern_edges)
    d_obj={}
    d_obj['border-dict'] = {e:True for e in all_edges}
    d_obj["non-border"] = intern_edges
    for e in intern_edges: d_obj['border-dict'][e] = False
    border_map_1 = [bdict1[e] for e in pcands1[0].fullset if d_obj['border-dict'][e]] # Locations of border edges
    border_map_2 = [bdict2[e] for e in pcands2[0].fullset if d_obj['border-dict'][e]] # Locations of border edges
    d1ord = [bdict1[e] for e in pcands1[0].fullset if not d_obj['border-dict'][e]] # locations of non-border edges
    d2ord = [bdict2[e] for e in pcands2[0].fullset if not d_obj['border-dict'][e]] # locations of non-border edges

    d1ord_inv = {border_map_1[i]:i for i in range(len(border_map_1))}
    d2ord_inv = {border_map_2[i]:i for i in range(len(border_map_2))}

    #bmap = {border_edges_1[i] : i for i in range(len(border_edges_1))}
    bmap = {**{e:d1ord_inv[bdict1[e]] for e in bdict1 if d_obj['border-dict'][e]}, **{e:(d2ord_inv[bdict2[e]]+len(border_map_1)) for e in bdict2 if d_obj['border-dict'][e]}}

    seen = dict()
    new_cost = cost.copy()
    new_cost[new_node] = dict()  # Restart this node
    #new_part_candidates = BeamList(params)
    # Join matrices for node 1
    j1 = [0]*len(dcands1)
    position_map_1 = {}
    it = 0

    #cprev = 0
    current_pos = 0
    for (k, bmat) in dcands1.items():
        bmat_in = bmat[d1ord, :]
        (r,c) = bmat_in.shape
        j1[it] = bmat_in
        position_map_1[k] = (current_pos, current_pos + c)
        it += 1
        
        current_pos = current_pos+ c
    
    # Join matrices for node 2
    j2 = [0]*len(dcands2)
    position_map_2 = {}
    it = 0
    #cprev = 0 
    current_pos = 0
    for (k, bmat) in dcands2.items():
        bmat_in = bmat[d2ord, :]
        (r,c) = bmat_in.shape
        j2[it] = bmat_in
        position_map_2[k] = (current_pos, current_pos + c)
        it += 1
        
        current_pos = current_pos + c

    
    # Stack matrices and multiply them
    M1 = np.hstack(j1)
    M2 = np.hstack(j2)
    R_total = M1.T @ M2
    (p1chosen_total, p2chosen_total) = np.where(R_total == -num_in)


    # Assumption - p1chosen is sorted - is wrong.
    # chosen_map_1 = {}
    # for (k, (s, e)) in position_map_1:
    #     chosen_map_1[k] = p1chosen_total[s:e]

    # chosen_map_2 = {}
    # for (k, (s, e)) in position_map_2:
    #     chosen_map_2[k] = p2chosen_total[s:e]

    # Iterate over partition candidates
    t0 = 0
    t1 = 0
    t2 =0
    t3 =0
    t4 = 0
    for partition1 in pcands1:
        uenc1 = dp.get_encoding(partition1)
        (s1, e1) = position_map_1[uenc1]
        bitmap1 = dcands1[uenc1]
        cond1 = np.logical_and(p1chosen_total >= s1, p1chosen_total < e1)
        for partition2 in pcands2:
            # Directed check

            t0 = time.time()

            uenc2 = dp.get_encoding(partition2)
            (s2, e2) = position_map_2[uenc2]
            bitmap2 = dcands2[uenc2]
            # R = R_total[s1:e1, s2:e2]

            t1 += time.time() - t0
            t0 = time.time()

            # (p1chosen, p2chosen) = np.where(R == -num_in)  # each of these indices corresponds to a cost
            cond = np.logical_and(cond1, np.logical_and(p2chosen_total >= s2, p2chosen_total < e2))
            p1chosen = p1chosen_total[cond] - s1
            p2chosen = p2chosen_total[cond] - s2
            chosenD1 = bitmap1[:, p1chosen]#.copy()
            chosenD2 = bitmap2[:, p2chosen]#.copy()
            outbits1 = chosenD1[border_map_1, :]
            outbits2 = chosenD2[border_map_2, :]
            # try:
            outbits = np.vstack((outbits1, outbits2))
            # except:
            #     print('Exception')

            t2 += time.time() - t0
            t0 = time.time()
            
            if outbits.size == 0: 
                # ujoin_test, ujoin2, t2 = dp.ujoin(partition1, partition2, d_obj)
                # if ujoin_test is not None:
                #     print('Here')
                #     pass
                continue

            # Undirected check
            ujoin, ujoin2, t2 = dp.ujoin(partition1, partition2, d_obj)


            t3 += time.time() - t0
            t0 = time.time()


            if ujoin is None: continue

            uenc_join = dp.get_encoding(ujoin) 
            (R, C) = outbits.shape
            for c in range(C):
                parents1 = (bdict1, chosenD1[:, c])
                parents2 = (bdict2, chosenD2[:, c])
                join = (ujoin, outbits[:,c])

                join_ind = "[" + uenc_join + ", " +  dp.get_encoding(ujoin, parents=(bmap, outbits[:,c]), sep=False, np_array=True, und=False) + "]"
                p1_ind = "[" + uenc1 + ", " + dp.get_encoding(partition1, parents1, sep=False, np_array=True, und=False) + "]"
                p2_ind = "[" + uenc2 + ", " + dp.get_encoding(partition2, parents2, sep=False, np_array=True, und=False) + "]"
        
                if join_ind not in seen:
                    seen[join_ind] = 1
                    back_track[-1][new_node][join_ind] = (-1, -1)

                if new_node not in new_cost:
                    new_cost[new_node] = {join_ind:float('inf')}
                elif join_ind not in new_cost[new_node]:
                    new_cost[new_node][join_ind] = float("inf")

                if new_cost[new_node][join_ind] > (cost[u][p1_ind]+cost[v][p2_ind]):
                    new_cost[new_node][join_ind] = (cost[u][p1_ind]+cost[v][p2_ind])
                    back_track[-1][new_node][join_ind] = {u:(partition1, parents1, bdict1), v:(partition2, parents2, bdict2)}
                    new_part_candidates.add(new_cost[new_node][join_ind], join, join_ind)
    
            t4 += time.time() - t0
            t0 = time.time()

    return back_track, new_cost, d_obj, bmap, border_map_1, border_map_2, (t1, t2, t3, t4)








