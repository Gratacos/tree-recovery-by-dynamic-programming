import networkx as nx
import numpy as np
from typing import List, Tuple
from functools import reduce
import scipy.sparse.csgraph as csg

from disjoint_set import DisjointSet # DS Code

#
# Class for handling partitions. Includes:
#   * static method generate all possible partitions of a given set. 
#   * a way to uniquely encode each partition
#   * join operator, given two partitions
#

PartitionT = List[List[int]]
#DPartitionT = List[(int, List[int])]

class Partition:
    """ Class for handling partition. """
    def __init__(self, partition: PartitionT, cost: float=0, max_inner_trees=0, root=-1):
        self.partition = partition
        
        # Change to set
        for i in range(len(self.partition)):
            self.partition[i] = set(self.partition[i])
        
        self.fullset = reduce(set.union, partition, set())

        self.encoding = None
        self.encoding = self.get_encoding()
        self.cost = cost
        self.num_inner_trees = 0
        self.max_inner_trees = 0
        self.root = root
    

    def get_encoding(self) -> int:
        """ 
        Returns a unique encoding for a partition 
        """
        if self.encoding is None:  # Avoid calculating twice
            items = self.fullset
            partition = self.partition

            sorted_items = sorted(items)
            self.item_map = dict()
            for i in range(len(sorted_items)): 
                self.item_map[sorted_items[i]] = i

            encoding_list = [0]*len(items)
            self.encoding_map = dict()  # useful data structure for user. Gets item, returns partition ind
            for i in range(len(partition)):
                pset = partition[i]
                for item in pset:
                    item_ind = self.item_map[item]
                    encoding_list[item_ind] = i
                    self.encoding_map[item] = i  # For future use
            
            base = len(items)
            encoding = 0
            part_ind = 0
            pmap = dict()

            for i in reversed(range(len(encoding_list))):
                if encoding_list[i] not in pmap:
                    pmap[encoding_list[i]] = part_ind
                    part_ind += 1
                encoding += base**i * pmap[encoding_list[i]]

        else:
            encoding = self.encoding

                    

        return encoding


    def join_um1(self,part2,djoin,root_edges=None):
        # p1, p2 = self.partition, part2.partition
        p1, p2 = [list(pset) for pset in self.partition], [list(pset) for pset in part2.partition]
        cutedge = set(djoin['border'])
        
        parent_diction = {e:e for pset in p2 for e in pset}
        for pset in p1:
            for e in pset:
                parent_diction[e]=pset[-1]
        ds = DisjointSet(parent_diction)
        [[ds.union(pset[i], pset[i+1]) for i in range(len(pset)-1)] for pset in p2]
        
        partition=[]
        
        # Root/connected check
        if root_edges is None:
            for s in ds.itersets():
                intersection=s & cutedge
                out_edges = [djoin['bmap'][e] == -1 for e in intersection]
                if len(intersection)==0 or (all(out_edges)):
                    return None
                else:
                    partition.append(intersection)
        else:
            root_found = False # To avoid checking for multiple roots
            for s in ds.itersets():
                intersection=s & cutedge
                out_edges = [djoin['bmap'][e] == -1 for e in intersection]
                if not root_found:
                    root_inter = root_edges & s
                    if len(root_inter) != 0:
                        root_found =True
                else:
                    root_inter = {}
                if len(intersection)==0 or (all(out_edges) and len(root_inter) == 0):
                    return None
                else:
                    partition.append(intersection)
                    

        return Partition(partition)

    
    def join_uf(self, partition2, djoin, root=None):
        p1 = self.partition
        p2 = partition2.partition
        cutedge = djoin['border-dict']
        interedge = djoin['non-border']
        partitions = p1 + p2
        cutedges = [e for e in cutedge if cutedge[e]]
        # edge2psets = {p1[i][j]:[i, 0] for i in range(len(p1)) for j in range(len(p1[i]))   }
        edge2psets = {e:[i, 0] for i in range(len(p1)) for e in p1[i]   }
        for i in range(len(p2)):
            for e in p2[i]:
                if e in edge2psets:
                    edge2psets[e][1] = i
            # for j in range(len(p2[i])):
            #     if p2[i][j] in edge2psets:
            #         edge2psets[p2[i][j]][1] = i

        outer_sets = [pset for pset in partitions if all(cutedge[e] for e in pset)]

        

        ds = DisjointSet()
        
        [ds.union(edge2psets[e][0], edge2psets[e][1]+len(p1)) for e in interedge]
        
        partition = []
        if root is None:
            for psets in ds.itersets():
                partition_set = [e for pind in psets for e in partitions[pind] if cutedge[e]]

                if partition_set == [] or all([djoin['bmap'][e] == -1 for e in partition_set]):
                    return None
                partition.append(partition_set)
        else:
            for psets in ds.itersets():
                partition_set = {e for pind in psets for e in partitions[pind] if cutedge[e]}
                if partition_set == {} or (all([djoin['bmap'][e] == -1 for e in partition_set]) and root & partition_set == {}):
                    return None
                partition.append(partition_set)

            
        return Partition(partition + outer_sets)


    def join(self, partition2, djoin):
        p1 = self.partition.copy()
        p2 = partition2.partition#.copy()
        cutedge = djoin['border-dict']
 
        edge2pset1 = {e:i for i in range(len(p1)) for e in p1[i]}
        p1_inner = [[e for e in pset if not cutedge[e]] for pset in p1]

        delete = set()
        for pset in p2:
            S = {e for e in pset if not cutedge[e]}
            inner_sets_s = {edge2pset1[e] for e in S}
            if len(inner_sets_s) != len(S):
                return None

            inner_sets = list(inner_sets_s)

            A = []
            for it in inner_sets:
                A += p1[it]

            E = (set(A)|pset)-S
            if len(E) == 0:
                return None

            A_inner = []
            for it in inner_sets:
                A_inner += p1_inner[it]

            E_inner = set(A)-S

            # Inner component
            delete |= inner_sets_s

            p1.append(E)
            p1_inner.append(E_inner)
        
            ind = len(p1_inner)-1
            for e in E_inner:
                edge2pset1[e] = ind

        return Partition([p1[i] for i in range(len(p1)) if i not in delete])

    def join_list(self, partition2, djoin):
        # TODO - p2 is wrong
        p1 = self.partition.copy()
        p2 = [{e + len(self.fullset) for e in pset} for pset in partition2.partition]#.copy()
        cutedge = djoin['border-list']
 
        edge2pset1 = [0]*len(cutedge)
        for i in range(len(p1)):
            for e in p1[i]:
                edge2pset1[e] = i
        p1_inner = [[e for e in pset if not cutedge[e]] for pset in p1]

        delete = set()
        for pset in p2:
            S = {e for e in pset if not cutedge[e]}
            inner_sets_s = {edge2pset1[e] for e in S}
            if len(inner_sets_s) != len(S):
                return None
            
            inner_sets = list(inner_sets_s)

            A = []
            for it in inner_sets:
                A += p1[it]

            E = (set(A)|pset)-S
            if len(E) == 0:
                return None

            A_inner = []
            for it in inner_sets:
                A_inner += p1_inner[it]

            E_inner = set(A)-S

            # Inner component
            delete |= inner_sets_s

            p1.append(E)
            p1_inner.append(E_inner)
        
            ind = len(p1_inner)-1
            for e in E_inner:
                edge2pset1[e] = ind

        return Partition([p1[i] for i in range(len(p1)) if i not in delete])


    def connected(self, partition2):
        p1 = self.partition.copy()
        p2 = partition2.partition.copy()
        # cutedge = djoin['border-dict']
        # interedge = djoin['non-border']
        edge2pset1 = {e:i for i in range(len(p1)) for e in p1[i]}
        p1_inner = [[e for e in pset] for pset in p1]

        # p1_cut_sets = [pset for pset in p1 if all([cutedge[e] for e in pset])]
        # start = len(p1)
        delete =set()
        for pset in p2:
            S = {e for e in pset}
            inner_sets_s = {edge2pset1[e] for e in S}
            if len(inner_sets_s) != len(S):
                #assert self.connected_old(partition2) == False, str(p1) + ' ' + str(p2) + ' connected check wrong'
                return False

            
        
            inner_sets = list(inner_sets_s)
            # print('Inner sets:', inner_sets)
        
            # A = [p1[it] for it in inner_sets]
            A = []
            for it in inner_sets:
                A += p1[it]
            # E = (set().union(*A)|pset)-S
            E = (set(A)|pset)-S
    
            # A_inner = [p1_inner[it] for it in inner_sets]
            # E_inner = set().union(*A_inner)-S
            A_inner = []
            for it in inner_sets:
                A_inner += p1_inner[it]
            # E = (set().union(*A)|pset)-S
            E_inner = set(A)-S

            delete |= inner_sets_s

        
            p1.append(E)
            p1_inner.append(E_inner)
        
            ind = len(p1_inner)-1
            for e in E_inner:
                edge2pset1[e] = ind
                # if loops[e] > 2:
                #     loops[e]+=1
                # else:
                #     assert self.connected_old(partition2) == False, str(p1) + ' ' + str(p2) + ' connected check wrong'
                #     return False
        comps = len([p1[i] for i in range(len(p1))]) - len(delete)

        #print([p1[i] for i in range(len(p1)) if i not in delete2])
        #assert ( (comps == 1) == self.connected_old(partition2)), str(p1) + ' ' + str(p2) + ' connected check wrong'
        return comps == 1



    #def join2(p1, p2, all_edges, cutedge, interedge, bmap, root=None):
    def join_uf(self, partition2, djoin, root=None):
        p1 = self.partition
        p2 = partition2.partition
        cutedge = djoin['border-dict']
        interedge = djoin['non-border']
        #bmap = djoin['bmap']
        # Create map: Edge -> (Set location in p1, Set location in p2)
        edge2psets = {e:[i, 0] for i in range(len(p1)) for e in p1[i]   }
        for i in range(len(p2)):
            for e in p2[i]:
                if e in edge2psets:
                    edge2psets[e][1] = i
                    
        #outer_sets = [pset for pset in partitions if all(cutedge[e] for e in pset)]

        # Perform Union Find
        ds = UF(len(p1) + len(p2))
        [ds.union(edge2psets[e][0], edge2psets[e][1]+len(p1)) for e in interedge]
        
        # partitions = [[] for i in range(max(ds._id)+1)]

        # ids: Old partition index -> New partition index
        ids = ds.itersets()
        ids2 = ds._id


        partitions = [[] for i in range(max(ids)+1)]
        inner_trees = [set() for i in range(max(ids)+1)]
        
        #print(ds._id, ':', ds._rank)
        for i in range(len(p1)):
            # # EXPERIMENTAL
            # print(inner_trees[ids[i]])
            # print('\t', {e for e in p1[i] if not cutedge[e]})
            # if len(inner_trees[ids[i]] & {e for e in p1[i] if not cutedge[e]}) <= 1:
            #     inner_trees[ids[i]] = inner_trees[ids[i]] | set(p1[i])
            # else:
            #     return 1#None
            # # END EXPERIMENTAL
            partitions[ids[i]] += [e for e in p1[i] if cutedge[e]]
            inner_trees[ids[i]] |= {e for e in p1[i] if not cutedge[e]}
            #partitions[i] += [e for e in p1[ds._id[i]] if cutedge[e]]
        for i in range(len(p1), len(p1) + len(p2)):
            # # EXPERIMENTAL
            # print(inner_trees[ids[i]])
            # print('\t', {e for e in p2[i-len(p1)] if not cutedge[e]})
            # if len(inner_trees[ids[i]] & {e for e in p2[i-len(p1)] if not cutedge[e]}) <= 1:
            #     inner_trees[ids[i]] = inner_trees[ids[i]] | set(p2[i-len(p1)])
            # else:
            #     return 2#None
            # # END EXPERIMENTAL
            partitions[ids[i]] += [e for e in p2[i-len(p1)] if cutedge[e]]
            inner_trees[ids[i]] |= {e for e in p2[i-len(p1)] if not cutedge[e]}
            #partitions[i] += [e for e in p2[ds._id[i]-len(p1)] if cutedge[e]]
        
        # Check for internal components
        for i in range(len(partitions)):
            pset = partitions[i]
            intern = inner_trees[i]
            verts = set([edge2psets[e][0] for e in intern]+[edge2psets[e][1]+len(p1) for e in intern])

            #if pset == []:
            if (intern != set() and len(verts)-len(intern) != 1) or pset == []:
                return None
            
        # else:
        #     for pset in partitions:
        #         if pset == [] or (all(bmap[e] == -1 for e in pset)  and root & pset == {}):
        #             return None 

        # print(p1, p2)
        # print('\t', partitions)
        part = Partition(partitions)
        if self.root != -1:
            part.root = ids[self.root]
        elif partition2.root != -1:
            part.root = ids[partition2.root + len(p1)]
        return part# + outer_sets

    
    def join_old(self, partition2, djoin=None, params=None):
        
        if params is not None and params.prototype:
            return self.join2(partition2, djoin)

        if djoin is None:
            inters_elems = self.fullset & partition2.fullset
            # out_elems = (self.fullset | partition2.fullset) - inters_elems
        else:
            inters_elems = set(djoin["non-border"])
            # out_elems = set(djoin["border"])

        conn_comp = []
        remaining_parts = self.partition.copy() + partition2.partition.copy()
        remaining_parts = [set(pset) for pset in remaining_parts]   # MOD
        num_inner_trees = 0
        while len(remaining_parts) != 0:
            accum = remaining_parts[0]

            del remaining_parts[0]

            updated = True
            while updated:
                
                to_rem = []
                updated = False
                for i in range(len(remaining_parts)):
                    # if len(accum & remaining_parts[i]) != 0:  # TODO: FIX: if == 1: do this. elif != 0: return none
                    #     accum = accum | remaining_parts[i]
                    #     to_rem.append(i)
                    #     updated = True
                    
                    if len(accum & remaining_parts[i]) == 1:
                        accum = accum | remaining_parts[i]
                        to_rem.append(i)
                        updated = True
                    elif len(accum & remaining_parts[i]) != 0:
                        #print("R1")
                        return None  # Self loop
                    
            
                for i in reversed(to_rem): 
                    del remaining_parts[i]

            if accum - inters_elems != set():
                conn_comp.append(accum - inters_elems)
            else:
                # how do we check num components?
                # return None
                num_inner_trees += 1
                if num_inner_trees > self.max_inner_trees:
                    return None
        
        res = Partition(conn_comp, self.cost + partition2.cost)
        res.num_inner_trees = num_inner_trees
        return res


    def join2(self, partition2, djoin=None):
        """ Join operation using networkx graphs """

        num_inner_trees = self.num_inner_trees + partition2.num_inner_trees
        if num_inner_trees > self.max_inner_trees:
            return None

        if djoin is None:
            inters_elems = self.fullset & partition2.fullset
            out_elems = (self.fullset | partition2.fullset) - inters_elems
        else:
            inters_elems = set(djoin["non-border"])
            out_elems = set(djoin["border"])
        
        # Plan out graph
        p1 = [list(pset) for pset in self.partition]
        p2 = [list(pset) for pset in partition2.partition]

        # Edgemap: edge name -> (pset1, pset2)
        #      or: edge name -> (pset1, new node)
        edgemap = dict()
        
        for i in range(len(p1)):
            for j in range(len(p1[i])):
                edgemap[p1[i][j]] = (i,)
        
        for i in range(len(p2)):
            for j in range(len(p2[i])):
                if p2[i][j] in edgemap:
                    edgemap[p2[i][j]] = (edgemap[p2[i][j]][0], i + len(p1))
                else:
                    edgemap[p2[i][j]] = (i+len(p1),)
        
        # Create graph
        G = nx.MultiGraph()
        newnode = len(p1) + len(p2)
        for key in edgemap:
            if len(edgemap[key]) == 1:
                k = G.add_edge(edgemap[key][0], newnode, name=key, inner = False)
                edgemap[key] = (edgemap[key][0], newnode, k)
                G.nodes[edgemap[key][0]]['inner'] = False
                G.nodes[newnode]['inner'] = False
                newnode += 1
            else:
                k = G.add_edge(edgemap[key][0], edgemap[key][1], name=key, inner=True)
                edgemap[key] = (edgemap[key][0], edgemap[key][1], k)
                G.nodes[edgemap[key][0]]['inner'] = True
                G.nodes[edgemap[key][1]]['inner'] = True

        #print(edgemap)
        # Check graph
        # (1) Tree check
        comps = list(nx.connected_components(G))
        numcomps = len(list(comps))
        
        if G.number_of_nodes() - G.number_of_edges() - numcomps != 0:
            #print(p1, "+", p2, "=", "None")
            return None
        else: 
            conn_comp = []
            # print(list(G.nodes()))
            # print(list(G.edges()))
            # print([G.edges[e]['name'] for e in G.edges(keys=True)])
            # print(list(comps))
            
            for comp in comps:
                conn_comp.append([])
                Gsub = G.subgraph(comp)
                if all([Gsub.edges[e]['inner'] for e in Gsub.edges(keys=True)]):
                    #print("Inner tree")
                    num_inner_trees += 1
                    if num_inner_trees > self.max_inner_trees:
                        #print(p1, "+", p2, "=", "None")
                        return None
                
                for e in Gsub.edges(keys=True):
                    if not G.edges[e]['inner']:
                        conn_comp[-1].append(G.edges[e]['name'])
                # for n in Gsub.nodes():
                #     if not G.nodes[n]['inner']:
                #         conn_comp[-1].append(G.nodes[n]['name'])

            res = Partition(conn_comp, self.cost + partition2.cost, max_inner_trees=self.max_inner_trees)
            res.num_inner_trees = num_inner_trees
            #print(p1, "+", p2, "=", conn_comp)
            if conn_comp == []:
                raise
            return res

    def join3(self, partition2, djoin=None):
        """ Join operation using scipy """
        num_inner_trees = self.num_inner_trees + partition2.num_inner_trees
        if num_inner_trees > self.max_inner_trees:
            return None

        if djoin is None:
            inters_elems = self.fullset & partition2.fullset
            out_elems = (self.fullset | partition2.fullset) - inters_elems
        else:
            inters_elems = set(djoin["non-border"])
            out_elems = set(djoin["border"])
        
        # Plan out adjacency list
        p1 = [list(pset) for pset in self.partition]
        p2 = [list(pset) for pset in partition2.partition]

        num_nodes = len(p1) + len(p2) + len(list(out_elems))

        all_elems = list(self.fullset)
        #edgeindex = {i : all_elems[i] for i in range(len(all_elems))}

        edgemap = dict()
        
        for i in range(len(p1)):
            for j in range(len(p1[i])):
                edgemap[p1[i][j]] = (i,)
        
        for i in range(len(p2)):
            for j in range(len(p2[i])):
                if p2[i][j] in edgemap:
                    edgemap[p2[i][j]] = (edgemap[p2[i][j]][0], i + len(p1))
                else:
                    edgemap[p2[i][j]] = (i+len(p1),)
        
        # Create adjacency list
        A = np.zeros((num_nodes, num_nodes))
        new_node = len(p1) + len(p2)
        for key in edgemap:
            if len(edgemap[key]) == 1:
                A[edgemap[key][0], new_node] = 1
                new_node += 1
            else:
                A[edgemap[key][0], edgemap[key][1]] = 2

        numcomps, labels = csg.connected_components(A)

        # Check cycles
        if num_nodes - len(all_elems) - numcomps != 0:
            return None

        # Check inner components
        unique = np.unique(labels)
        #inner_components = 0
        for u in unique:
            subG = A[log, log]
            
            log = labels == u
            if np.all(np.logical_not(subG == 2)):
                #inner_components += 1
                return None
            # Add pset


    def connected_uf_old(self, partition2):
        p1 = self.partition
        p2 = partition2.partition
        #cutedge = djoin['border-dict']
        interedge = self.fullset#djoin['non-border']
        # partitions = p1 + p2
        # cutedges = [e for e in cutedge if cutedge[e]]
        # edge2psets = {p1[i][j]:[i, 0] for i in range(len(p1)) for j in range(len(p1[i]))   }
        edge2psets = {e:[i, 0] for i in range(len(p1)) for e in p1[i]   }

        # for i in range(len(p2)):
        #     for j in range(len(p2[i])):
        #         if p2[i][j] in edge2psets:
        #             edge2psets[p2[i][j]][1] = i
        for i in range(len(p2)):
            for e in p2[i]:
                if e in edge2psets:
                    edge2psets[e][1] = i

        ds = DisjointSet()
        [ds.union(edge2psets[e][0], edge2psets[e][1]+len(p1)) for e in interedge]
        
        if len(list(ds.itersets())) == 1:
            return True
        else:
            return False


    def connected_uf(self, partition2):
        p1 = self.partition
        p2 = partition2.partition
        # if len(p1) + len(p2) > 30:
        #     return self.connected_old(partition2)
        interedge = self.fullset
        #partitions = p1 + p2
        #cutedges = [e for e in cutedge if cutedge[e]]
        #edge2psets = {p1[i][j]:[i, 0] for i in range(len(p1)) for j in range(len(p1[i]))   }
        edge2psets = {e:[i, 0] for i in range(len(p1)) for e in p1[i]   }
        for i in range(len(p2)):
            for e in p2[i]:
                if e in edge2psets:
                    edge2psets[e][1] = i
                    
        #outer_sets = [pset for pset in partitions if all(cutedge[e] for e in pset)]

        ds = UF(len(p1) + len(p2))
        [ds.union(edge2psets[e][0], edge2psets[e][1]+len(p1)) for e in interedge]
        
        # partitions = [[] for i in range(max(ds._id)+1)]
        
        ids = ds.itersets()
        
        # #print(ds._id, ':', ds._rank)
        # for i in range(len(p1)):
        #     partitions[ids[i]] += [e for e in p1[i] if cutedge[e]]
        #     #partitions[i] += [e for e in p1[ds._id[i]] if cutedge[e]]
        # for i in range(len(p1), len(p1) + len(p2)):
        #     partitions[ids[i]] += [e for e in p2[i-len(p1)] if cutedge[e]]
        #     #partitions[i] += [e for e in p2[ds._id[i]-len(p1)] if cutedge[e]]

        # DEBUGCHECK

        # END DEBUGCHECK
        inner_trees = [set() for i in range(max(ids)+1)]
        
        #print(ds._id, ':', ds._rank)
        for i in range(len(p1)):
            # # EXPERIMENTAL
            # print(inner_trees[ids[i]])
            # print('\t', {e for e in p1[i] if not cutedge[e]})
            # if len(inner_trees[ids[i]] & {e for e in p1[i] if not cutedge[e]}) <= 1:
            #     inner_trees[ids[i]] = inner_trees[ids[i]] | set(p1[i])
            # else:
            #     return 1#None
            # # END EXPERIMENTAL
            #partitions[ids[i]] += [e for e in p1[i] if cutedge[e]]
            inner_trees[ids[i]] |= {e for e in p1[i]}
            #partitions[i] += [e for e in p1[ds._id[i]] if cutedge[e]]
        for i in range(len(p1), len(p1) + len(p2)):
            # # EXPERIMENTAL
            # print(inner_trees[ids[i]])
            # print('\t', {e for e in p2[i-len(p1)] if not cutedge[e]})
            # if len(inner_trees[ids[i]] & {e for e in p2[i-len(p1)] if not cutedge[e]}) <= 1:
            #     inner_trees[ids[i]] = inner_trees[ids[i]] | set(p2[i-len(p1)])
            # else:
            #     return 2#None
            # # END EXPERIMENTAL
            #partitions[ids[i]] += [e for e in p2[i-len(p1)] if cutedge[e]]
            inner_trees[ids[i]] |= {e for e in p2[i-len(p1)]}
            #partitions[i] += [e for e in p2[ds._id[i]-len(p1)] if cutedge[e]]
        
        # Check for internal components
        for i in range(len(inner_trees)):
            intern = inner_trees[i]
            verts = set([edge2psets[e][0] for e in intern]+[edge2psets[e][1]+len(p1) for e in intern])

            #if pset == []:
            if (intern != set() and len(verts)-len(intern) != 1):
                # assert  False == self.connected_old(partition2), 'Connected old not false' + str(self.partition) + ' ' + str(partition2.partition)
                return False

        if max(ids) == min(ids):
            res = True
        else:
            res = False
        
        #print(ds._id)
        # assert res == self.connected_old(partition2), 'Connected not equal' + str(self.partition) + ' ' + str(partition2.partition)
        return res


        
    def connected_old_uf(self,part2):
        # Check that there is exactly 1 component
        # No cycles are possible, right?
        p1, p2 = [list(pset) for pset in self.partition], [list(pset) for pset in part2.partition]
        #cutedge = set(djoin['border'])
        
        parent_diction = {e:e for pset in p2 for e in pset}
        for pset in p1:
            for e in pset:
                parent_diction[e]=pset[-1]
        ds = DisjointSet(parent_diction)
        [[ds.union(pset[i], pset[i+1]) for i in range(len(pset)-1)] for pset in p2]
        
        # partition=[]
        # 
        # Root/connected check
        # if root_edges is None:
        #     for s in ds.itersets():
        #         intersection=s & cutedge
        #         out_edges = [djoin['bmap'][e] == -1 for e in intersection]
        #         if all(out_edges): # len(intersection)==0
        #             return None
        #         else:
        #             partition.append(intersection)
        # else:
        #     root_found = False # To avoid checking for multiple roots
        #     for s in ds.itersets():
        #         intersection=s & cutedge
        #         out_edges = [djoin['bmap'][e] == -1 for e in intersection]
        #         if not root_found:
        #             root_inter = root_edges & s
        #             if len(root_inter) != 0:
        #                 root_found =True
        #         else:
        #             root_inter = {}
        #         if  all(out_edges) and len(root_inter) == 0: # len(intersection)==0
        #             return None
        #         else:
        #             partition.append(intersection)
                    
        if len(list(ds.itersets())) == 1:
            return True
        else:
            return False

    def connected_old(self, partition2):
        inters_elems = self.fullset & partition2.fullset
        out_elems = (self.fullset | partition2.fullset) - inters_elems

        conn_comp = []
        remaining_parts = self.partition.copy() + partition2.partition.copy()
        remaining_parts = [set(pset) for pset in remaining_parts]   # MOD
        while len(remaining_parts) != 0:
            accum = remaining_parts[0]

            del remaining_parts[0]

            updated = True
            while updated:
                
                to_rem = []
                updated = False
                for i in range(len(remaining_parts)):
                    # if len(accum & remaining_parts[i]) != 0:  # TODO: FIX: if == 1: do this. elif != 0: return none
                    #     accum = accum | remaining_parts[i]
                    #     to_rem.append(i)
                    #     updated = True
                    
                    if len(accum & remaining_parts[i]) == 1:
                        accum = accum | remaining_parts[i]
                        to_rem.append(i)
                        updated = True
                    elif len(accum & remaining_parts[i]) != 0:
                        return False  # Self loop
                    
            
                for i in reversed(to_rem): 
                    del remaining_parts[i]

            if accum - inters_elems != set():
                conn_comp.append(accum - inters_elems)
            elif accum == inters_elems:
                if len(conn_comp) <= 1:
                    return True
                else:
                    return False
            else:
                #assert 2==3, 'Self loop on final step'
                return False
        
        if len(conn_comp) == 1:
            return True
        else:
            return False
        #return Partition(conn_comp, self.cost + partition2.cost)



    @staticmethod
    def all_partitions(in_list: List[int]):
        """ Given a single list of integers, return all possible partitions of that list. """
        in_list_sorted = sorted(in_list)

        partitions: List[PartitionT] = [[[in_list_sorted[0]]]]

        for elem in in_list_sorted[1:]:
            new_partitions: List[PartitionT] = []
            for partition in partitions:
                # New partitions with the new element at each index
                for i in range(len(partition)):
                    new_partitions.append(partition[:i] + [partition[i] + [elem]] + partition[(i+1):])
                # The old partitions will now have a set with the new element only
                partition.append([elem])

            partitions.extend(new_partitions)
        
        # Convert to class:
        class_partitions: List[Partition] = []
        for partition in partitions:
            class_partitions.append(Partition(partition))

        return class_partitions


# class DirectedPartition:
#     """ 
#     Class for handling directed partition. Directed partitions have
#     partition sets with format (parent, Set), the parent is included in Set.
#     """
#     def __init__(self, partition: DPartitionT, cost: float=0):
#         self.partition = [0]*len(partition)
#         self.parents = [0]*len(partition)

#         # Change to set
#         for i in range(len(partition)):
#             (parent, Set) = partition[i]
#             self.partition[i] = set(Set)
#             self.parents[i] = set(parent)
        
#         self.fullset = reduce(set.union, partition, set())

#         self.encoding = None
#         self.encoding = self.get_encoding()
#         self.cost = cost
    

#     def get_encoding(self) -> int:
#         """ 
#         Returns a unique encoding for a partition 
#         """
#         if self.encoding is None:  # Avoid calculating twice
#             items = self.fullset
#             partition = self.partition

#             sorted_items = sorted(items)
#             self.item_map = dict()
#             for i in range(len(sorted_items)): 
#                 self.item_map[sorted_items[i]] = i

#             encoding_list = [0]*len(items)
#             self.encoding_map = dict()  # useful data structure for user. Gets item, returns partition ind
#             for i in range(len(partition)):
#                 pset = partition[i]
#                 for item in pset:
#                     item_ind = self.item_map[item]
#                     encoding_list[item_ind] = i
#                     self.encoding_map[item] = i  # For future use
            
#             base = len(items)
#             encoding = 0
#             part_ind = 0
#             pmap = dict()

#             for i in reversed(range(len(encoding_list))):
#                 if encoding_list[i] not in pmap:
#                     pmap[encoding_list[i]] = part_ind
#                     part_ind += 1
#                 encoding += base**i * pmap[encoding_list[i]]
            
#             # Now deal with parents, since this is directed
#             for i in range(len(self.parents)):
#                 p_encoding += base**i * pmap[self.parents[i]]
            
#             encoding = (p_encoding, encoding)

#         else:
#             encoding = self.encoding

                    

#         return encoding

    
#     def join(self, partition2):
#         inters_elems = self.fullset & partition2.fullset
#         out_elems = (self.fullset | partition2.fullset) - inters_elems

#         conn_comp = []
#         remaining_parts = self.partition.copy() + partition2.partition.copy()
#         while len(remaining_parts) != 0:
#             accum = remaining_parts[0]

#             del remaining_parts[0]

#             updated = True
#             while updated:
                
#                 to_rem = []
#                 updated = False
#                 for i in range(len(remaining_parts)):
#                     # if len(accum & remaining_parts[i]) != 0:  # TODO: FIX: if == 1: do this. elif != 0: return none
#                     #     accum = accum | remaining_parts[i]
#                     #     to_rem.append(i)
#                     #     updated = True
                    
#                     if len(accum & remaining_parts[i]) == 1:
#                         accum = accum | remaining_parts[i]
#                         to_rem.append(i)
#                         updated = True
#                     elif len(accum & remaining_parts[i]) != 0:
#                         #print("R1")
#                         return None  # Self loop
                    
            
#                 for i in reversed(to_rem): 
#                     del remaining_parts[i]

#             if accum - inters_elems != set():
#                 conn_comp.append(accum - inters_elems)
#             else:
#                 #print("R2")
#                 return None
        
#         return Partition(conn_comp, self.cost + partition2.cost)

#     def connected(self, partition2):
#         inters_elems = self.fullset & partition2.fullset
#         out_elems = (self.fullset | partition2.fullset) - inters_elems

#         conn_comp = []
#         remaining_parts = self.partition.copy() + partition2.partition.copy()
#         while len(remaining_parts) != 0:
#             accum = remaining_parts[0]

#             del remaining_parts[0]

#             updated = True
#             while updated:
                
#                 to_rem = []
#                 updated = False
#                 for i in range(len(remaining_parts)):
#                     # if len(accum & remaining_parts[i]) != 0:  # TODO: FIX: if == 1: do this. elif != 0: return none
#                     #     accum = accum | remaining_parts[i]
#                     #     to_rem.append(i)
#                     #     updated = True
                    
#                     if len(accum & remaining_parts[i]) == 1:
#                         accum = accum | remaining_parts[i]
#                         to_rem.append(i)
#                         updated = True
#                     elif len(accum & remaining_parts[i]) != 0:
#                         return False  # Self loop
                    
            
#                 for i in reversed(to_rem): 
#                     del remaining_parts[i]

#             if accum - inters_elems != set():
#                 conn_comp.append(accum - inters_elems)
#             elif accum == inters_elems:
#                 if len(conn_comp) <= 1:
#                     return True
#                 else:
#                     return False
#             else:
#                 return False
        
#         if len(conn_comp) == 1:
#             return True
#         else:
#             return False
#         #return Partition(conn_comp, self.cost + partition2.cost)



#     @staticmethod
#     def all_partitions(in_list: List[int]):
#         """ Given a single list of integers, return all possible partitions of that list. """
#         # Do all combinations of the sets 
#         in_list_sorted = sorted(in_list)

#         partitions: List[PartitionT] = [[[in_list_sorted[0]]]]

#         for elem in in_list_sorted[1:]:
#             new_partitions: List[PartitionT] = []
#             for partition in partitions:
#                 # New partitions with the new element at each index
#                 for i in range(len(partition)):
#                     new_partitions.append(partition[:i] + [partition[i] + [elem]] + partition[(i+1):])
#                 # The old partitions will now have a set with the new element only
#                 partition.append([elem])

#             partitions.extend(new_partitions)
        
#         # Convert to class:
#         class_partitions: List[Partition] = []
#         for partition in partitions:
#             class_partitions.append(Partition(partition))
        
#         # Now do all combinations of the parents
#         all_parents

#         return class_partitions






"""This module implements an union find or disjoint set data structure.

An union find data structure can keep track of a set of elements into a number
of disjoint (nonoverlapping) subsets. That is why it is also known as the
disjoint set data structure. Mainly two useful operations on such a data
structure can be performed. A *find* operation determines which subset a
particular element is in. This can be used for determining if two
elements are in the same subset. An *union* Join two subsets into a
single subset.

The complexity of these two operations depend on the particular implementation.
It is possible to achieve constant time (O(1)) for any one of those operations
while the operation is penalized. A balance between the complexities of these
two operations is desirable and achievable following two enhancements:

1.  Using union by rank -- always attach the smaller tree to the root of the
    larger tree.
2.  Using path compression -- flattening the structure of the tree whenever
    find is used on it.

complexity:
    * find -- :math:`O(\\alpha(N))` where :math:`\\alpha(n)` is
      `inverse ackerman function
      <http://en.wikipedia.org/wiki/Ackermann_function#Inverse>`_.
    * union -- :math:`O(\\alpha(N))` where :math:`\\alpha(n)` is
      `inverse ackerman function
      <http://en.wikipedia.org/wiki/Ackermann_function#Inverse>`_.

"""


class UF:
    """An implementation of union find data structure.
    It uses weighted quick union by rank with path compression.
    """

    def __init__(self, N):
        """Initialize an empty union find object with N items.

        Args:
            N: Number of items in the union find object.
        """

        self._id = list(range(N))
        self._count = N
        self._rank = [0] * N

    def find(self, p):
        """Find the set identifier for the item p."""

        id = self._id
        while p != id[p]:
            p = id[p] = id[id[p]]   # Path compression using halving.
        return p

    def count(self):
        """Return the number of items."""

        return self._count

    def connected(self, p, q):
        """Check if the items p and q are on the same set or not."""

        return self.find(p) == self.find(q)

    def union(self, p, q):
        """Combine sets containing p and q into a single set."""

        id = self._id
        rank = self._rank

        i = self.find(p)
        j = self.find(q)
        if i == j:
            return

        self._count -= 1
        if rank[i] < rank[j]:
            id[i] = j
        elif rank[i] > rank[j]:
            id[j] = i
        else:
            id[j] = i
            rank[i] += 1

    def __str__(self):
        """String representation of the union find object."""
        return " ".join([str(x) for x in self._id])

    def __repr__(self):
        """Representation of the union find object."""
        return "UF(" + str(self) + ")"
    
    def itersets(self):
        ids = [self.find(self._id[i]) for i in range(len(self._id))]
        mapping = {}
        nex = 0
        for id in ids:
            if id not in mapping:
                mapping[id] = nex
                nex = nex + 1
        return [mapping[id] for id in ids]



