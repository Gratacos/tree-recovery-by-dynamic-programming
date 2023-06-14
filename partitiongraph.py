import networkx as nx
import partitions as prt
from beamlist import BeamList
from queue import PriorityQueue
import os
import json as js
import random
import mpmath as mp
from typing import List
from scipy.spatial import Delaunay
from itertools import product
from random import random, seed
import numpy as np
import dpartitions as dp
import math


DIRECTED_TIMING = 0


def get_time():
    return getattr(os.times(), 'user')


mp.mp.dps = 30

def number_directed_partitions(n):
    nCk = lambda n, k: math.factorial(n) / math.factorial(k) / math.factorial(n-k)
    idempotent = lambda n, k: nCk(n, k) * k**(n-k)
    return sum([idempotent(n, k) for k in range(1,int(n)+1)])

stats = {
    'all': [],
    'chosen' : [],
    'per join cost' : []
}

class Params:
    def __init__(self,
             beam=True,
             variable_beam=False,
             init_width=float("inf"), 
             beam_heap=True,
             qcost="outdegree",
             merged_cost=999999,
             beam_type='fixed',
             min_width=10,
             verbose=True,
             prototype=False,
             histogram=False,
             minnorm=True,
             norm='none'):
             
             self.beam = beam
             self.variable_beam = variable_beam
             self.init_width = init_width
             self.qcost = qcost
             self.merged_cost = merged_cost
             self.beam_heap = beam_heap
             self.min_width = min_width
             self.prototype = prototype
             self.verbose = verbose
             self.histogram = histogram
             self.minnorm = minnorm 
             self.beam_type = beam_type
             self.norm = norm

             if self.beam_type == 1:
                 self.beam_type = 'fixed'
             elif self.beam_type == 2:
                 self.beam_type = 'percent'
             elif self.beam_type == 3:
                 self.beam_type = 'probability'

params = Params()

class PGraph(nx.MultiGraph):
    def __init__(self, *args, **kwargs):
        super(PGraph, self).__init__(*args, **kwargs)
        
        self.edge_map = {self.edges[edge]['name']:edge for edge in self.edges(keys=True)}
        
        # self.next_edge_name = 0
        # self.edge_map = dict()
        # for edge in self.edges(keys=True):
        #     self.edges[edge]["name"] = self.next_edge_name
        #     self.edge_map[self.next_edge_name] = edge
        #     self.next_edge_name += 1

        # if len(self.nodes()) == 0: self.next_vertex_name = 0    
        # else: self.next_vertex_name = max([int(n) for n in self.nodes()])+1
        
        #for node in self.nodes():
        #    neighboring_edges = {self.edges[e]["name"] for e in self.edges(node)}
        #    #self.nodes[node]["partition"] = prt.Partition([neighboring_edges])
    
    # def add_edge(self, *args, **kwargs): #  GUSMOD
    #     k = super(PGraph, self).add_edge(*args, **kwargs)
    #     edge = (args[0], args[1], k)

    #     self.next_vertex_name = max(max(self.next_vertex_name, int(args[0])+1), int(args[1])+1)

        
    #     if "name" not in self.edges[edge]:
    #         self.edges[edge]["name"] = self.next_edge_name
    #         self.next_edge_name = self.next_edge_name + 1
        
    #     self.next_edge_name = max(self.next_edge_name, self.edges[edge]["name"] + 1)
    #     self.edge_map[ self.edges[edge]["name"] ] = edge
        
        
    #     return k

    # def add_node(self, *args, **kwargs): # GUSMOD
    #     super(PGraph, self).add_node(*args, **kwargs)
    #     # if "partition" not in self.nodes[args[0]]:
    #     #     neighboring_edges = {self.edges[e]["name"] for e in self.edges(args[0], keys=True)}
    #     #     #self.nodes[args[0]]["partition"] = prt.Partition([neighboring_edges])

    #     # self.next_vertex_name = max(self.next_vertex_name, int(args[0]) + 1)

    
    def partition_graph(self):
        """ Returns the partitioned graph """
        #return_graph = PGraph(self)
        return_graph = PGraph()
        for e in self.edges(keys=True):
            return_graph.add_edge(*e)
            return_graph.edges[e]["name"] = self.edges[e]["name"]
            
        return_graph.edge_map = self.edge_map.copy()


        for node in self.nodes():
            if len(self.nodes[node]["partition"].partition) != 1:
                for pset in self.nodes[node]["partition"].partition:#.partition[1:]:
                    vp = return_graph.copy_vertex(node)
                    for edge_name in pset:
                        (e1, e2, k) = self.edge_map[edge_name]
                        (e1, e2, k) = return_graph.edge_map[edge_name]

                        if e1 == node:
                            e1, e2 = e2, e1

                        return_graph.move_edge((e1, e2, k), (e1, vp))
                
                return_graph.remove_node(node)
        
        return return_graph
    
    def copy_vertex(self, v):
        """utility function to copy a vertex and all of its attributes
           without copying its adjacent edges. 
           also, return the new vertex's name"""

        node = self.next_vertex_name
        self.add_node(node)
        
        for (key, val) in self.nodes[v].items():
            self.nodes[node][key] = val
        
        self.nodes[node]["cluster"] = v
        self.nodes[v]["cluster"] = v

        return node
    
    def move_edge(self, e1, e2):
        k= self.add_edge(*e2, name=self.edges[e1]["name"])
        e2 = (e2[0], e2[1], k)  # In case user didn't specify k for second edge

        self.edge_map[self.edges[e1]["name"]] = e2

        for (key, val) in self.edges[e1].items():
            self.edges[e2][key] = val
        
        self.remove_edge(*e1)
        return k
    
    def join_vertices(self, v1, v2):
        """ Metadata of v1 will be kept """

        if True:# (v1, v2) in self.edges():
            # Deal with partition and cost
            ### MOD
            # if "partition" in self.nodes[v1] and "partition" in self.nodes[v2]:
            #     p1 = self.nodes[v1]["partition"]
            #     p2 = self.nodes[v2]["partition"]
            # else:
            #     p1 = None
            #     p2 = None
            # if p1 is not None and p2 is not None: 
            #     p3 = p1.join(p2, params=params)
            # else:
            #     p3 = None
            # #if p3 is None: return False

            # self.nodes[v1]["partition"] = p3

            # if p3 is not None:
            #     self.nodes[v1]["partition"].cost = p1.cost + p2.cost
            ### MOD

            # Deal with vertices and edges:
            #print("---")
            neigh_edges = list(self.edges(v2, keys=True)).copy()
            for (u1, u2, k) in neigh_edges:
                if u1 == v2:
                    u1, u2 = u2, u1  # Make it so that (u1, u2) = (u1, v2)
                if u1 != v1:
                    #print(self.edges[(u1, u2, k)]["name"], " | ", (u1, u2, k), "->", (v1, u1))
                    self.move_edge((u1, u2, k), (v1, u1))

            self.remove_node(v2)

            return True
        
        else:
            return False


    def save(self, fname: str, genNames=True):
        """ Save (1) The graph's adjacency list and (2) The graph's *chosen* partitions """
        
        # if has an extension, replace it, we will only support one extension
        fname = fname.split(".")
        fname = ".".join(fname[:-1]) + ".pgraph"

        # If file exists, delete first
        if os.path.exists(fname):
            os.remove(fname)

        # Save the file
        f = open("demofile2.txt", "a")
        f.write("Now the file has more content!")
        f.close()

        # Creates a new file 
        with open(fname, 'w') as graphFile: 

            print("writing to file", fname)
            graphFile.write("pgraph\n")

            # Extract some data for header
            numEdges = str(self.number_of_edges())
            numNodes = str(self.number_of_nodes())

            graphFile.write("edges " + numEdges + "\n" + "nodes " + numNodes + "\n")

            # Write edges:
            if genNames:
                next_name = 1
                for edge in self.edges(keys=True):
                    self.edges[edge]["name"] = str(next_name)
                    next_name += 1

            for (u, v, k) in self.edges(keys=True):
                name = self.edges[(u, v, k)]["name"] 
                graphFile.write(str(u) + " " + str(v) + " " + name + "\n")
            
            # Write vertices:
            for node in self.nodes():
                if 'x' not in self.nodes[node] or 'y' not in self.nodes[node] or 'z' not in self.nodes[node]:
                    x, y, z = str(-1), str(-1), str(-1)
                else: 
                    x, y, z = str(self.nodes[node]['x']), str(self.nodes[node]['y']), str(self.nodes[node]['z'])
                
                partition_string = "_"
                if "partition" in self.nodes[node]:
                    for pset in self.nodes[node]["partition"].partition:
                        for edge_name in pset:
                            partition_string += str(edge_name) + ","
                        partition_string = partition_string[:-1] + "_"  # Currently, [-1] has a trailing comma
                
                graphFile.write(x + " " + y + " " + z + " " + partition_string + "\n")

    def to_json(self, costs = None, parents=None, names=True, fname=None, directed=False):
        acceptable_node_keys = {'x', 'y', 'z', 'partition', 'partition_candidates', "name", "num_cands", "parent_lists", "directed_chosen", "num_theo_cands", "bell_num", "num_undirected", "parent_bool", "num_join_ops"}
        json = {"nodes" : [], "edges" : [], "costs":costs}

        for key in self.graph:
            json[key] = self.graph[key]
       
        # Format nodes
        snodes = self.nodes()
        list(snodes).sort()
        for node in snodes:
            node_attr = self.nodes[node].copy()
            node_attr["name"] = node

            # "listify" the partitions for json format
            if "partition" in node_attr and node_attr["partition"] is not None:
                node_attr["partition"] = [list(pset) for pset in node_attr["partition"].partition]
            if "partition_candidates" in node_attr:
                node_attr["partition_candidates"] = []
                for cand in self.nodes[node]["partition_candidates"]:

                    if costs is not None:
                        if cand.get_encoding() in costs[node] and not directed:
                            cost = costs[node][cand.get_encoding()]
                        else:
                            cost = -1
                        # else:
                        #     print("")
                        #     print(dp.get_encoding(cand), "->",self.nodes[node]["parent_lists"][dp.get_encoding(cand)])
                        #     print(dp.get_encoding(cand, self.nodes[node]["parent_lists"][dp.get_encoding(cand)]))
                        #     print("")
                        #     #cost = costs[node][dp.get_encoding(cand, self.nodes[n[dp.get_encoding(cand)])]
                        #     cost = costs[node][dp.get_encoding(cand)]
                    else:
                        cost = -1
                    
                    if parents is not None and node in parents and cand.get_encoding() in parents[node]:
                        pars = parents[node][cand.get_encoding()][1:]

                        for i in range(len(pars)):

                            new_pars = dict()
                            for pnode in pars[i]:
                                new_pars[pnode] = [list(pset) for pset in pars[i][pnode].partition]
                            
                            pars[i] = new_pars
                    else:
                        pars = -1

                    node_attr["partition_candidates"].append({"partition":[list(pset) for pset in cand.partition], "cost":cost, "parents":pars})
                
            # Only include the acceptable keys
            node_attr = {key : node_attr[key] for key in node_attr if key in acceptable_node_keys}
            json["nodes"].append(node_attr)
        
        # Format edges
        names = 1
        for (u, v, k) in self.edges(keys=True):
            if names:
                edge = {"link" : [int(u), int(v)], "name" : int(self.edges[(u, v, k)]["name"])}
            else:
                edge = {"link" : [int(u), int(v)], "name" : names}
                names += 1
            json["edges"].append(edge)

        if fname is not None:
            with open(fname, "w") as fp:
                #print(json)
                js.dump(json, fp, indent=4)

        #print("JSON FILE\n", json)
        return js.dumps(json, indent=4)
    
    def from_json(self, json):
        for i in range(len(json["nodes"])):
            node = json["nodes"][i]["name"]
            self.add_node( node )
            for attribute in json["nodes"][i]:
                self.nodes[node][attribute] = json["nodes"][i][attribute]
        
        for edge in json["edges"]:
            [u, v] = edge["link"]
            k = self.add_edge(u, v)
            for attribute in edge:
                self.edges[(u, v, k)][attribute] = edge[attribute]
        
        for attribute in [attr for attr in json if attr not in ["edges", "nodes"]]:
            self.graph[attribute] = json[attribute]


    def load(self, fname: str):
        """ load a pgraph file """
        with open(fname, 'r') as graphFile:
            line = graphFile.readline()  #  "plyfile"
            line = graphFile.readline() 
            line = line.split(" ")
            if line[0] == "nodes":
                num_nodes = int(line[1])
                line = graphFile.readline()
                line = line.split(" ")
                num_edges = int(line[1])
            else:
                num_edges = int(line[1])
                line = graphFile.readline()
                line = line.split(" ")
                num_nodes = int(line[1])
            
            for i in range(num_edges):
                line = graphFile.readline()
                line = line.split(" ")
                
                u = int(line[0])
                v = int(line[1])
                name = line[2]

                k = self.add_edge(u, v)
                self.edges[(u, v, k)]["name"] = name

            for i in range(num_nodes):
                line = graphFile.readline()
                line = line.split(" ")

                x = float(line[0])
                y = float(line[1])
                z = float(line[2])
                partition_raw = line[3]

                # parse partition
                partition = []
                partition_raw = partition_raw.split("_")
                for pset_raw in partition_raw:
                    pset = {e for e in pset_raw.split(",")}
                    partition.append(pset)
                
                self.add_node(i, x=x, y=y, z=z, partition=prt.Partition(partition))


# def beamlist_add(beam_list, item, cost):
#     global params
#     beam_list.add((cost, item))
#     if params.beam:
#         if beam_list >= params.init_width:
#             beam_list.pop()


# Erosion Algorithm
import time

NUMBER_MINS = 0
#@jit(parallel=True)
def joined_cost(G, partition_candidates_1, partition_candidates_2, new_node, n1, n2, cost, back_track):
    """ Given partition candidates and a cost function, what is the cost of the new joined partition? """
    #global NUMBER_MINS
    global params

    seen = set()

    contributions=dict()
    contributions[-1] = (new_node, n1, n2)  # "metadata"
    new_cost = cost.copy()

    new_cost[new_node] = dict()  # Restart this node

    new_part_candidates = BeamList(params)

    timings = []
    gran_times = []
    
    numnotnones = 0

    all_granulated_timings = []

    # print("Join:")
    # print("candidates 1", n1)
    # for p in partition_candidates_1:
    #     print("    ", p.partition)
    # print("candidates 2", n2)
    # for p in partition_candidates_2:
    #     print("    ", p.partition)

    for partition1 in partition_candidates_1:
        for partition2 in partition_candidates_2:
            t0 = time.time()

            join = partition1.join(partition2, params=params)


            if join is not None:
                join_ind = join.get_encoding()

                if join_ind not in seen:
                    seen.add(join_ind)
                    back_track[-1][new_node][join_ind] = (-1, -1)

                p1_ind = partition1.get_encoding()
                p2_ind = partition2.get_encoding()

                if new_node not in new_cost:
                    new_cost[new_node] = dict()
                
                if join_ind not in new_cost[new_node]:
                    new_cost[new_node][join_ind] = float("inf")

                # 
                if new_cost[new_node][join_ind] > cost[n1][p1_ind] + cost[n2][p2_ind]:
                    new_cost[new_node][join_ind] = cost[n1][p1_ind] + cost[n2][p2_ind]
                    back_track[-1][new_node][join_ind] = {n1:partition1, n2:partition2}
                    new_part_candidates.add(new_cost[new_node][join_ind], join, join_ind)

                timings.append(time.time() - t0)

    G.nodes[new_node]["partition_candidates"] = new_part_candidates.get_list()
    if params.histogram:
        global stats
        #print('----', new_part_candidates.per_join_cost)
        stats['per join cost'].append(new_part_candidates.per_join_cost)

    # print("Result")
    # for p in G.nodes[new_node]["partition_candidates"]:
    #     print("    ", p.partition)
    # print(".....................")

    # fift_in_e = False
    # print("")
    # print("All p candidates:")
    # for node in G.nodes():
    #     print("    Node", node)
    #     print("    with edges:", [G.edges[e]["name"] for e in G.edges(node, keys=True)])
        
    #     for p in G.nodes[node]["partition_candidates"]:
    #         print("        ", p.partition)
    #         if 15 in [G.edges[e]["name"] for e in G.edges(node, keys=True)] and 15 not in p.fullset:
    #             fift_in_e = True
    # print("")

    return (new_cost, contributions, timings, numnotnones, all_granulated_timings)

time_mult = 0
time_enc = 0
time_getenc = 0
import dpartitions as dp
def joined_cost_directed(G, partition_candidates_1, pcands1, partition_candidates_2, pcands2, new_node, n1, n2, cost, back_track):
    """ Given partition candidates and a cost function, what is the cost of the new joined partition? """
    #global NUMBER_MINS
    global params
    global time_mult
    global time_enc
    global time_getenc

    seen = set()

    contributions=dict()
    contributions[-1] = (new_node, n1, n2)  # "metadata"
    new_cost = cost.copy()

    new_cost[new_node] = dict()  # Restart this node

    new_part_candidates = BeamList(params)

    timings = []
    gran_times = []
    
    numnotnones = 0

    all_granulated_timings = []
    num_join_ops = 0

    for partition1 in partition_candidates_1:
        for partition2 in partition_candidates_2:
            # Faster to check undirected first:
            ujoin = dp.ujoin(partition1, partition2, params=params)
            if ujoin is None: continue

            for parents1 in pcands1[dp.get_encoding(partition1)]:
                #print("pcands2:", pcands2)
                for parents2 in pcands2[dp.get_encoding(partition2)]:
                    num_join_ops += 1
                    t0 = time.time()

                    join = dp.djoin(ujoin, partition1, partition2, parents1, parents2)
                    if join == False: join = None
                    time_mult += time.time() - t0

                    t0 = time.time()
                    if join is not None:
                        join_undirected, join_parents = join

                        t01 = time.time()
                        join_ind = dp.get_encoding(join_undirected, join_parents)
                        p1_ind = dp.get_encoding(partition1, parents1)
                        p2_ind = dp.get_encoding(partition2, parents2)
                        time_getenc += time.time() - t01

                        if join_ind not in seen:
                            seen.add(join_ind)
                            back_track[-1][new_node][join_ind] = (-1, -1)

                        if new_node not in new_cost:
                            new_cost[new_node] = dict()
                        
                        if join_ind not in new_cost[new_node]:
                            new_cost[new_node][join_ind] = float("inf")

                        if new_cost[new_node][join_ind] > cost[n1][p1_ind] + cost[n2][p2_ind]:
                            new_cost[new_node][join_ind] = cost[n1][p1_ind] + cost[n2][p2_ind]
                            back_track[-1][new_node][join_ind] = {n1:(partition1, parents1), n2:(partition2, parents2)}
                            new_part_candidates.add(new_cost[new_node][join_ind], join, join_ind)

                        timings.append(time.time() - t0)
                    time_enc += time.time() - t0 
    partcands, parentcands = zip(*new_part_candidates.get_list())

    if params.histogram:
        global stats
        #print('----', new_part_candidates.per_join_cost)
        stats['per join cost'].append(new_part_candidates.per_join_cost)

    partcands = list(partcands)
    parentcands = list(parentcands)
    parentcands_dict = {dp.get_encoding(cand):[] for cand in partcands}
    new_part_cands = []  # Only store unique in here
    seen_cand = set()
    for i in range(len(partcands)):
        sorted_part, sorted_par = dp.sort_partition(partcands[i], parentcands[i])
        part_cand_enc = dp.get_encoding(partcands[i])
        parentcands_dict[part_cand_enc].append(sorted_par)
        #partcands[i] = prt.Partition(sorted_part)  # Replace partition with sorted version
        if part_cand_enc not in seen_cand:
            new_part_cands.append(prt.Partition(sorted_part))
            seen_cand.add(part_cand_enc)


    G.nodes[new_node]["partition_candidates"] = new_part_cands#partcands
    G.nodes[new_node]["parent_lists"] = parentcands_dict
    G.nodes[new_node]["num_cands"] = sum([len(parentcands_dict[key]) for key in parentcands_dict])
    G.nodes[new_node]["num_undirected"] = len(new_part_cands)
    G.nodes[new_node]["num_join_ops"] = num_join_ops

    return (new_cost, contributions, timings, numnotnones, all_granulated_timings)


def QCost(G, item):
    cost = params.qcost
    (u, v) = item

    if isinstance(cost, tuple):
        qcost_ret = []
        for c in cost:
            #print(c)
            params.qcost = c
            qcost_ret.append(QCost(G, item))
        params.qcost = cost
        return tuple(qcost_ret)

    if cost == "maxdegree":
        return max(nx.degree(G, u), nx.degree(G, v))
    elif cost == "prod":
        #return nx.degree(G, u) * nx.degree(G, v) # TODO: Partition candidates?
        if "partition_candidates" in G.nodes[u] and "partition_candidates" in G.nodes[v]:
            return len(G.nodes[u]["partition_candidates"]) * len(G.nodes[v]["partition_candidates"])
        else:
            return nx.degree(G, u) * nx.degree(G, v)
    elif cost == "outdegree":
        uedges = {G.edges[e]["name"] for e in G.edges(u, keys=True)}
        vedges = {G.edges[e]["name"] for e in G.edges(v, keys=True)}
        inters_elems = uedges & vedges
        out_elems = (uedges | vedges) - inters_elems
        return len(out_elems) 
    elif cost == "indegree":
        uedges = {G.edges[e]["name"] for e in G.edges(u, keys=True)}
        vedges = {G.edges[e]["name"] for e in G.edges(v, keys=True)}
        inters_elems = uedges & vedges
        return len(inters_elems)
    elif cost == "boundary":
        # if "partition_candidates" in G.nodes[u] and "partition_candidates" in G.nodes[v]:
        #     prod = len(G.nodes[u]["partition_candidates"]) * len(G.nodes[v]["partition_candidates"])
        # else:
        #     prod = nx.degree(G, u) * nx.degree(G, v)
        uedges = {G.edges[e]["name"] for e in G.edges(u, keys=True)}
        vedges = {G.edges[e]["name"] for e in G.edges(v, keys=True)}
        inters_elems = uedges & vedges
        out_elems = (uedges | vedges) - inters_elems
        prod = len(out_elems)

        if (G.nodes[u]["merged"] == G.nodes[v]["merged"]):
            return prod + params.merged_cost
        else:
            return prod
    elif cost == "diff":
        # Outdegree
        uedges = {G.edges[e]["name"] for e in G.edges(u, keys=True)}
        vedges = {G.edges[e]["name"] for e in G.edges(v, keys=True)}
        inters_elems = uedges & vedges
        out_elems = len((uedges | vedges) - inters_elems)
        bell_out = np.log(float(mp.bell(out_elems)))

        # Prod:
        if "partition_candidates" in G.nodes[u] and "partition_candidates" in G.nodes[v]:
            prod = np.log(len(G.nodes[u]["partition_candidates"])) + np.log(len(G.nodes[v]["partition_candidates"]))
        else:
            prod = nx.degree(G, u) * nx.degree(G, v)

        def bell(n): pass

        return prod - bell_out

    elif cost == "diffinv":
        # Outdegree
        uedges = {G.edges[e]["name"] for e in G.edges(u, keys=True)}
        vedges = {G.edges[e]["name"] for e in G.edges(v, keys=True)}
        inters_elems = uedges & vedges
        out_elems = len((uedges | vedges) - inters_elems)
        bell_out = np.log(float(mp.bell(out_elems)))

        # Prod:
        if "partition_candidates" in G.nodes[u] and "partition_candidates" in G.nodes[v]:
            prod = np.log(len(G.nodes[u]["partition_candidates"])) + np.log(len(G.nodes[v]["partition_candidates"]))
        else:
            prod = nx.degree(G, u) * nx.degree(G, v)

        def bell(n): pass

        return bell_out - prod
        


def QPut(G, pq, node, max_degree):
    # Find neighbor with smallest degree (pq is modified)
    best_deg = float("inf")
    best_deg_neigh = -1
    for neigh in nx.neighbors(G, node):
        if neigh != node and nx.degree(G, neigh) < best_deg:
            best_deg = nx.degree(G, neigh)
            best_deg_neigh = neigh
    # Put neighbor node pair in queue
    qcost = QCost(G, (node, best_deg_neigh))
    #if qcost <= max_degree:
    if min(nx.degree(G, node), nx.degree(G, best_deg_neigh)) <= max_degree:
        pq.put( (qcost, (node, best_deg_neigh)) )

def QGet(G, pq, max_degree):
    # Get next item with accurate degree
    item_obtained = False
    item = None
    while not item_obtained and not pq.empty():
        (deg, (node, neigh)) = pq.get()
        if (node in G.nodes() and neigh in G.nodes()) and QCost(G, (node, neigh)) != deg:
            pq.put((QCost(G, (node, neigh)), (node, neigh)))
        elif node in G.nodes() and neigh in G.nodes():
            item = (QCost(G, (node, neigh)), (node, neigh))
            item_obtained = True
        elif node in G.nodes():
            QPut(G, pq, node, max_degree)
        elif neigh in G.nodes():
            QPut(G, pq, neigh, max_degree)
    
    return item 
            

def erosion(G: PGraph, max_degree: int = 2, max_iter=float("inf"), cost: dict=None, store_graphs=True, is_directed=False, bits=False, nodict=False):

    print("\n\nStarting erosion\n")

    for node in G.nodes():
        G.nodes[node]["merged"] = False

    # Some plot data
    times = []
    cand_prod = []
    degrees = []
    ndegrees = []
    vjointimes = []
    timings = []
    numnotnones = []
    gran = []
    all_cands = []
    graphs = []
    contributions = []
    back_track = []
    sumlogs = []
    max_cands = []  # The maximum candidate accross iterations up to time point

    # Calculate sum of logs:
    sumlogs.append(0)
    for node in G.nodes():
        sumlogs[0] += np.log(len(G.nodes[node]["partition_candidates"]))

    max_cands.append(max( [len(G.nodes[v]["partition_candidates"]) for v in G.nodes()] ))

    if cost is not None: 
        new_cost = cost.copy()
    else: 
        new_cost = dict()

    # Initialize priority queue with cost: max(degree(node), degree(neigh)) and item (node, neigh)
    next_pair = PriorityQueue()
    for (u, v) in G.edges():
        deg = QCost(G, (u, v))
        if min(nx.degree(G, u), nx.degree(G, v)) <= max_degree and u != v:
            next_pair.put( (deg, (u, v)) )

    # Begin main loop
    T = 1
    while not next_pair.empty() and T <= max_iter and G.number_of_nodes() > 2:

        t0 = time.time()

        item = QGet(G, next_pair, max_degree)
        if item is None:
            break
        (deg, (node, neighbor)) = item
        if node not in G.nodes():
            continue

        sumlogs.append(sumlogs[-1])
        sumlogs[-1] -= np.log(len(G.nodes[node]["partition_candidates"])) + np.log(len(G.nodes[neighbor]["partition_candidates"])) 

        ndegs = (nx.degree(G, node), nx.degree(G, neighbor))

        if "partition_candidates" in G.nodes[node]:
            ncand = (len(G.nodes[node]["partition_candidates"]), len(G.nodes[neighbor]["partition_candidates"]))
            all_cands.append((G.nodes[node]["partition_candidates"], G.nodes[neighbor]["partition_candidates"]))
        else:
            ncand = None


        # GUS MOD
        print(':', len(G.nodes[node]["partition_candidates"]), "<>", len(G.nodes[neighbor]["partition_candidates"]), '=', end=' ')

        # If cost dict is provided, then add the cost of the new merged vertex, requires "partition_candidates" property
        if cost is not None:
            back_track.append(dict())
            back_track[-1][neighbor] = dict()
            if not is_directed and not bits:
                (new_cost, contrib, timing, numnotnone, granu) = joined_cost(G, G.nodes[node]["partition_candidates"], \
                    G.nodes[neighbor]["partition_candidates"], neighbor, node, neighbor, new_cost, back_track)
            elif not bits:
                (new_cost, contrib, timing, numnotnone, granu) = joined_cost_directed(G, G.nodes[node]["partition_candidates"], G.nodes[node]["parent_lists"], \
                    G.nodes[neighbor]["partition_candidates"], G.nodes[neighbor]["parent_lists"], neighbor, node, neighbor, new_cost, back_track)
            else:
                if not nodict:
                    new_cost = joined_cost_directed_2(G, neighbor, node, neighbor, new_cost, back_track)
                else:
                    (new_cost, contrib, timing, numnotnone, granu) = joined_cost_directed_no_dict(G, G.nodes[node]["partition_candidates"], G.nodes[node]["parent_lists"], \
                        G.nodes[neighbor]["partition_candidates"], G.nodes[neighbor]["parent_lists"], neighbor, node, neighbor, new_cost, back_track)
        flg = G.join_vertices(neighbor, node)
        G.nodes[neighbor]["merged"] = True

        print(len(G.nodes[neighbor]["partition_candidates"]))

        sumlogs[-1] += np.log(len(G.nodes[neighbor]["partition_candidates"]))
        max_cands.append(max(max_cands[-1], len(G.nodes[neighbor]["partition_candidates"])))

        if store_graphs:#G.number_of_nodes() < 20 and store_graphs:
            Gp = G.copy()
            graphs.append(Gp)

        QPut(G, next_pair, neighbor, max_degree)



        # if G.number_of_nodes() % 100 == 0 or G.number_of_nodes() < 100:
        #     print("Iteration", T, "with", G.number_of_nodes(), "nodes. Degree", ndegs, ncand)

        if ncand is not None and store_graphs:
            cand_prod.append(ncand[0] * ncand[1])
            times.append(time.time() - t0)

        T += 1

    v1 = list(G.nodes())[0]
    v2 = list(G.nodes())[1]
    cand_prod.append( len(G.nodes[v1]["partition_candidates"])*len(G.nodes[v2]["partition_candidates"]) )


    return (new_cost, contributions, back_track, times, cand_prod, degrees, ndegrees, timings, numnotnones, gran, all_cands, vjointimes, graphs, sumlogs, max_cands)

GLOBAL_CANDS = []
#def erosion_merge(G: PGraph, merge_seq: List[List[int]], max_degree: int = 2, max_iter=float("inf"), cost: dict=None, store_graphs=True, is_directed=False, bits=False, nodict=False):
def erosion_merge(G: PGraph, merge_seq: List[List[int]], cost: dict=None):#, is_directed=False, bits=False, nodict=False):
    global DIRECTED_TIMING
    global GLOBAL_CANDS
    for n in G:
        GLOBAL_CANDS.append(G.nodes[n]['parent_lists'])

    mlist_map = {n : n for n in G.nodes()}
    N = G.number_of_nodes()

    for node in G.nodes():
        G.nodes[node]["merged"] = False

    back_track = []

    if cost is not None: 
        new_cost = cost.copy()
    else: 
        new_cost = dict()

    # Begin main loop
    T = 1
    while T <= len(merge_seq) and G.number_of_nodes() > 2:#T < len(merge_seq):# and T <= max_iter and G.number_of_nodes() > 2:

        item = merge_seq[T-1]

        if item is None:
            break
        node, neighbor = mlist_map[item[0]], mlist_map[item[1]]

        if node not in G.nodes():
            continue

        # GUS MOD
        # if T-1 == 395:
        #     nc1 = sum([v.shape[1] for (k,v) in G.nodes[node]['parent_lists'].items()])
        #     nc2 = sum([v.shape[1] for (k,v) in G.nodes[neighbor]['parent_lists'].items()])
        #     print(node, '<>', neighbor, '\t:\t', nc1, "<>", nc2, '=', end=' ')

        # If cost dict is provided, then add the cost of the new merged vertex, requires "partition_candidates" property
        if cost is not None:
            back_track.append(dict())
            back_track[-1][neighbor] = dict()

            #new_cost = joined_cost_directed_2(G, neighbor, node, neighbor, new_cost, back_track)
            new_cost = joined_cost_directed_old(G, neighbor, node, neighbor, new_cost, back_track)

        
        flg = G.join_vertices(neighbor, node)


        mlist_map[N + T] = neighbor

        G.nodes[neighbor]["merged"] = True  # For some heuristics

        T += 1


    print('DIRECTED TIMING:', DIRECTED_TIMING)

    return new_cost, back_track


def calculate_optimal(G, costs):
    """ Calculate the optimal solution for an erroded G """
    assert G.number_of_nodes() == 2, "Graph not erroded"
    global debug_timings
    # t0 = time.time()
    
    best_cost = float("inf")
    
    u, v = list(G.nodes())
    bmap_u = G.nodes[u]["parent_dict"]
    bmap_v = G.nodes[v]["parent_dict"]
    found = False

    # Debug:
    # ninternal = len(set(bmap_u) & set(bmap_v))
    # ntotal = len(set(bmap_u) | set(bmap_v))
    
    for partition1 in G.nodes[u]["partition_candidates"]:
        for partition2 in G.nodes[v]["partition_candidates"]:
            # Check if connected component with no cycles:
            # t0 = time.time()

            # p1 = partition1.partition
            # p2 = partition2.partition
            # p1nparts =len(p1)
            # p2nparts =len(p2)
            # p1size = sum([len(pset) for pset in p1])
            # p2size = sum([len(pset) for pset in p2])

            

            # is_connected = partition1.connected(partition2)
            # t1 = time.time() - t0
            # t0 = time.time()
            # is_connected = partition1.connected_old(partition2)
            # t2 = time.time() - t0 
            # t0 = time.time()
            # is_connected = partition1.connected(partition2)
            # debug_timings.append((ninternal, ntotal, p1nparts, p2nparts, p1size, p2size, time.time() - t0, -1))
            # if not partition1.connected(partition2): continue

            # assert is_c == is_connected, str(p1.partition)+'+'+str(p2.partition)+'='+str(is_c)+' vs '+str(is_connected)


            # if is_connected:
            outbits_u = G.nodes[u]["parent_lists"][dp.get_encoding(partition1)]
            (R, C) = outbits_u.shape
            outbits_v = G.nodes[v]["parent_lists"][dp.get_encoding(partition2)]
            (R, D) = outbits_v.shape

            for c in range(C):
                for d in range(D):

                    #enc1 = dp.get_encoding(partition1, parents1)
                    enc1 = dp.get_encoding(partition1, parents=(bmap_u, outbits_u[:,c]), sep=False, np_array=True)
                    #enc2 = dp.get_encoding(partition2, parents2)
                    enc2 = dp.get_encoding(partition2, parents=(bmap_v, outbits_v[:,d]), sep=False, np_array=True)
                    
                    parents1 = (bmap_u, outbits_u[:,c])
                    parents2 = (bmap_v, outbits_v[:,d])

                    if type(costs[u][enc1]) == tuple:
                        sumcosts = costs[u][enc1][1] + costs[v][enc2][1]
                    else:
                        sumcosts = costs[u][enc1] + costs[v][enc2]

                    #assert partition1.connected(partition2) == partition1.connected_old(partition2), str(partition1.partition)+'vs'+str(partition2.partition)

                    if sumcosts < best_cost and \
                        dp.dconnected(partition1, partition2, parents1, parents2, bits=True)  and \
                           partition1.connected(partition2):
                        # t0 = time.time()
                        # is_connected = partition1.connected(partition2)
                        # debug_timings.append((ninternal, ntotal, p1nparts, p2nparts, p1size, p2size, time.time() - t0, -1))
                        # if not is_connected: continue
                        # if is_connected:
                        best_cost = sumcosts
                        G.nodes[u]["partition"] = partition1
                        G.nodes[u]["parents"] = outbits_u[:,c]
                        G.nodes[v]["partition"] = partition2
                        G.nodes[v]["parents"] = outbits_v[:,d]
                        found = True 

    # debug_timings.append(time.time() - t0)
    if not found:
        raise NoSolution()


def backtrack_erosion(graph, backtrack, directed=False, og_graph=None, save_graphs=False, bits=False, delta=False):

    if delta:
        global stats
        stats['chosen'] = []
        stats['all'] = []
    graphs = []  # Keep track of all intermediate steps

    # Create map from edge_name -> edge
    edgemap = dict()
    for edge in graph.edges(keys=True):
        edgemap[graph.edges[edge]["name"]] = edge

    # Go from last erosion iteration to the first in outer loop
    backtrack_rev = backtrack.copy()
    backtrack_rev.reverse()
    inner_iter=0
    outer_iter = 0
    for backtrack_t in backtrack_rev:

        #if outer_iter % 100 == 0:
        #print("Iteration", outer_iter, "of", len(backtrack_rev))

        for node in backtrack_t:
            inner_iter += 1

            if not directed:
                backtrack_partition_choice = backtrack_t[node][graph.nodes[node]["partition"].get_encoding()]
            else: # TODO: How to handle backtracking?
                chosen_partition = graph.nodes[node]["partition"]
                chosen_parents = graph.nodes[node]["parents"]
                if not bits:
                    enc = dp.get_encoding(chosen_partition, chosen_parents)
                else:
                    chosen_dict = graph.nodes[node]["parent_dict"]
                    #print("//....//")
                    #print((chosen_dict, chosen_parents))
                    enc = dp.get_encoding(chosen_partition, parents=(chosen_dict, chosen_parents), sep=False, np_array=True)
                
                # GET STATS
                if delta:
                    #print(len(stats['per join cost']), '--', -outer_iter - 1)
                    costs = stats['per join cost'][-outer_iter - 1]
                    stats['chosen'].append(costs[enc])
                    stats['all'].append([costs[key] for key in costs])

                backtrack_partition_choice = backtrack_t[node][enc]

            if not directed:
                node_previous_partition = [backtrack_partition_choice[n] for n in backtrack_partition_choice if n == node][0]
            else:
                if bits:
                    (node_previous_partition, node_previous_parents, node_dict) = [backtrack_partition_choice[n] for n in backtrack_partition_choice if n == node][0]
                else:
                    (node_previous_partition, node_previous_parents) = [backtrack_partition_choice[n] for n in backtrack_partition_choice if n == node][0]
            other_previous_node = [n for n in backtrack_partition_choice if n != node][0]
            if not directed:
                other_previous_partition = [backtrack_partition_choice[n] for n in backtrack_partition_choice if n != node][0]
            else: #####################
                if bits:
                    (other_previous_partition, other_previous_parents, other_dict) = [backtrack_partition_choice[n] for n in backtrack_partition_choice if n != node][0]
                else:
                    (other_previous_partition, other_previous_parents) = [backtrack_partition_choice[n] for n in backtrack_partition_choice if n != node][0]

            new_edge_names = {e for e in other_previous_partition.fullset if e not in edgemap}
            other_previous_edges = [edgemap[e] for e in other_previous_partition.fullset if e not in new_edge_names]

            graph.add_node(other_previous_node, partition=other_previous_partition)
            if directed:
                #print("Parents 0: ", np.array(other_previous_parents[1]))
                if bits:
                    graph.nodes[other_previous_node]["parents"] = np.array(other_previous_parents[1])
                    graph.nodes[other_previous_node]["parent_dict"] = other_dict
                else:
                    graph.nodes[other_previous_node]["parents"] = other_previous_parents
            if og_graph is not None:
                for key in og_graph.nodes[other_previous_node]:
                    if key not in graph.nodes[other_previous_node]:
                        graph.nodes[other_previous_node][key] = og_graph.nodes[other_previous_node][key]

            for (u, v, k) in other_previous_edges:
                if u != node:
                    u, v = v, u  # make it so that u = node, v = other thing

                name = graph.edges[(u, v, k)]["name"]

                k2 = graph.move_edge((u, v, k), (v, other_previous_node))

                edgemap[name] = (v, other_previous_node, k2)

            for e in new_edge_names:
                k3 = graph.add_edge(node, other_previous_node, name=e)
                edgemap[e] = (node, other_previous_node, k3)
            
            graph.nodes[node]["partition"] = node_previous_partition
            if directed:
                #print("Parents 0 (2): ", np.array(other_previous_parents[1]))
                if bits:
                    graph.nodes[node]["parents"] = np.array(node_previous_parents[1])
                    graph.nodes[node]["parent_dict"] = node_dict
                else:
                    graph.nodes[node]["parents"] = node_previous_parents
            if og_graph is not None:
                for key in og_graph.nodes[node]:
                    if key not in graph.nodes[node]:
                        graph.nodes[node][key] = og_graph.nodes[node][key]

            if save_graphs:
                graphs.append((graph.copy(), []))
        outer_iter += 1

    return graphs


#seed(3)
#np.random.seed(3)
            
def brute_force(G, cost, first_viable=False):
    """ Covers all possible combinations of partitions """
    current_best_cost = float('inf')
    current_best_choice = []
    partition_candidates = [G.nodes[n]["partition_candidates"] for n in G.nodes()]
    node_order = [n for n in G.nodes()]
    all_comb = product(*partition_candidates)

    #print("We have to iterate over", len(list(product(*partition_candidates))), "choices.")
    
    for choice in all_comb:
        it = 0
        for partition in choice:
            G.nodes[node_order[it]]["partition"] = partition
            it += 1
        
        # If a tree, then update the current best cost and choice of partitions
        G_part = G.partition_graph()
        if G_part.number_of_nodes() - G_part.number_of_edges() == 1 and nx.number_connected_components(G_part) == 1:
            #print("Node", n, "encoding", G.nodes[n]["partition"].get_encoding(), "for partition", G.nodes[n]["partition"].partition)
            if first_viable:
                current_best_cost = 0
                current_best_choice = choice
                break
            total_cost = sum([cost[n][G.nodes[n]["partition"].get_encoding()] for n in G.nodes()])
            if total_cost < current_best_cost:
                current_best_cost = total_cost
                current_best_choice = choice
    
    # Update G with the best choice
    it = 0
    for partition in current_best_choice:
        G.nodes[node_order[it]]["partition"] = partition
        it += 1

    return current_best_cost, current_best_choice


def generate_test_graph(init_points, max_verts = float("inf"), random_seed=None, width_mult=1, directed=False):
    
    if random_seed is not None:
        seed(random_seed)
        np.random.seed(random_seed)

    # Generate triangles
    num_triangles = float("inf")
    while num_triangles >= max_verts:
        points = np.random.normal(3, 2.5, size=(init_points, 2))
        points[:,1] = points[:,1]*width_mult

        xmap=dict()
        ymap=dict()
        for i in range(points[:,1].size):
            xmap[i+1] = points[i,0]
            ymap[i+1] = points[i,1]

        tri = Delaunay(points)
        triangles = tri.simplices
        num_triangles = triangles.shape[0]

    triangle_graph = nx.Graph()

    next_edge_name = 1
    edges_map = dict()
    for i in range(num_triangles):
        u = triangles[i, 0]
        v = triangles[i, 1]
        w = triangles[i, 2]

        if u not in triangle_graph.nodes():
            triangle_graph.add_node(u, x=points[u,0], y=points[u,1])
        if v not in triangle_graph.nodes():
            triangle_graph.add_node(v, x=points[v,0], y=points[v,1])
        if w not in triangle_graph.nodes():
            triangle_graph.add_node(w, x=points[w,0], y=points[w,1])
        if (u, v) not in triangle_graph.edges():
            triangle_graph.add_edge(u, v)
        if (v, w) not in triangle_graph.edges():
            triangle_graph.add_edge(v, w)
        if (w, u) not in triangle_graph.edges():
            triangle_graph.add_edge(w, u)

        # Get centroid:
        xmap[i+1] = (points[u,0] + points[v,0] + points[w,0])/3
        ymap[i+1] = (points[u,1] + points[v,1] + points[w,1])/3

        if (u, v) not in edges_map: 
            edges_map[(u, v)] = [i]
            edges_map[(v, u)] = [i]
        else: 
            edges_map[(u, v)].append(i)
            edges_map[(v, u)].append(i)
        if (v, w) not in edges_map: 
            edges_map[(w, v)] = [i]
            edges_map[(v, w)] = [i]
        else: 
            edges_map[(w, v)].append(i)
            edges_map[(v, w)].append(i)
        if (u, w) not in edges_map: 
            edges_map[(u, w)] = [i]
            edges_map[(w, u)] = [i]
        else: 
            edges_map[(u, w)].append(i)
            edges_map[(w, u)].append(i)

    edges = []
    seen = set()
    for edge in edges_map:
        if len(edges_map[edge]) == 2:
            t1 = edges_map[edge][0]+1
            t2 = edges_map[edge][1]+1

            if (t2, t1) not in seen:
                edges.append((t1, t2, next_edge_name))
                next_edge_name += 1
                seen.add((t2, t1))
        

    # Costs and partitions
    G = PGraph()

    for (u, v, k) in edges:
        G.add_node(u, x=xmap[u], y=ymap[u])
        G.add_node(v, x=xmap[v], y=ymap[v])
        G.add_edge(u, v, name=k)

    # This time, assign all partitions and their respective costs:
    cost = dict()
    for node in G.nodes():
        edge_names = [G.edges[e]["name"] for e in G.edges(node, keys=True)]
        G.nodes[node]["partition_candidates"] = prt.Partition.all_partitions(edge_names)
        G.nodes[node]["partition"] = G.nodes[node]["partition_candidates"][len(G.nodes[node]["partition_candidates"])-1]

        # Create parent list
        if directed:
            if node != 1:
                G.nodes[node]["parent_lists"] = dp.generate_directed(G.nodes[node]["partition_candidates"])
            else:
                G.nodes[node]["parent_lists"] = dp.generate_directed(G.nodes[node]["partition_candidates"], True)
            cost[node] = dp.generate_costs(G.nodes[node]["partition_candidates"], \
                G.nodes[node]["parent_lists"])
            
            G.nodes[node]["num_cands"] = sum([len(G.nodes[node]["parent_lists"][key]) for key in G.nodes[node]["parent_lists"]])
            G.nodes[node]["num_theo_cands"] = number_directed_partitions(nx.degree(G, node))
        else:
            cost[node] = dict()
            for part in G.nodes[node]["partition_candidates"]:
                cost[node][part.get_encoding()] = random()
    
    return G, cost, triangle_graph







# BIT JOINS
def partition_convert(partition, parentlists, bitdict0=None):
    """ Temp convert function """
    edges = list(partition.fullset)
    if bitdict0 is None:
        bitdict = {edges[i] : i for i in range(len(edges))}
    else:
        bitdict = bitdict0
    bitmaps = -np.ones((len(edges), len(parentlists)))
    for i in range(len(parentlists)):
        parentlist = parentlists[i]
        if None in parentlist: parentlist.remove(None)
        bitmaps[[bitdict[e] for e in parentlist], i] = 1
    if bitdict0 is not None:
        return bitmaps
    else:
        return (bitdict, bitmaps)

def partition_revert(partition, diction, bitmap):
    (R, C) = bitmap.shape 

    parentlists = []
    for c in range(C):
        parentlist = []
        # Not very viable WHERE is the parent? In which set?


import dpartitions2 as dp2 

class NoSolution(Exception):
    def __init__(self, message="No solution was found with the given beam width and partition candidates."):
        self.message = message


def add_costs(c1, c2):
    if type(c1) == tuple:
        return (not( (not c1[0]) and (not c2[0]) ), c1[1] + c2[1])
    else:
        return c1 + c2


debug_timings = []

# def joined_cost_directed_2(G, partition_candidates_1, pcands1, partition_candidates_2, pcands2, new_node, n1, n2, cost, back_track):
def joined_cost_directed_old(G, new_node, n1, n2, cost, back_track):
    """ Given partition candidates and a cost function, what is the cost of the new joined partition? """

    global params
    # global debug_timings # DEBUG CODE
    
    seen = dict()

    # contributions=dict()
    # contributions[-1] = (new_node, n1, n2)  # "metadata"

    new_cost = cost.copy()
    new_cost[new_node] = dict()  # Restart this node

    new_part_candidates = BeamList(params)

    pl_bits = dict()
    pl_bits[n1] = G.nodes[n1]["parent_lists"]
    pl_bits[n2] = G.nodes[n2]["parent_lists"]
    bdict1 = G.nodes[n1]["parent_dict"]
    bdict2 = G.nodes[n2]["parent_dict"]


    d_obj = dp2.generate_intersection_map(G.nodes[n1]["partition_candidates"], bdict1, 
                             G.nodes[n2]["partition_candidates"], bdict2)

    

    # (3) Undirected join, and extract directed join object
    ######################1
    (D1ord, D2ord) = d_obj["ords"]
    (border_map_1, border_edges_1, border_map_2, border_edges_2) = d_obj["border_edges"]

    num_in = d_obj["num_edges"]
    # outbits_dict = {}
    bmap = {border_edges_1[i] : i for i in range(len(border_edges_1))}
    bmap = {**bmap, **{border_edges_2[i] : i + len(border_edges_1) for i in range(len(border_edges_2))}}

    #print("Parent Map:", bmap)
    
    ######################2
    j_obj = []
    cands1 = G.nodes[n1]["partition_candidates"]
    cands2 = G.nodes[n2]["partition_candidates"]

    if n1 == G.graph['root']:
        root_edges = G.graph['root_edges']
    elif n2 == G.graph['root']:
        root_edges = G.graph['root_edges']
    else:
        root_edges=None

    d_obj['bmap'] = bmap
    d_obj['border-dict'] = {e:1 for e in d_obj['border']}
    for e in d_obj['non-border']: d_obj['border-dict'][e] = 0

    (outd, bmap) = dp2.direction_join(j_obj, d_obj)
    
    # t1 = 0
    # t2 = 0
    # t3 = 0
    # t4 = 0
    global DIRECTED_TIMING
    
    for partition1 in cands1:#partition_candidates_1: NEW MOD
        uenc1 = dp.get_encoding(partition1)

        t0 = get_time()

        bitmap1 = pl_bits[n1][uenc1] # DBUG
        bitmap1 = np.array(bitmap1)

        DIRECTED_TIMING += get_time() - t0 

        for partition2 in cands2:#partition_candidates_2:
            

            t0 = get_time()
            # DIRECTED CHECK
            #################1 multiplication operation and outbits obtain
            # t0 = time.time()
            uenc2 = dp.get_encoding(partition2)
            bitmap2 = pl_bits[n2][uenc2]
            #partition_encoding = ujoin.get_encoding()
            
            # Temp: not sure why converting to np.array
            bitmap2 = np.array(bitmap2)
            R = bitmap1[D1ord, :].T @ bitmap2[D2ord, :]
            (p1chosen, p2chosen) = np.where(R == -num_in)  # each of these indices corresponds to a cost
            # if p1chosen.size == 0: continue 
            chosenD1 = bitmap1[:, p1chosen]
            chosenD2 = bitmap2[:, p2chosen]

            outbits1 = chosenD1[border_map_1, :]
            outbits2 = chosenD2[border_map_2, :]
            outbits = np.vstack((outbits1, outbits2))

            DIRECTED_TIMING += get_time() - t0

            if outbits.size == 0: continue

            # # # # END CYCLE CHECK
            

            # # DEBUG
            # if ujoin is None or outbits.size == 0:
            #     enc1 = None
            # else:
            #     enc1 = ujoin.get_encoding()
            # if ujoin2 is None or s0 == 0:
            #     enc2 = None
            # else:
            #     enc2 = ujoin2.get_encoding()

            # #assert enc1 == enc2, str(partition1.partition) + ' + ' + str(partition2.partition)
            # if enc1 != enc2:
            #     print('Stop')


            #outbits_dict[partition_encoding] = outbits
            ######################2 cost check


            # ENCODING MATCH
            # t0 = time.time()
            
            # if s0 != 0:
            # t0 = time.time()
            ujoin = dp.ujoin(partition1, partition2, d_obj,root_edges=root_edges)#, params=params) # DS Code
            # debug_timings.append((ninternal, ntotal, p1nparts, p2nparts, p1size, p2size, time.time() - t0, -1))
            # t2 += time.time() - t0
            if ujoin is None: continue
            # if ujoin is not None: # fill join object
            # t0 = time.time()
            uenc_join = dp.get_encoding(ujoin) 

            
            (R, C) = outbits.shape

            for c in range(C):
                # Renaming for similarity
                parents1 = (bdict1, chosenD1[:, c])
                parents2 = (bdict2, chosenD2[:, c])
                # join_undirected, join_parents = ujoin, outbits[:,c]

                join = (ujoin, outbits[:,c])

                # t02 = time.time()

                join_ind = "[" + uenc_join + ", " +  dp.get_encoding(ujoin, parents=(bmap, outbits[:,c]), sep=False, np_array=True, und=False) + "]"
                p1_ind = "[" + uenc1 + ", " + dp.get_encoding(partition1, parents1, sep=False, np_array=True, und=False) + "]"
                p2_ind = "[" + uenc2 + ", " + dp.get_encoding(partition2, parents2, sep=False, np_array=True, und=False) + "]"
                
                # time_getenc += time.time() - t02

                if join_ind not in seen:
                    #seen.add(join_ind)
                    seen[join_ind] = 1
                    back_track[-1][new_node][join_ind] = (-1, -1)

                # hash -> O(1)
                if new_node not in new_cost:
                    new_cost[new_node] = dict()
                
                # hash -> O(1)
                if join_ind not in new_cost[new_node]:
                    # if type(cost[n1][p1_ind]) != tuple:
                    new_cost[new_node][join_ind] = float("inf")
                    # else:
                    #     new_cost[new_node][join_ind] = (float("inf"), float("inf"))

                #print(new_cost[new_node][join_ind], join, join_ind)
                if new_cost[new_node][join_ind] > add_costs(cost[n1][p1_ind], cost[n2][p2_ind]):
                    new_cost[new_node][join_ind] = add_costs(cost[n1][p1_ind], cost[n2][p2_ind])
                    back_track[-1][new_node][join_ind] = {n1:(partition1, parents1, bdict1), n2:(partition2, parents2, bdict2)}
                    new_part_candidates.add(new_cost[new_node][join_ind], join, join_ind)
            # t3 += time.time() - t0

    # t0 = time.time()
    # Execute join ops
    (outd, bmap) = dp2.direction_join(j_obj, d_obj)

    # Persistent Partition
    if 'greedy_partition' in G.nodes[n1] and 'greedy_partition' in G.nodes[n2] \
        and 'greedy_parents' in G.nodes[n1] and 'greedy_parents' in G.nodes[n2]:
        # Rather than making n^2 comparisons, we can do this once at the cost
        # of a redundant partition join
       
        #ujoin = dp.ujoin(G.nodes[n1]['greedy_partition'], G.nodes[n2]['greedy_partition'], d_obj, params=params)
        ujoin = dp.ujoin(G.nodes[n1]['greedy_partition'], G.nodes[n2]['greedy_partition'], d_obj, params=params)
        # ninternal = len([e for (e, b) in d_obj['border-dict'].items() if b==1])
        # ntotal = len(d_obj['border-dict'])
        # p1 = G.nodes[n1]['greedy_partition'].partition
        # p2 = G.nodes[n2]['greedy_partition'].partition
        # p1nparts =len(p1)
        # p2nparts =len(p2)
        # p1size = sum([len(pset) for pset in p1])
        # p2size = sum([len(pset) for pset in p2])
        # debug_timings.append((ninternal, ntotal, p1nparts, p2nparts, p1size, p2size, t1, t2))
        
        # bmap is dictionary
        # if dp.ujoin(G.nodes[n1]['greedy_partition'], G.nodes[n2]['greedy_partition'], d_obj, params=params) is None:
        #     print("None")
            
        chosenD1 = G.nodes[n1]['greedy_parents']
        chosenD2 = G.nodes[n2]['greedy_parents']
        outbits1 = chosenD1[border_map_1]
        outbits2 = chosenD2[border_map_2]

        outbits = np.hstack((outbits1, outbits2))
        # if ujoin is None:
        #     2+2
        #     ujoin = dp.ujoin(G.nodes[n1]['greedy_partition'], G.nodes[n2]['greedy_partition'], d_obj, params=params)
    
        enc = dp.get_encoding(ujoin, (bmap, outbits), np_array=True)

        # enc1 = dp.get_encoding(G.nodes[n1]['greedy_partition'], (G.nodes[n1]['parent_dict'], outbits1), np_array=True)
        # enc2 = dp.get_encoding(G.nodes[n2]['greedy_partition'], (G.nodes[n2]['parent_dict'], outbits2), np_array=True)
        # gred_cost = add_costs(cost[n1][enc1], cost[n2][enc2])

        # print(partition1.partition, '+', partition2.partition, '=', ujoin.partition)print
        # p1 = sorted([(e, chosenD1[bdict1[e]]) for e in border_edges_1])
        # p2 = sorted([(e, chosenD2[bdict2[e]]) for e in border_edges_2])
        # re = sorted([(e, outbits[bmap[e]]) for e in border_edges_1 + border_edges_2])
        # print(p1, '+', p2, '=', re)
        G.nodes[new_node]['greedy_partition'] = ujoin
        G.nodes[new_node]['greedy_parents'] = outbits

        # if n1 == 588 or n2 == 588:
        #     print()
        #     print('n1', n1, ', n2', n2, (nx.degree(G,n1), nx.degree(G, n2)))
        #     print(G.nodes[n1]['greedy_partition'].partition, '+', G.nodes[n2]['greedy_partition'].partition, '\n\t=', ujoin.partition)
        #     print()
        unzipped = list(zip(*new_part_candidates.get_list( enc )))
    else:
        unzipped = list(zip(*new_part_candidates.get_list()))

    if unzipped == []:
        raise NoSolution()
    partcands, parentcands = unzipped

    # if params.histogram:
    #     #global stats
    #     #print('----', new_part_candidates.per_join_cost)
    #     stats['per join cost'].append(new_part_candidates.per_join_cost)

    partcands = list(partcands)
    parentcands = list(parentcands)
    parentcands_dict = {dp.get_encoding(cand):[] for cand in partcands}
    new_part_cands = []  # Only store unique in here
    seen_cand = set()
    for i in range(len(partcands)):
        #sorted_part, sorted_par = dp.sort_partition(partcands[i], parentcands[i])
        sorted_part, sorted_par, root_ind = dp.sort_partition(partcands[i], (bmap, parentcands[i]), True, partcands[i].root)
        part_cand_enc = dp.get_encoding(partcands[i])
        #parentcands_dict[part_cand_enc].append(sorted_par)
        parentcands_dict[part_cand_enc].append(parentcands[i])
        #partcands[i] = prt.Partition(sorted_part)  # Replace partition with sorted version
        if part_cand_enc not in seen_cand:
            partition = prt.Partition(sorted_part)
            partition.root = root_ind
            new_part_cands.append(partition)
            seen_cand.add(part_cand_enc)
        # else:
        #     raise Exception('Debug Exception')

    parentcands_dict = {key : np.array(parentcands_dict[key]).T for key in parentcands_dict}
    # print(parentcands)

    G.nodes[new_node]["partition_candidates"] = new_part_cands#partcands
    G.nodes[new_node]["parent_lists"] = parentcands_dict
    #G.nodes[new_node]["num_cands"] = sum([len(parentcands_dict[key]) for key in parentcands_dict])
    G.nodes[new_node]["num_undirected"] = len(new_part_cands)
    # G.nodes[new_node]["num_join_ops"] = num_join_ops
    G.nodes[new_node]["parent_dict"] = bmap

    # t4 = time.time() - t0

    # debug_timings.append((t1, t2, t3, t4))

    return new_cost#, timings, numnotnones, all_granulated_timings)


def joined_cost_directed_2(G, new_node, n1, n2, cost, back_track):
    """ Given partition candidates and a cost function, what is the cost of the new joined partition? """

    global params
    global debug_timings

    new_part_candidates = BeamList(params)
    back_track, new_cost, d_obj, bmap, border_map_1, border_map_2, timings = dp2.construct_matrices(G, n1, n2, new_node, back_track, cost, new_part_candidates)

    debug_timings.append(timings)
    #(outd, bmap) = dp2.direction_join([], d_obj)

    # Persistent Partition
    if 'greedy_partition' in G.nodes[n1] and 'greedy_partition' in G.nodes[n2] \
        and 'greedy_parents' in G.nodes[n1] and 'greedy_parents' in G.nodes[n2]:

        ujoin = dp.ujoin(G.nodes[n1]['greedy_partition'], G.nodes[n2]['greedy_partition'], d_obj, params=params)
            
        chosenD1 = G.nodes[n1]['greedy_parents']
        chosenD2 = G.nodes[n2]['greedy_parents']
        outbits1 = chosenD1[border_map_1]
        outbits2 = chosenD2[border_map_2]

        outbits = np.hstack((outbits1, outbits2))
    
        enc = dp.get_encoding(ujoin, (bmap, outbits), np_array=True)
        G.nodes[new_node]['greedy_partition'] = ujoin
        G.nodes[new_node]['greedy_parents'] = outbits
        unzipped = list(zip(*new_part_candidates.get_list( enc )))
    else:
        unzipped = list(zip(*new_part_candidates.get_list()))

    if unzipped == []:
        raise NoSolution()
    partcands, parentcands = unzipped

    partcands = list(partcands)
    parentcands = list(parentcands)
    parentcands_dict = {dp.get_encoding(cand):[] for cand in partcands}
    new_part_cands = []  # Only store unique in here
    seen_cand = set()
    for i in range(len(partcands)):
        sorted_part, sorted_par, root_ind = dp.sort_partition(partcands[i], (bmap, parentcands[i]), True, partcands[i].root)
        part_cand_enc = dp.get_encoding(partcands[i])
        parentcands_dict[part_cand_enc].append(parentcands[i])

        if part_cand_enc not in seen_cand:
            partition = prt.Partition(sorted_part)
            partition.root = root_ind
            new_part_cands.append(partition)
            seen_cand.add(part_cand_enc)

    parentcands_dict = {key : np.array(parentcands_dict[key]).T for key in parentcands_dict}
    
    G.nodes[new_node]["partition_candidates"] = new_part_cands
    G.nodes[new_node]["parent_lists"] = parentcands_dict
    G.nodes[new_node]["num_undirected"] = len(new_part_cands)
    G.nodes[new_node]["parent_dict"] = bmap

    return new_cost





######################################
# Trying it without the dictionaries #
######################################
time_sort = 0
def joined_cost_directed_no_dict(G, partition_candidates_1, pcands1, partition_candidates_2, pcands2, new_node, n1, n2, cost, back_track):
    """ Given partition candidates and a cost function, what is the cost of the new joined partition? """
    global time_mult
    global time_enc
    #global NUMBER_MINS
    global params
    global time_sort

    #seen = set()
    seen = dict()

    contributions=dict()
    contributions[-1] = (new_node, n1, n2)  # "metadata"
    new_cost = cost.copy()
    new_cost[new_node] = dict()  # Restart this node
    new_part_candidates = BeamList(params)

    timings = []
    gran_times = []
    numnotnones = 0

    all_granulated_timings = []
    num_join_ops = 0

    cost_l = []

    # (2) Extract direction object
    pl_bits = dict()
    pl_bits[n1] = G.nodes[n1]["parent_lists"]
    pl_bits[n2] = G.nodes[n2]["parent_lists"]
    bdict1 = G.nodes[n1]["parent_dict"]
    bdict2 = G.nodes[n2]["parent_dict"]

    d_obj = dp2.generate_intersection_map(G.nodes[n1]["partition_candidates"], bdict1, 
                             G.nodes[n2]["partition_candidates"], bdict2)

    # (3) Undirected join, and extract directed join object
    ######################1
    (D1ord, D2ord) = d_obj["ords"]
    (border_map_1, border_edges_1, border_map_2, border_edges_2) = d_obj["border_edges"]

    num_in = d_obj["num_edges"]
    outbits_dict = {}
    bmap = {border_edges_1[i] : i for i in range(len(border_edges_1))}
    bmap = {**bmap, **{border_edges_2[i] : i + len(border_edges_1) for i in range(len(border_edges_2))}}

    #print("Parent Map:", bmap)
    
    ######################2
    j_obj = []
    for partition1 in partition_candidates_1:
        for partition2 in partition_candidates_2:
            # Faster to check undirected first:
            ujoin = dp.ujoin(partition1, partition2, d_obj, params=params)
            #if ujoin is None: continue
            if ujoin is not None: # fill join object
                # j_obj.append( (ujoin.get_encoding(),\
                #     pl_bits[partition1.get_encoding()],\
                #     pl_bits[partition2.get_encoding()] ))
                
                #################1 multiplication operation and outbits obtain
                t0 = time.time()
                bitmap1 = pl_bits[n1][dp.get_encoding(partition1)]
                bitmap2 = pl_bits[n2][dp.get_encoding(partition2)]
                partition_encoding = ujoin.get_encoding()
                
                # Temp: not sure why...
                bitmap1 = np.array(bitmap1)
                bitmap2 = np.array(bitmap2)
                R = bitmap1[D1ord, :].T @ bitmap2[D2ord, :]
                (p1chosen, p2chosen) = np.where(R == -num_in)  # each of these indices corresponds to a cost
                chosenD1 = bitmap1[:, p1chosen]
                chosenD2 = bitmap2[:, p2chosen]

                outbits1 = chosenD1[border_map_1, :]
                outbits2 = chosenD2[border_map_2, :]
                outbits = np.vstack((outbits1, outbits2))

                outbits_dict[partition_encoding] = outbits
                ######################2 cost check
                (R, C) = outbits.shape
                time_mult += time.time() - t0
                t0 = time.time()

                # Match encodings and compare costs:
                for c in range(C):
                    # Renaming for similarity
                    parents1 = (bdict1, chosenD1[:, c])
                    parents2 = (bdict2, chosenD2[:, c])
                    join_undirected, join_parents = ujoin, outbits[:,c]

                    join = (ujoin, outbits[:,c])

                    join_ind = dp.get_encoding(ujoin, parents=(bmap, outbits[:,c]), sep=False, np_array=True)
                    p1_ind = dp.get_encoding(partition1, parents1, sep=False, np_array=True)
                    p2_ind = dp.get_encoding(partition2, parents2, sep=False, np_array=True)

                    cost_l.append((join_ind, cost[n1][p1_ind] + cost[n2][p2_ind], \
                        (partition1, parents1, bdict1), (partition2, parents2, bdict2), join, join_ind))

                    # if join_ind not in seen:
                    #     #seen.add(join_ind)
                    #     seen[join_ind] = 1
                    #     back_track[-1][new_node][join_ind] = (-1, -1)

                    # # hash -> O(1)
                    # if new_node not in new_cost:
                    #     new_cost[new_node] = dict()
                    
                    # # hash -> O(1)
                    # if join_ind not in new_cost[new_node]:
                    #     new_cost[new_node][join_ind] = float("inf")

                    # if new_cost[new_node][join_ind] > cost[n1][p1_ind] + cost[n2][p2_ind]:
                    #     new_cost[new_node][join_ind] = cost[n1][p1_ind] + cost[n2][p2_ind]
                    #     back_track[-1][new_node][join_ind] = {n1:(partition1, parents1, bdict1), n2:(partition2, parents2, bdict2)}
                    #     new_part_candidates.add(new_cost[new_node][join_ind], join, join_ind)

                    #timings.append(time.time() - t0)
                time_enc += time.time() - t0
    
    # Resolve duplicates:
    t0 = time.time()
    cost_l.sort()
    chosen_costs = []

    for i in range(len(cost_l) - 1):
        
        if cost_l[i][0] == cost_l[i+1][0] and cost_l[i][1] < cost_l[i+1][1]:
            cost_l[i+1] = cost_l[i]
        elif cost_l[i][0] != cost_l[i+1][0]:
            new_part_candidates.add(cost_l[i][1], cost_l[i][4], cost_l[i][5])
            chosen_costs.append(cost_l[i])

    new_part_candidates.add(cost_l[-1][1], cost_l[-1][4], cost_l[-1][5])
    chosen_costs.append(cost_l[-1])

    time_sort += time.time() - t0

    # Put back to format:
    new_cost[new_node] = {
        item[0] : item[1] for item in chosen_costs
    }
    back_track[-1][new_node] = {
        item[0] : {n1:item[2], n2:item[3]} for item in chosen_costs
    }


    # Execute join ops
    (outd, bmap) = dp2.direction_join(j_obj, d_obj)

    # Question: How to backtrack and handle cost function? 

    #print(new_part_candidates.get_list())
    partcands, parentcands = zip(*new_part_candidates.get_list())
    partcands = list(partcands)
    parentcands = list(parentcands)
    parentcands_dict = {dp.get_encoding(cand):[] for cand in partcands}
    new_part_cands = []  # Only store unique in here
    seen_cand = set()
    for i in range(len(partcands)):
        #sorted_part, sorted_par = dp.sort_partition(partcands[i], parentcands[i])
        sorted_part, sorted_par = dp.sort_partition(partcands[i], (bmap, parentcands[i]), True)
        part_cand_enc = dp.get_encoding(partcands[i])
        #parentcands_dict[part_cand_enc].append(sorted_par)
        parentcands_dict[part_cand_enc].append(parentcands[i])
        #partcands[i] = prt.Partition(sorted_part)  # Replace partition with sorted version
        if part_cand_enc not in seen_cand:
            new_part_cands.append(prt.Partition(sorted_part))
            seen_cand.add(part_cand_enc)

    parentcands_dict = {key : np.array(parentcands_dict[key]).T for key in parentcands_dict}
    # print(parentcands)

    G.nodes[new_node]["partition_candidates"] = new_part_cands#partcands
    G.nodes[new_node]["parent_lists"] = parentcands_dict
    #G.nodes[new_node]["num_cands"] = sum([len(parentcands_dict[key]) for key in parentcands_dict])
    G.nodes[new_node]["num_undirected"] = len(new_part_cands)
    G.nodes[new_node]["num_join_ops"] = num_join_ops
    G.nodes[new_node]["parent_dict"] = bmap

    return (new_cost, contributions, timings, numnotnones, all_granulated_timings)


def compute_cost(beampg, cost_dict, pname='parents', all_vertices=False):
    # cb = nx.cycle_basis(nx.Graph(beampg))
    # cycle_vertices = set()
    # for l in cb:
    #     cycle_vertices = cycle_vertices | set(l)
    if all_vertices:
        verts = [n for n in beampg]
    else:
        verts = [n for n in beampg if nx.degree(beampg, n) > 2 and n != beampg.graph['root']]
    beam_cost = 0
    for node in verts:
        parbits = beampg.nodes[node][pname]
        bdict = beampg.nodes[node]["parent_dict"]
        beampg.nodes[node]["cost"] = cost_dict[node][dp.get_encoding(beampg.nodes[node]["partition"], (bdict, parbits), sep=False, np_array=True)]
        beam_cost += beampg.nodes[node]["cost"]

    return beam_cost


def partition_convert(partition, parentlists, bitdict0=None):
    """ Convert integer partition to bit partition """
    edges = list(partition.fullset)
    if bitdict0 is None:
        bitdict = {edges[i] : i for i in range(len(edges))}
    else:
        bitdict = bitdict0
    bitmaps = -np.ones((len(edges), len(parentlists)))
    for i in range(len(parentlists)):
        parentlist = parentlists[i]
        if None in parentlist: 
            #ind = parentlist.index(None)
            #new_parents = partition.partition[ind]
            parentlist.remove(None) 
            #parentlist += new_parents
            
        bitmaps[[bitdict[e] for e in parentlist], i] = 1
    if bitdict0 is not None:
        return bitmaps
    else:
        return (bitdict, bitmaps)