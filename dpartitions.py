import partitions as prt
import random
from itertools import product
import time

# Partition, Parents are lists of the same size

def u_is_joinable(p1: prt.Partition, parents_1, p2: prt.Partition, parents_2):
    if p1.join(p2) is None:
        #print("Undirected check fail")
        return False

def d_is_joinable(p1: prt.Partition, parents_1, p2: prt.Partition, parents_2):
    all_items_1 = {edge for pset in p1.partition for edge in pset}
    all_items_2 = {edge for pset in p2.partition for edge in pset}
    children_1 = all_items_1 - set(parents_1)
    children_2 = all_items_2 - set(parents_2)
    # Children intersect check:
    if (children_1 & children_2) != set():  # Check if set not empty
        #print("Children check fail")
        return False
    elif (set(parents_1) & set(parents_2)) != set():
        #print("Parent check fail")
        return False
    else:
        return True


def ujoin(p1: prt.Partition, p2: prt.Partition, d_obj=None, params=None, root_edges=None):
    # if len(p1.partition) + len(p2.partition) > 30:
    #     return p1.join(p2, d_obj, root_edges)#, params=params)#  DS Code
    # else:
    #     return p1.join_old(p2, d_obj, params=params) # DS Code
    # t0 = time.time()
    
    # time1 = time.time() - t0
    # t1 = time.time()

    j1 = p1.join(p2, d_obj)#, root_edges)
    # j1 = p1.join_old(p2, d_obj, params=params)
    # j1 = p1.join_old(p2, d_obj)#, root_edges)

    # # time2 = time.time() - t1

    # if j1 is None or j2 is None:
    #     assert j1 == j2, str(p1.partition)+'+'+str(p2.partition)+' - Nones not equal'
    # else:
    #     #print(j1.partition, j2.partition)
    #     assert j1.get_encoding() == j2.get_encoding(), str(p1.partition)+'+'+str(p2.partition)+'='+str(j1.partition)+' vs '+str(j2.partition)

    # time1 = 0
    # time2 = 0
    return j1

def djoin(undirected_join: prt.Partition, p1: prt.Partition, p2: prt.Partition, parents_1, parents_2):
    """
    The directed poriton of the join operator. This function is not used when
    using bits (numpy arrays)
    """
    #print(parents_1, "*", parents_2, " -> ")
    #print("     ", is_joinable(p1, parents_1, p2, parents_2))
    #if None in parents_1 or None in parents_2:
    #    print(parents_1, "*", parents_2, " -> ", is_joinable(p1, parents_1, p2, parents_2))

    all_items_1 = {edge for pset in p1.partition for edge in pset}
    if None in parents_1: 
        all_items_1.add(None)
    all_items_2 = {edge for pset in p2.partition for edge in pset}
    if None in parents_2: 
        all_items_2.add(None)
    children_1 = all_items_1 - set(parents_1)
    children_2 = all_items_2 - set(parents_2)
    in_items = all_items_1 & all_items_2
    out_items = (all_items_1 | all_items_2) - in_items

    if (children_1 & children_2) != set():  # Check if set not empty
        #print("Children check fail")
        return None
    elif (set(parents_1) & set(parents_2)) != set():
        #print("Parent check fail")
        return None

    new_parents = [-1]*len(undirected_join.partition)
    for i in range(len(undirected_join.partition)):
        pset = undirected_join.partition[i].copy()
        
        #print(p1.partition, ",", parents_1, " <> ", p2.partition, ",", parents_2)
        #print(set(pset), "& (", set(parents_1), "|", set(parents_2), ")")
        
        results = list((set(pset) & (set(parents_1) | set(parents_2))))
        if len(results)==1:
            #print("    ->", results[0])
            new_parents[i] = results[0]
        elif None in parents_1 or None in parents_2 and len(results) == 0:
            #print("    ->", None)
            new_parents[i] = None
        else:
            print("Error")
            

    #if None in parents_1 or None in parents_2:
    #    print(parents_1, "*", parents_2, "=", new_parents)
    #print(undirected_join.partition, new_parents, ":")
    return undirected_join, new_parents

import numpy as np
import dpartitions2 as dp2
def dconnected(p1: prt.Partition, p2: prt.Partition, parents_1, parents_2, bits=False):
    """
    Check whether a pair of partitions is a connected tree (directed portion)
    """
    if not bits:
        all_items_1 = {edge for pset in p1.partition for edge in pset}
        if None in parents_1: 
            all_items_1.add(None)
        all_items_2 = {edge for pset in p2.partition for edge in pset}
        if None in parents_2: 
            all_items_2.add(None)
        children_1 = all_items_1 - set(parents_1)
        children_2 = all_items_2 - set(parents_2)
        in_items = all_items_1 & all_items_2
        out_items = (all_items_1 | all_items_2) - in_items

        if (children_1 & children_2) != set():  # Check if set not empty
            #print("Children check fail")
            return False
        elif (set(parents_1) & set(parents_2)) != set():
            #print("Parent check fail")
            return False
        return True
    else:
        # Do bitwise version here
        (bdict1, bits1) = parents_1
        (bdict2, bits2) = parents_2
        d_obj = dp2.generate_intersection_map([], bdict1, 
                             [], bdict2)
        (D1ord, D2ord) = d_obj["ords"]
        (border_map_1, border_edges_1, border_map_2, border_edges_2) = d_obj["border_edges"]
        num_in = d_obj["num_edges"]
        #uj = ujoin(p1, p2, d_obj)
        #if ujoin is None: continue
        #if uj is not None:
        bitmap1 = np.array(bits1)
        bitmap2 = np.array(bits2)
        R = bitmap1[D1ord].T @ bitmap2[D2ord]
        #print(";;;")
        #print(bitmap1[D1ord].T, bitmap2[D2ord])
        #print(R)
        if R == -num_in:
            return True
        else:
            return False
        #(p1chosen, p2chosen) = np.where(R == -num_in)


def sort_partition(p: prt.Partition, parents, np_array=False, root_ind=None):
    """ 
    Sort the partitions so that you end up with a unique encoding. 
    Input:  p: undirected partition
            parents: directed portion
            np_array: Whether you're using bits
    """
    # Sort general list
    if np_array:
        if parents is not None:
            part = [list(pset) for pset in p.partition]
            for i in range(len(part)):
                part[i].sort()
            if root_ind is not None:
                if root_ind != -1:
                    #root_ind = np.argsort([p[0] for p in part])[root_ind]
                    root_ind = np.where(np.argsort([p[0] for p in part]) == root_ind)[0][0]
                    # print('Hello')
                part.sort()
                (pdict, pbits) = parents
                pars = pbits
                return part, pars, root_ind
            else:
                part.sort()
                (pdict, pbits) = parents
                pars = pbits
                return part, pars
        else:
            part = [list(pset) for pset in p.partition]
            for i in range(len(part)):
                part[i].sort()
            part.sort()
            return part, None
    else:
        if parents is not None:
            pl = [list(cand) for cand in p.partition]
            # Sort individual
            
            for i in range(len(pl)):
                pl[i] = list(pl[i])
                pl[i].sort()

            zipped = list(zip(pl, parents))
            zipped.sort()
            part, pars = zip(*zipped)
            pl = list(part)
            pars = list(pars)

            #print(pl, "; ", pars)
            #print("")
            
            return pl, pars
        else:
            #print(p)
            part = [list(pset) for pset in p.partition]
            for i in range(len(part)):
                part[i].sort()
            part.sort()
            return part, None


def get_encoding(p: prt.Partition, parents=None, sep=False, np_array=False, und=True):
    """ 
    Generate a string encoding of partitions
    Input:  p: The undirected partition
            parents: the directions
            sep: Whether you want the undirected and directed portion separately encoded
            np_array: Whether you're using bits or parent lists
            und: Whether you want to return the undirected portion at all
    """
    if not np_array:
        p2, parents2 = sort_partition(p, parents, False)
        if not sep and parents is not None:
            return str([p2, parents2])
        elif sep and parents is not None:
            return str(p2), str(parents2)
        else:
            return str(p2)
    else:
        if und:
            p2, parents2 = sort_partition(p, parents, True)
            if not sep and parents is not None:
                return str([p2, parents2.tolist()])
            elif sep and parents is not None:
                return str(p2), str(parents2.tolist())
            else:
                return str(p2)
        else:
            #p2, parents2 = sort_partition(p, parents, True)  # REM
            (pdict, parents2) = parents
            return str(parents2.tolist())



def generate_directed(partition_candidates, with_root=False, only_full=True):
    """ 
    Given a list of undirected partition candidates, generate directed candidates
    Input: List[partition]
    Output: Dictionary[partition encoding -> List[parent list]]
     """
    parent_lists = dict()
    for partition in partition_candidates:
        part = [list(pset) for pset in partition.partition]
        
        if not with_root:
            all_parents = product(*part)
            all_parents = [list(p) for p in all_parents]
        elif with_root and not only_full:
            all_parents = []
            for i in range(len(part)):
                partnull = part.copy()
                partnull[i] = [None]
                all_parents += [list(p) for p in product(*partnull)]
        elif with_root and len(part) == 1:
            all_parents = [[None]]
        else:
            all_parents = []

        parent_lists[get_encoding(partition)] = all_parents
    
    #print(parent_lists)
    return parent_lists


def generate_costs(partition_candidates, parent_lists, np_array=False):
    """ Generate random costs. Each node will have one of these dictionaries """
    costs = dict()
    for partition in partition_candidates:
        for plist in parent_lists[get_encoding(partition)]:
            if not np_array:
                costs[get_encoding(partition, plist)] = random.random()
            else:
                costs[get_encoding(partition, plist)] = random.random()
    
    return costs
