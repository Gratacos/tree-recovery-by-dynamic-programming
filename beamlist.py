import numpy as np

#from merge_seq import priority_cost

class BeamList:
    def __init__(self, params):
        self.limit = params.init_width
        self.type = params.beam_type
        self.hist = params.histogram
        self.locations = dict()  # partition_encoding -> location in list
        self.beam_nodes = []
        self.heap = params.beam_heap
        self.min_width = params.min_width
        self.verbose = params.verbose
        self.minnorm = params.minnorm
        self.norm = params.norm

    def heapify(self, ind=None):
        # This is a heap-insertion, not an actual "heapify"
        key = self.beam_nodes[ind][2]
        if ind is not None and key is not None:
            # Note, this only goes up, assumming costs only go down
            self.locations[key] = ind  # Initialize, in case loop is not executed
            while self.beam_nodes[ind][0] > self.beam_nodes[int((ind - 1)/2)][0]:
                par_ind = int((ind - 1)/2)
                (par_cost, par_val, par_key) = self.beam_nodes[par_ind]

                self.beam_nodes[ind], self.beam_nodes[par_ind] = self.beam_nodes[par_ind], self.beam_nodes[ind]
                self.locations[par_key], self.locations[key] = ind, par_ind

                ind = par_ind

    def add(self, cost, value, key, prioritize=False):
        """ 
        cost:       Value to sort by in the beam-search
        value:      Partition
        key:        Unique encoding to identify the value
        prioritize: Whether this value should always be at the start of the list
        """ 
        if key in self.locations:
            # If value exists, update cost, will still need to heapify afterwards
            ind = self.locations[key]
            self.beam_nodes[ind] = (cost, self.beam_nodes[ind][1], self.beam_nodes[ind][2])
        else:
            self.beam_nodes.append( (cost, value, key) )
            ind = len(self.beam_nodes)-1
            self.locations[key] = ind

        if self.limit < float("inf") and self.heap:
            self.heapify(ind=ind)
            for i in reversed(range(self.limit, len(self.beam_nodes))):
                del self.locations[self.beam_nodes[i][2]]
                del self.beam_nodes[i]

    def get_list(self, priority_key=None):
        #print(self.beam_nodes)
        if not self.heap:
            if self.type == 'fixed' and self.limit < len(self.beam_nodes):
                # Fixed beam: Same width throughout entire erosion process
                if priority_key is not None:
                    #print(self.locations.keys(), priority_key )
                    priority_index = self.locations[priority_key]
                    priority_value = self.beam_nodes[priority_index]
                    # del self.beam_nodes[priority_index]
                    self.beam_nodes.sort(key=lambda a: a[0])
                    # self.beam_nodes.insert(0, priority_value)
                    self.beam_nodes = self.beam_nodes[0:self.limit]
                    if priority_key not in [k for (n, v, k) in self.beam_nodes]:
                        #self.beam_nodes.pop()
                        #self.beam_nodes.insert(0, priority_value)
                        self.beam_nodes[-1] = priority_value
                else:
                    self.beam_nodes.sort(key=lambda a: a[0])
                    self.beam_nodes = self.beam_nodes[0:self.limit]
            elif self.type == 'percent' and self.limit <= 1.0 and self.limit < max(self.min_width, int(self.limit*len(self.beam_nodes))):
                # Percent beam: Cut the same percent of nodes throughout entire erosion (with minimum width)
                self.beam_nodes.sort(key=lambda a: a[0])
                cutoff = max(self.min_width, int(self.limit*len(self.beam_nodes)))
                if self.verbose:
                    print(len(self.beam_nodes), "cut to", cutoff, "(" + str(int(self.limit*100)) + "%)")

                self.beam_nodes = self.beam_nodes[0:cutoff]
            elif self.type == 'probability':
                # Cumulative likelihood cutoff
                self.beam_nodes.sort(key=lambda a: a[0])
                summ = 1
                it = len(self.beam_nodes)
                for i in range(len(self.beam_nodes)):
                    summ += self.beam_nodes[i][0]
                    if 1 - 1/summ > self.limit:
                        it = max(self.min_width, i)
                        break
                if self.verbose:
                    print(len(self.beam_nodes), "cut to", it, "(" + str(1 - 1/summ) + " cost reached)")

                self.beam_nodes = self.beam_nodes[0:it]
            elif (self.type == 'delta min' and self.min_width < len(self.beam_nodes)) or self.type == 'delta':
                # Fixed percent of likelihood
                costs = [bm[0] for bm in self.beam_nodes]
                maxcost, mincost = max(costs), min(costs)
                epsilon = 0.00001

                if self.norm == 'division':
                    costs = [((costs[i])/(maxcost + epsilon)) for i in range(len(costs))]
                elif self.norm == 'minmax':
                    costs = [((costs[i] - mincost)/(maxcost - mincost + epsilon)) for i in range(len(costs))]
                
                bestcost = min(costs)
                self.beam_nodes = [self.beam_nodes[i] for i in range(len(self.beam_nodes)) if (costs[i] - bestcost) <= self.limit]

            # elif self.type == 'delta min' and self.min_width < len(self.beam_nodes):
            #     # Fixed percent of likelihood, always have at least 500 elems
            #     self.beam_nodes.sort(key=lambda a: a[0])
            #     costs = [bm[0] for bm in self.beam_nodes]
            #     maxcost, mincost = max(costs), min(costs)
            #     epsilon = 0.00001

            #     if self.norm == 'division':
            #         costs = [((costs[i])/(maxcost + epsilon)) for i in range(len(costs))]
            #     elif self.norm == 'minmax':
            #         costs = [((costs[i] - mincost)/(maxcost - mincost + epsilon)) for i in range(len(costs))]

            #     bestcost = min(costs)
            #     cutoff = len(self.beam_nodes)
            #     for i in range(self.min_width+1, len(costs)):
            #         if (costs[i] - bestcost) > self.limit:
            #             cutoff = i
            #             break
                
            #     self.beam_nodes = self.beam_nodes[0:cutoff]
                
        if self.hist:
            self.per_join_cost = {
                self.beam_nodes[i][2] : self.beam_nodes[i][0] for i in range(len(self.beam_nodes))
            }
        return [self.beam_nodes[i][1] for i in range(len(self.beam_nodes))]