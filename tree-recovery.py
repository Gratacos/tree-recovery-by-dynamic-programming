import sys
import networkx
import partitiongraph
import merge_seq
import skel_help
import greedy_solution

if __name__ == "__main__":
    if len(sys.argv) > 3:
        graph_fname = sys.argv[1]
        costs_fname = sys.argv[2]
        beam_width  = int(sys.argv[3])
    elif len(sys.argv) > 2:
        graph_fname = sys.argv[1]
        costs_fname = sys.argv[2]
        beam_width = float('inf')
    else:
        graph_fname = "../../2d_datasets/rice/graph_01.ply"
        costs_fname = "../../2d_datasets/rice/graph_01.txt"
        beam_width = 1#float('inf')

    # Load graph (command line argument 1)
    G = skel_help.ply_to_graph(graph_fname)
    nodes = sorted(G)
    networkx.relabel_nodes(G, {nodes[i]:(i+1) for i in range(len(nodes))}, copy=False)
    G = partitiongraph.PGraph(G)
    

    # Load costs (command line argument 2)
    costs =  skel_help.load_costs(costs_fname, G)
    G.graph['cost'] = costs 

    # Compute merge sequence
    sequence =  merge_seq.compute_merge_sequence(G)

    # Greedy solution
    if beam_width != float('inf'):
        timings =  greedy_solution.get_greedy(G)

        for n in G:
            G.nodes[n]['partition'] = G.nodes[n]['greedy_partition']
            G.nodes[n]['parents'] = G.nodes[n]['greedy_parents']
    else:
        for n in G:
            if 'greedy_partition' in G.nodes[n]:
                del G.nodes[n]['greedy_partition']
            if 'greedy_parents' in G.nodes[n]:
                del G.nodes[n]['greedy_parents']
    
    G.graph['root_edges'] = {G.edges[e]['name'] for e in G.edges(G.graph['root'], keys=True)}
    # for n in G.nodes():
    #     for partition in G.nodes[n]['partition_candidates']:
    #         partition.root = -1
    # for i in range(len(G.nodes[G.graph['root']]['partition_candidates'])):
    #     G.nodes[G.graph['root']]['partition_candidates'][i].root = 0

    # G.nodes[G.graph['root']]['greedy_partition'].root = 0

    # Perform DP algorithm
    G.graph['root_edges'] = [G.edges[e]['name'] for e in G.edges(G.graph['root'], keys=True)]
    [new_cost, backtrack] =  partitiongraph.erosion_merge(G, sequence, costs)
    partitiongraph.calculate_optimal(G, new_cost)
    backtrackeds =  partitiongraph.backtrack_erosion(G, backtrack, directed=True, bits=True)

