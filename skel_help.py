import networkx as nx 
import numpy as np
import math
from plyfile import PlyData, PlyElement
import os.path 
import h5py
import scipy as scp 
from scipy import ndimage
import partitions
import dpartitions


def load_costs(file_name, G):
    # This loads both the costs and the partitions
    costs = {n:{} for n in G}
    dmap = {0:-1.0, 1:1.0}
    ign = {' ', '\n'}

    with open(file_name, 'r') as f:
        next_node = 1
        for line in f.readlines():
            pcands = []
            dir_cands = {}
            parent_dict = {}
            for candidate in [c for c in line.split(';') if c not in ign]:
                partition = [cand for cand in candidate.split(',') if cand not in ign]
                
                edge_names = [int(n) for n in partition[0].split(' ')]
                undirected = [int(n) for n in partition[1].split(' ')]
                directed   = [int(n) for n in partition[2].split(' ')]
                cost       = float(partition[-1])

                if len(parent_dict) == 0:
                    parent_dict = {edge_names[i]:i for i in range(len(edge_names))}

                part = [[] for i in range(max(undirected)+1)] # []*(max(undirected)+1)
                dirs = [0]*len(directed)

                for i in range(len(edge_names)):
                    part[undirected[i]].append(edge_names[i])
                    dirs[parent_dict[edge_names[i]]] = dmap[directed[i]]
                
                part = partitions.Partition(part)
                pcands.append(part)
                part_key = dpartitions.get_encoding(part)
                cand_key = dpartitions.get_encoding(part, parents=(parent_dict, np.array(dirs)), sep=False, np_array=True)

                if part_key not in dir_cands:
                    dir_cands[part_key] = [dirs]
                else:
                    dir_cands[part_key].append(dirs)

                if all([d == -1 for d in dirs]):
                    G.graph['root'] = next_node
                costs[next_node][cand_key] = cost

            G.nodes[next_node]['partition_candidates'] = pcands
            G.nodes[next_node]['parent_dict'] = parent_dict
            G.nodes[next_node]['parent_lists'] = {k:np.array(v).T for (k, v) in dir_cands.items()}

            next_node += 1
    return costs 
        



def ply_to_graph(file_name, normalize_width_radius=False):
    """ Load a NetworkX graph from a PLY file. """
    print("Loading graph from", file_name, "...")

    plydata = PlyData.read(file_name)

    elements_with_type = plydata.elements
    for element in elements_with_type:
        el_name = element.name 
        
        if el_name == 'vertex':
            num_verts = element.count
            keys = list(range(num_verts))

            v_property_list = dict()
            for prop in element.properties:
                prop_value = plydata[el_name][prop.name]

                # We store width and radius later if we want to normalize them. In a couple of lines,
                # you'll see that we need both the radius and the width to normalize. That way we
                # preserve the ratio. 
                if normalize_width_radius and prop.name == "radius":
                    tmp_rad = prop_value
                elif normalize_width_radius and prop.name == "bt2":
                    tmp_wid = prop_value
                else:
                    v_property_list[prop.name] = dict(zip(keys,list(prop_value))) 

            if normalize_width_radius:
                rad = tmp_rad
                wid = tmp_wid
                max_val = max(rad.max(), wid.max())
                min_val = min(rad.min(), wid.min())
                rad = (rad-min_val)/(max_val-min_val)
                wid = (wid-min_val)/(max_val-min_val)
                v_property_list["radius"] = dict(zip(keys, list(rad)))
                v_property_list["bt2"] = dict(zip(keys, list(wid)))

        elif el_name == 'edge':
            v1 = plydata[el_name]['vertex1']
            v2 = plydata[el_name]['vertex2']
            edges = list(zip(list(v1), list(v2)))

            e_property_list = dict()
            for prop in element.properties:
                if prop.name != 'vertex1' and prop.name != 'vertex2':
                    prop_value = plydata[el_name][prop.name]
                    e_property_list[prop.name] = dict(zip(edges,list(prop_value)))
    
    g = nx.Graph(edges)
    for name in v_property_list:
        nx.set_node_attributes(g,v_property_list[name],name)
    for name in e_property_list:
        nx.set_edge_attributes(g,e_property_list[name],name)

    print("Graph with", g.number_of_nodes(), "nodes and", g.number_of_edges(), " edges loaded.")
    
    return g


# TODO: Make sure every array has ordering from 0 to num_vertex
def graph_to_ply(graph, file_name):
    """ Saves a NetworkX graph to a ply file """
    g = nx.convert_node_labels_to_integers(graph)

    print("Saving ply file...")

    # Extract the properties of the nodes
    sample = list(g.nodes())[0]
    vertex_properties = list(g.nodes[sample].keys())

    # Extract the values of said properties:
    vertex_prop_values = []
    dt = []
    for prop in vertex_properties:
        # Obtain the type of property, we want to know, because if
        # the property is a list, we are not able to store it

        # Not ordered:
        v_props = []
        v_order = []
        for v in g.nodes():
            v_order.append(v)
            v_props.append(g.nodes[v][prop])
        to_sort = np.array([v_order, v_props]).transpose()
        to_sort = to_sort[np.argsort(to_sort[:, 0])]
        v_props = list(to_sort[:, 1])
        # v_props = list(nx.get_node_attributes(g, name=prop).values())

        if np.array(v_props[0]).dtype == np.dtype('i8'):
            dt_tmp = np.dtype('i4')
        else:
            dt_tmp = np.array(v_props[0]).dtype

        if type(v_props[0]) != list:
            dtype = (prop, dt_tmp)
        else:
            dtype = (prop, dt_tmp,(len(v_props),))

        vertex_prop_values.append(v_props)

        dt.append(dtype)
    vertex_prop_values = list(zip(*vertex_prop_values))
    vertex_prop_values = np.array(vertex_prop_values, dtype=dt)
    
    # Extract properties for the edges
    if type(g) == nx.MultiGraph:
        sample = list(g.edges(keys=True))[0]
    else:
        sample = list(g.edges())[0]
    edge_properties = list(g.edges[sample].keys())

    # Extract the values of said properties:
    edges = list(g.edges())
    edges = zip(*edges) # I want it [[source1,source2,...],[target1,target2,...]]
    edge_prop_values = list(map(list,edges)) # Remove tuples
    dt = [('vertex1', 'i4'),('vertex2', 'i4')]
    for prop in edge_properties:
        if prop == "extremes":
            continue # TODO: this is bandaid code. Better way is to just not store lists in general
        e_props = list(nx.get_edge_attributes(g, name=prop).values())

        if np.array(e_props[0]).dtype == np.dtype('i8'):
            dt_tmp = np.dtype('i4')
        else:
            dt_tmp = np.array(e_props[0]).dtype

        if type(e_props[0]) != list:
            dtype = (prop, dt_tmp)
        else:
            dtype = (prop, dt_tmp, (len(e_props),))

        edge_prop_values.append(e_props)
        dt.append(dtype)
    edge_prop_values = list(zip(*edge_prop_values))
    edge_prop_values = np.array(edge_prop_values, dtype=dt)
    
    # Now for PlyFile stuff:
    v = PlyElement.describe(vertex_prop_values, 'vertex')
    e = PlyElement.describe(edge_prop_values, 'edge')
    PlyData([v,e], text=True).write(file_name)

    print("Graph saved to", file_name)


def graph_to_csv(g, file_name_vertex, file_name_edges):
    """ Saves a NetworkX graph to a CSV with the vertex properties, and a csv with the edges. """
    # Extract the properties of the nodes:
    if g.number_of_nodes() > 0:
        sample = list(g.nodes())[0]
        vertex_properties = list(g.nodes[sample].keys())

        # Extract the values of said properties:
        vertex_prop_values = []
        for prop in vertex_properties:
            # Not ordered:
            v_props = []
            v_order = []
            for v in g.nodes():
                v_order.append(v)
                v_props.append(g.nodes[v][prop])
            to_sort = np.array([v_order,v_props]).transpose()
            to_sort = to_sort[np.argsort(to_sort[:, 0])]
            v_props = list(to_sort[:,1])

            vertex_prop_values.append(v_props)

        vertex_prop_values = np.array(vertex_prop_values).transpose()
        vertex_prop_values = np.append(np.array([vertex_properties]), vertex_prop_values, axis=0)
    
    if g.number_of_edges() > 0:
        if type(g) == nx.MultiGraph or type(g) == nx.classes.graphviews.SubMultiGraph:
            sample = list(g.edges(keys=True))[0]
            edge_properties = ["vertex1", "vertex2", "key"] + list(g.edges[sample].keys())
            edges = list(g.edges(keys=True))
            prop_start=3
        else: 
            sample = list(g.edges())[0]
            edge_properties = ["vertex1", "vertex2"] + list(g.edges[sample].keys())
            edges = list(g.edges())
            prop_start=2

        # Extract the values of said properties:
        edges = zip(*edges) # I want it [[source1,source2,...],[target1,target2,...]]
        edge_prop_values = list(map(list,edges)) # Remove tuples
        for prop in edge_properties[prop_start:]:
            e_props = list(nx.get_edge_attributes(g, name=prop).values())

            edge_prop_values.append(e_props)

        edge_prop_values = np.array(edge_prop_values).transpose()
        edge_prop_values = np.append(np.array([edge_properties]), edge_prop_values, axis=0)
    
    # Now for PlyFile stuff:
    
    np.savetxt(file_name_vertex, vertex_prop_values, delimiter=",", fmt='%s')
    np.savetxt(file_name_edges, edge_prop_values, delimiter=",", fmt='%s')


def graph_to_cmp(g, file_name):
    """ Saves a NetworkX graph into a cmp file for chimera viewing. """
    x = nx.get_node_attributes(g,'x').values()
    y = nx.get_node_attributes(g,'y').values()
    z = nx.get_node_attributes(g,'z').values()
    mx = max(x)
    my = max(y)
    mz = max(z)
    vol = np.zeros([mx+1, my+1, mz+1])
    vol[x,y,z] = 1

    # Save:
    save_vol = np.ascontiguousarray(vol, dtype=np.uint8)
    out = h5py.File(file_name, 'w')
    out.create_group('map0').create_dataset('data0',data=save_vol,chunks=True,compression='gzip')
    out.close()


def node_filter(graph, nodes, attribute_filter=lambda g, node: True):
    """ Filter out nodes """
    nodes =[n for n in nodes if attribute_filter(graph, n)]

    return nx.subgraph(graph, nodes).copy()


def edge_filter_out(graph, attribute_filter=lambda g, edge: False, remove_nodes=False):
    """ Filter out edges. Set remove_nodes to true to remove the associated nodes """
    if type(graph) == nx.MultiGraph or type(graph) == nx.classes.graphviews.SubMultiGraph:
        edges = graph.edges(keys=True)
    else:
        edges = graph.edges()
    if remove_nodes: 
        nodes_to_remove = set()
        nodes_to_include = set()
        for edge in edges:
            if attribute_filter(graph, edge):
                nodes_to_remove.add(edge[0])
                nodes_to_remove.add(edge[1])
            else:
                nodes_to_include.add(edge[0])
                nodes_to_include.add(edge[1])
        nodes_to_include = nodes_to_include.difference(nodes_to_remove)
        return nx.subgraph(graph, list(nodes_to_include) ).copy()
    else:
        g_copy = graph.copy()
        edges = [edge for edge in edges if attribute_filter(graph, edge)]
        g_copy.remove_edges_from(edges)
        return g_copy


def skel_to_vol(skeleton, attr="radius", res=0.2, node_at=False):
    """ Convert a skeleton to a volume. """
    x = nx.get_node_attributes(skeleton, "x")
    y = nx.get_node_attributes(skeleton, "y")
    z = nx.get_node_attributes(skeleton, "z")

    if not node_at:
        at = nx.get_edge_attributes(skeleton, attr)
    else:
        at = nx.get_node_attributes(skeleton, attr)

    # Find the size of the volume from the coordinates
    xmax = max(x.values())
    ymax = max(y.values())
    zmax = max(z.values())

    vol = np.zeros((int(xmax)+1,int(ymax)+1,int(zmax)+1))

    if type(skeleton) == nx.Graph:
        edges = skeleton.edges()
    elif type(skeleton) == nx.MultiGraph:
        edges = skeleton.edges(keys=True)

    for edge in skeleton.edges():
        u = edge[0]
        v = edge[1]

        if not node_at:
            w = int(at[edge]*255)
        else:
            w = int(at[u]*255)

        # Include edge in volume:
        u = np.array([x[u], y[u], z[u]])
        v = np.array([x[v], y[v], z[v]])

        vol[int(u[0]), int(u[1]), int(u[2])] = w
        vol[int(v[0]), int(v[1]), int(v[2])] = w

        v_dist = v - u

        for i in np.arange(0, 1, res):
            ncoord = u + i*v_dist
            vol[int(ncoord[0]), int(ncoord[1]), int(ncoord[2])] = w

    return vol

