import networkx as nx
##############################################################################
#                            CONSENSUS FUNCTIONS                             #
##############################################################################

KM_SIZE = 16
THRESHOLD = 0.75


def find_consensus(g, nb_run, threshold=0.75):
    path = []
    best_node = max(g.nodes, key=lambda x: g.nodes[x]['weight'])
    path.append(best_node)
    g.nodes[best_node]['path'] = 'consensus'
    left = best_node
    right = best_node
    search = True
    while(search):
        if(left):
            l_list = list(g.predecessors(left))
            left = l_list[0] if l_list else None
        if(right):
            r_list = list(g.successors(right))
            right = r_list[0] if r_list else None
        if(left and g.nodes[left]['weight'] / nb_run >= threshold):
            path = [left] + path
            g.nodes[left]['path'] = 'consensus'
        # Adding
        if(right and g.nodes[right]['weight'] / nb_run >= threshold):
            path.append(right)
            g.nodes[right]['path'] = 'consensus'
        # If no k-mer are added to the path, stop.
        if(not (path[0] == left or path[-1] == right)):
            search = False
    return(path)


def build_consensus_graph(adp_dict):
    # print("Working on ",method, which_end, file=sys.stderr)
    g = nx.DiGraph()
    for adp, w in adp_dict.items():
        prev = ""
        for i in range(len(adp) - KM_SIZE + 1):
            km = adp[i:i + KM_SIZE]
            if(km not in g.nodes):
                g.add_node(km, weight=w)
            else:
                g.nodes[km]["weight"] += w
            if(prev != ""):
                g.add_edge(prev, km)
            prev = km
    return(g)