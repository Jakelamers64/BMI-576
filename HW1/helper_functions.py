def overlap_length(left, right):
    """Returns the length of the longest suffix of left that is a prefix of right

    Args:
        left: a string
        right: a string
    Returns:
        An integer length of the longest overlap (0 if there is no overlap)
    """
    if left == right:
        return len(left)
    else:
        for i in range(0,len(left))[::-1]:
            if right.startswith(left[(len(left) - i - 1):(len(left))]):
                return len(left[(len(left) - i - 1):(len(left))])

        return 0

def read_strings_from_file(filename):
    return [line.rstrip() for line in open(filename)]

def test_greedy_assemble_with_files(reads_filename, superstring_list_filename):
    reads = read_strings_from_file(reads_filename)
    superstring_list = read_strings_from_file(superstring_list_filename)
    assert greedy_assemble(reads) == superstring_list

def get_all_ordered_ham_path(ham_path):
    ordered_ham_paths = []

    if len(ham_path) > 0:

        edges_out = [i[1] for i in ham_path]
        starts_of_paths = [i for i in ham_path if i[0] not in edges_out]

        if len(starts_of_paths) > 1:
            for edge in starts_of_paths:
                ordered_ham_paths.append(order_ham_path(ham_path,edge))
        else:
            ordered_ham_paths.append(order_ham_path(ham_path))

        return ordered_ham_paths
    else:
        return ordered_ham_paths


def order_ham_path(ham_path,edge=None):
    assert ham_path != None

    ordered = []
    edges_out = [i[1] for i in ham_path]
    edges_in = [i[0] for i in ham_path]

    if edge == None:
        if len([i for i in ham_path if i[0] not in edges_out]) > 0:
            edge = [i for i in ham_path if i[0] not in edges_out][0]
        else:
            # ham_path contains no possible start locations
            return ham_path

        [ordered.append(i) for i in order_ham_path(ham_path,edge)]

        return ordered
    elif edge[1] in edges_in:
        ordered.append(edge)

        edge = [i for i in ham_path if i[0] == edge[1]][0]

        [ordered.append(i) for i in order_ham_path(ham_path,edge)]

        return ordered
    else:
        ordered.append(edge)
        return ordered

# refrenced code https://algocoding.wordpress.com/2015/04/02/detecting-cycles-in-a-directed-graph-with-dfs-python/
def adds_cycle(G,edge):
    # adds next edge to the ham_path
    G.append(edge)

    # calcs the number of nodes with the potential of starting a ham_path
    edges_out = [i[1] for i in G]
    starts_of_paths = [i for i in G if i[0] not in edges_out]

    if len(G) > 0 and not len(starts_of_paths) > 0:
        return True

    if len(starts_of_paths) > 1:
        ordered_G = get_all_ordered_ham_path(G)
    else:
        ordered_G = [G]

    G = ordered_G

    contains_cycle = [False]

    for path in G:

        visited = { node : False for node in path  }

        for node in visited:
            if not visited[node]:
                dfs_visit(G,node,visited,contains_cycle)
            if contains_cycle[0]:
                break

        if contains_cycle[0]:
            break

    return contains_cycle[0]

def dfs_visit(G,node,visited,contains_cycle):
    if contains_cycle[0]:
        return

    visited[node] = True

    next_node = [i for i in G if i[0] == node[1]]

    if len(next_node) > 0:
        next_node = next_node[0]
    else:
        return

    if visited[next_node]:
        contains_cycle[0] = True
        return
    if not visited[next_node]:
        dfs_visit(G,next_node,visited,contains_cycle)

#print(adds_cycle([],(1,0,0)))
