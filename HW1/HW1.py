import helper_functions as help

# Code for PROBLEM 1
# You are welcome to develop your code as a separate Python module
# and import it here if that is more convenient for you.
def greedy_assemble(reads, min_overlap=0):
    """Assembles a set of reads using the graph-based greedy algorithm.

    Args:
        reads: a list of strings
        min_overlap: the minimum length of an allowed overlap between two reads
    Returns:
        A list of strings (contigs) that collectively contain all input reads
    """
    edges = []
    ham_path = []
    super_str = []
    visited = []

    for i in range(len(reads)):
        for j in [x for x in range(len(reads)) if x != i]:
            overlap = help.overlap_length(reads[i],reads[j])
            if overlap >= min_overlap:
                edges.append((i,j,overlap))

    edges.sort(key=lambda y: (-y[2],reads[y[0]] + reads[y[1]]))

    #print(edges)

    for i in range(len(edges)):
        edge = edges[i]

        if i == 0:
            ham_path.append(edge)
            visited.append(edge[0])
        else:
            if sum([edge[0] == i[0] for i in ham_path]) == 0 and sum([edge[1] == i[1] for i in ham_path]) == 0 and sum([edge[1] == i for i in visited]) == 0:
                ham_path.append(edge)
                visited.append(edge)

    #print(ham_path)

    for i in range(len(ham_path)):
        if i < len(ham_path) - 1:
            super_str.append(reads[ham_path[i][0]])
        else:
            super_str.append(reads[ham_path[i][0]])
            super_str.append(reads[ham_path[i][1]])

    if len(super_str)  > 0:
        super_str = [''.join(super_str)]

    nodes_in_ham_path = [edge[0] for edge in ham_path]
    [nodes_in_ham_path.append(edge[1]) for edge in ham_path]

    for i in range(len(reads)):
        if i not in nodes_in_ham_path:
            super_str.append(reads[i])

    super_str.sort()

    #print(super_str)

    return super_str

print(str(greedy_assemble(["C","A","T","G"])) == str(['ACGT']))

print(str(greedy_assemble(["CGAAG", "ATCGA", "AGAG", "GGG"],4)) == str(['AGAG', 'ATCGA', 'CGAAG', 'GGG']))
