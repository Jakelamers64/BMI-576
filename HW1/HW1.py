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

    for i in range(len(reads)):
        for j in [x for x in range(len(reads)) if x != i]:
            overlap = help.overlap_length(reads[i],reads[j])
            if overlap >= min_overlap:
                edges.append((i,j,overlap))

    edges.sort(key=lambda y: (-y[2],reads[y[0]] + reads[y[1]]))

    return super_str

print(str(greedy_assemble(["C","A","T","G"])) == str(['ACGT']))
