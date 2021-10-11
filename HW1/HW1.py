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

    #print(reads)

    for i in range(len(reads)):
        for j in [x for x in range(len(reads)) if x != i]:
            overlap = help.overlap_length(reads[i],reads[j])
            if overlap >= min_overlap:
                edges.append((i,j,overlap))

    edges.sort(key=lambda y: (-y[2],reads[y[0]], reads[y[1]]))

    #print(edges)

    for i in range(len(edges)):
        edge = edges[i]

        if i == 0:
            ham_path.append(edge)
        else:
            if sum([edge[0] == i[0] for i in ham_path]) == 0 and sum([edge[1] == i[1] for i in ham_path]) == 0 and not help.adds_cycle([i for i in ham_path],edge):
                ham_path.append(edge)




    #print(ham_path)

    if len(ham_path) > 0:
        ham_path = help.get_all_ordered_ham_path(ham_path)

    #print(ham_path)

    for sub_path in ham_path:
        contig = []

        for i in range(len(sub_path)):
            if i < len(sub_path) - 1:
                contig.append(reads[sub_path[i][0]][0:len(reads[sub_path[i][0]])-sub_path[i][2]])
            else:
                contig.append(reads[sub_path[i][0]][0:len(reads[sub_path[i][0]])-sub_path[i][2]])
                contig.append(reads[sub_path[i][1]])

        super_str.append(''.join(contig))

    #print(ham_path)

    nodes_in_ham_path = [edge[0] for sub_path in ham_path for edge in sub_path]
    [nodes_in_ham_path.append(edge[1]) for sub_path in ham_path for edge in sub_path]

    #print(nodes_in_ham_path)

    for i in range(len(reads)):
        if i not in nodes_in_ham_path:
            super_str.append(reads[i])

    if min_overlap == 0:
        super_str = [''.join(super_str)]

    super_str.sort()

    #print(super_str)

    return super_str

def read_strings_from_file(filename):
    return [line.rstrip() for line in open(filename)]

def test_greedy_assemble_with_files(reads_filename, superstring_list_filename):
    reads = read_strings_from_file(reads_filename)
    superstring_list = read_strings_from_file(superstring_list_filename)
    assert greedy_assemble(reads) == superstring_list

# TEST: greedy_assemble returns a list of strings
sanity_test_reads = read_strings_from_file("HW1/tests/test_reads.txt")
sanity_test_assembly = greedy_assemble(sanity_test_reads)
assert isinstance(sanity_test_assembly, list), "Return value of greedy_assemble is not a list"
assert all(isinstance(s, str) for s in sanity_test_assembly), "Elements of list are not all strings"
print("SUCCESS: greedy_assemble returns a list of strings passed!")

# TEST: greedy_assemble returns a superstring list
def check_is_superstring_list(assembly, reads):
    for read in reads:
        assert any(read in contig for contig in assembly), f"read '{read}' is not contained in assembly"

sanity_test_assembly = greedy_assemble(sanity_test_reads)
check_is_superstring_list(sanity_test_assembly, sanity_test_reads)
print("SUCCESS: greedy_assemble returns a superstring list passed!")

# TEST: sanity_test assembly min_overlap=0
sanity_test_assembly = greedy_assemble(sanity_test_reads, 0)
assert sanity_test_assembly == ['the_quick_brown_fox_jumps_over_the_lazy_dog']
print("SUCCESS: sanity_test_assembly_min_overlap_0 passed!")

# TEST: sanity_test assembly min_overlap=2
sanity_test_assembly = sorted(greedy_assemble(sanity_test_reads, 2))
assert sanity_test_assembly == ['ck_brown_fox_ju', 'er_the_lazy_dog', 'the_quic', 'umps_ove']
print("SUCCESS: sanity_test_assembly_min_overlap_2 passed!")

# TEST: greedy_assemble_small_test_1
small_test1_reads = ["GTT", "ATCTC", "CTCAA"]
assert greedy_assemble(small_test1_reads) == ["ATCTCAAGTT"]
print("SUCCESS: greedy_assemble_small_test_1 passed!")

# TEST: greedy_assemble_small_test_2_0
small_test2_reads = ["CGAAG", "ATCGA", "AGAG", "GGG"]
assert greedy_assemble(small_test2_reads, 0) == ["ATCGAAGAGGG"]
print("SUCCESS: greedy_assemble_small_test_2_0 passed!")

# TEST: greedy_assemble_small_test_2_1
small_test2_reads = ["CGAAG", "ATCGA", "AGAG", "GGG"]
assert greedy_assemble(small_test2_reads, 1) == ["ATCGAAGAGGG"]
print("SUCCESS: greedy_assemble_small_test_2_1 passed!")

# TEST: greedy_assemble_small_test_2_2
small_test2_reads = ["CGAAG", "ATCGA", "AGAG", "GGG"]
assert sorted(greedy_assemble(small_test2_reads, 2)) == ['ATCGAAGAG', 'GGG']
print("SUCCESS: greedy_assemble_small_test_2_2 passed!")

# TEST: greedy_assemble_small_test_2_3
small_test2_reads = ["CGAAG", "ATCGA", "AGAG", "GGG"]
assert sorted(greedy_assemble(small_test2_reads, 3)) == ['AGAG', 'ATCGAAG', 'GGG']
print("SUCCESS: greedy_assemble_small_test_2_3 passed!")

# TEST: greedy_assemble_small_test_2_4
small_test2_reads = ["CGAAG", "ATCGA", "AGAG", "GGG"]
assert sorted(greedy_assemble(small_test2_reads, 4)) == ['AGAG', 'ATCGA', 'CGAAG', 'GGG']
print("SUCCESS: greedy_assemble_small_test_2_4 passed!")

# TEST: greedy_assemble_small_test_2_4
small_test2_reads = ["CGAAG", "ATCGA", "AGAG", "GGG"]
assert sorted(greedy_assemble(small_test2_reads, 4)) == ['AGAG', 'ATCGA', 'CGAAG', 'GGG']
print("SUCCESS: greedy_assemble_small_test_2_4 passed!")

# TEST: greedy_assemble_small_test_3
small_test3_reads = ["C", "A", "T", "G"]
assert greedy_assemble(small_test3_reads) == ["ACGT"]
print("SUCCESS: greedy_assemble_small_test_3 passed!")

# TEST: greedy_assemble large test 1
test_greedy_assemble_with_files("HW1/tests/large_test1_reads.txt", "HW1/tests/large_test1_superstring.txt")
print("SUCCESS: greedy_assemble large test 1 passed!")

#print(str(greedy_assemble(["C","A","T","G"])) == str(['ACGT']))

#print(str(greedy_assemble(["CGAAG", "ATCGA", "AGAG", "GGG"],4)) == str(['AGAG', 'ATCGA', 'CGAAG', 'GGG']))

#print(str(greedy_assemble(["GGG","CGAAG", "ATCGA", "AGAG"],3)) == str(['AGAG', 'ATCGAAG', 'GGG']))

#print(str(greedy_assemble(["GGG","CGAAG", "ATCGA", "AGAG"],2)) == str(['ATCGAAGAG','GGG']))

#print(str(greedy_assemble(["GGG","CGAAG", "ATCGA", "AGAG"],1)) == str(['ATCGAAGAGGG']))

#(str(greedy_assemble(["GGG","CGAAG", "ATCGA", "AGAG", "ATC"],3)) == str(['AGAG', 'ATCGAAG', 'GGG']))

#sanity_test_reads = help.read_strings_from_file("HW1/tests/test_reads.txt")

#print(str(greedy_assemble(sanity_test_reads)) == str(['the_quick_brown_fox_jumps_over_the_lazy_dog']))

#print(greedy_assemble(['the_qu','mps_ov','e_qui','ps_ove','ck_bro']))

#print(greedy_assemble(['the_qu','e_qui','ck_bro','wn_f','n_fox_','ox_ju','mps_ov']))

#print(help.order_ham_path(ham_path=[(4,7,3),(0,2,3),(2,4,3),(7,8,3)]))

#print(help.adds_cycle([(0, 1, 4), (2, 6, 0), (1, 2, 0)], (6, 0, 0)))
