import RNA
import math


seq   = "ATAATGATACTTCCGTCCAGTCACAGTCCGAGCGTGAAGCGGCAAGCCTGCCGCAATCCGCCGCAGGCAATGAGGGCGGCCCGGTGAGATCCACGGGGTGGTGAAATTCCAATGA"
kT    = 0.61632077549999997
pairs = [[ 0. for j in range(0, len(seq) + 1) ] for j in range(0, len(seq) + 1)]
pf    = 0.

fc = RNA.fold_compound(seq)
# add soft constraint for base pair (1,10)
# fc.sc_add_bp(1,10,-3)

for s in fc.subopt(150):
    print(s.structure)
    # # add contribution of structure to partition function
    # pf = pf + math.exp(-s.energy / kT)
    # # extract base pairs from structure
    # pairtable = RNA.ptable(s.structure)
    # for i in range(1, len(seq) + 1):
    #     if pairtable[i] > i:
    #         pairs[i][pairtable[i]] = pairs[i][pairtable[i]] + math.exp(-s.energy / kT)

# print ("Free energy of ensemble: %g\n" % ( -kT * math.log(pf)))

# for i in range(1, len(seq) + 1):
#     for j in range(i + 1, len(seq) + 1):
#         pairs[i][j] = pairs[i][j] / pf
#         if pairs[i][j] > 0:
#             print ("pr(%d, %d) = %g" % (i, j, pairs[i][j]))

# sequence = "GGGGAAAACCCC"
# # Set global switch for unique ML decomposition
# RNA.cvar.uniq_ML = 1

# subopt_data = { 'counter' : 1, 'sequence' : sequence }

# # Print a subopt result as FASTA record
# def print_subopt_result(structure, energy, data):
#     if not structure == None:
#         print (">subopt %d" % data['counter'])
#         print ("%s" % data['sequence'])
#         print ("%s [%6.2f]" % (structure, energy))
#         # increase structure counter
#         data['counter'] = data['counter'] + 1

# # Create a 'fold_compound' for our sequence
# a = RNA.fold_compound(sequence)

# # Enumerate all structures 500 dacal/mol = 5 kcal/mol arround
# # the MFE and print each structure using the function above
# a.subopt_cb(500, print_subopt_result, subopt_data)
