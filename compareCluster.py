import multiprocessing
from functools import partial
import pickle
from singleRunCompareCluster import *


def compare_clusters(all_family_list, num_of_structures, num_of_sequences):
    pool = multiprocessing.Pool()
    func = partial(run_all_clusters, all_family_list)
    # results = pool.map(func, [(x, y, z) for x in range(1, num_of_structures) for y in range(1, num_of_sequences) for z in range(4, len(all_family_list)+1)])
    results = pool.map(func, [(x, y) for x in range(3, num_of_structures+1) for y in range(3, num_of_sequences+1)])

    pool.close()
    pool.join()
    
    with open('result_50familiy.data', 'wb') as filehandle:
        # store the data as binary data stream
        pickle.dump(results, filehandle)


def run_all_clusters(all_family_list, zip_structure_sequence_family):
    # num_of_structures, num_of_sequences, num_of_family= zip_structure_sequence_family
    num_of_structures, num_of_sequences = zip_structure_sequence_family


    # family_list = all_family_list[:num_of_family]
    family_list = all_family_list
    return singleRunCompareCluster(family_list, num_of_structures, num_of_sequences)
    


