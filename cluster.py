from super_n_motif.project import *
from subopt import *
import json


def generate_cluster(family_name, num_of_clusters, num_of_sequences):
    dict_RNA = compute_structure(family_name, num_of_clusters, num_of_sequences)

    # transform dictionary to dataframe
    dataframe_RNA = pd.DataFrame(data=dict_RNA)
    print(dataframe_RNA)

    for n in range(num_of_clusters):
        id = ">subopt %d" % n
        subopt = dataframe_RNA.loc[["sequence", id]]
        
        file_id = '_' + family_name + '_subopt_%d' % n

        project = Project(subopt, file_id, n, family_name)
        project.runSNM()

        # project.buildTree(score="variance", percentage=motif_percentage)
        # print (json.dumps(project.tree.tree_format('json')))