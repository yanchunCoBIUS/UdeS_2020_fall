from Clusters import *
from NewClusters import *


def main():
    family_list = ["RF02913", "RF02914", "RF02924", "RF03064"]
    
    num_of_structures = 5
    num_of_sequences = 5

    clusters = Clusters(family_list, num_of_structures, num_of_sequences)
    dataFrameSNM = clusters.get_dataFrameSNM()
    # print(dataFrameSNM)

    ###
    # dict_allFamilies['RF02913']['URS0000D66AA1_12908/1-115']: dict_keys([
    # 'sequence', 
    # '>RF02913#URS0000D66AA1_12908/1-115#subopt0', 
    # '>RF02913#URS0000D66AA1_12908/1-115#subopt1', 
    # '>RF02913#URS0000D66AA1_12908/1-115#subopt2', 
    # '>RF02913#URS0000D66AA1_12908/1-115#subopt3',
    # '>RF02913#URS0000D66AA1_12908/1-115#subopt4'
    # ])
    ###
    dict_allFamilies = clusters.get_dict_allFamilies()

    x, y = clusters.getFeaturesAndTarget()

    # label_birch = clusters.birch(x, y)
    # new_clusters, new_x, new_y = get_labeled_cluster(dict_allFamilies, dataFrameSNM, family_list, num_of_structures, num_of_sequences, label_birch)
    # new_clusters.birch(new_x, new_y)
    # new_clusters.agglomerativeClustering(new_x, new_y)
    # new_clusters.k_means(new_x, new_y)
    # new_clusters.miniBatchKMeans(new_x, new_y)
    # new_clusters.dbscan(new_x, new_y)
    # new_clusters.meanShift(new_x, new_y)
    # new_clusters.spectralClustering(new_x, new_y)
    # new_clusters.gaussianMixture(new_x, new_y)

    label_agglomerative = clusters.agglomerativeClustering(x, y)
    new_clusters, new_x, new_y = get_labeled_cluster(dict_allFamilies, dataFrameSNM, family_list, num_of_structures, num_of_sequences, label_agglomerative)
    new_clusters.birch(new_x, new_y)
    new_clusters.agglomerativeClustering(new_x, new_y)
    new_clusters.k_means(new_x, new_y)
    new_clusters.miniBatchKMeans(new_x, new_y)
    new_clusters.dbscan(new_x, new_y)
    new_clusters.meanShift(new_x, new_y)
    new_clusters.spectralClustering(new_x, new_y)
    new_clusters.gaussianMixture(new_x, new_y)

    # label_kmeans = clusters.k_means(x, y)
    # new_clusters, new_x, new_y = get_labeled_cluster(dict_allFamilies, dataFrameSNM, family_list, num_of_structures, num_of_sequences, label_kmeans)
    # new_clusters.birch(new_x, new_y)
    # new_clusters.agglomerativeClustering(new_x, new_y)
    # new_clusters.k_means(new_x, new_y)
    # new_clusters.miniBatchKMeans(new_x, new_y)
    # new_clusters.dbscan(new_x, new_y)
    # new_clusters.meanShift(new_x, new_y)
    # new_clusters.spectralClustering(new_x, new_y)
    # new_clusters.gaussianMixture(new_x, new_y)

    # label_miniBatchKMeans = clusters.miniBatchKMeans(x, y)
    # new_clusters, new_x, new_y = get_labeled_cluster(dict_allFamilies, dataFrameSNM, family_list, num_of_structures, num_of_sequences, label_miniBatchKMeans)
    # new_clusters.birch(new_x, new_y)
    # new_clusters.agglomerativeClustering(new_x, new_y)
    # new_clusters.k_means(new_x, new_y)
    # new_clusters.miniBatchKMeans(new_x, new_y)
    # new_clusters.dbscan(new_x, new_y)
    # new_clusters.meanShift(new_x, new_y)
    # new_clusters.spectralClustering(new_x, new_y)
    # new_clusters.gaussianMixture(new_x, new_y)

    # label_dbscan = clusters.dbscan(x, y)
    # new_clusters, new_x, new_y = get_labeled_cluster(dict_allFamilies, dataFrameSNM, family_list, num_of_structures, num_of_sequences, label_dbscan)
    # new_clusters.birch(new_x, new_y)
    # new_clusters.agglomerativeClustering(new_x, new_y)
    # new_clusters.k_means(new_x, new_y)
    # new_clusters.miniBatchKMeans(new_x, new_y)
    # new_clusters.dbscan(new_x, new_y)
    # new_clusters.meanShift(new_x, new_y)
    # new_clusters.spectralClustering(new_x, new_y)
    # new_clusters.gaussianMixture(new_x, new_y)

    # label_meanShift = clusters.meanShift(x, y)
    # new_clusters, new_x, new_y = get_labeled_cluster(dict_allFamilies, dataFrameSNM, family_list, num_of_structures, num_of_sequences, label_meanShift)
    # new_clusters.birch(new_x, new_y)
    # new_clusters.agglomerativeClustering(new_x, new_y)
    # new_clusters.k_means(new_x, new_y)
    # new_clusters.miniBatchKMeans(new_x, new_y)
    # new_clusters.dbscan(new_x, new_y)
    # new_clusters.meanShift(new_x, new_y)
    # new_clusters.spectralClustering(new_x, new_y)
    # new_clusters.gaussianMixture(new_x, new_y)

    # label_spectralClustering = clusters.spectralClustering(x, y)
    # new_clusters, new_x, new_y = get_labeled_cluster(dict_allFamilies, dataFrameSNM, family_list, num_of_structures, num_of_sequences, label_spectralClustering)
    # new_clusters.birch(new_x, new_y)
    # new_clusters.agglomerativeClustering(new_x, new_y)
    # new_clusters.k_means(new_x, new_y)
    # new_clusters.miniBatchKMeans(new_x, new_y)
    # new_clusters.dbscan(new_x, new_y)
    # new_clusters.meanShift(new_x, new_y)
    # new_clusters.spectralClustering(new_x, new_y)
    # new_clusters.gaussianMixture(new_x, new_y)

    # label_gaussianMixture = clusters.gaussianMixture(x, y)
    # new_clusters, new_x, new_y = get_labeled_cluster(dict_allFamilies, dataFrameSNM, family_list, num_of_structures, num_of_sequences, label_gaussianMixture)
    # new_clusters.birch(new_x, new_y)
    # new_clusters.agglomerativeClustering(new_x, new_y)
    # new_clusters.k_means(new_x, new_y)
    # new_clusters.miniBatchKMeans(new_x, new_y)
    # new_clusters.dbscan(new_x, new_y)
    # new_clusters.meanShift(new_x, new_y)
    # new_clusters.spectralClustering(new_x, new_y)
    # new_clusters.gaussianMixture(new_x, new_y)
  


def get_labeled_cluster(dict_allFamilies, dataFrameSNM, family_list, num_of_structures, num_of_sequences, clustered_label):
    """
    sequence[i] example:
        RF03064#CP001107.1/991894-991966#subopt3
    splited_family_sequence_structure = sequence[i].split('#')
        splited_family_sequence_structure[0]: 'RF03064'
        splited_family_sequence_structure[1]: 'CP001107.1/991894-991966'
        splited_family_sequence_structure[2]: 'subopt3'
    """
    dict_labeled = dict()
    sequence = dataFrameSNM.iloc[:, 0]
    list_cluster_count = []
    rna_sequence = []
    
    for i in range(len(sequence)):
        splited_family_sequence_structure = sequence[i].split('#')
        ###
        # splited_family_sequence_structure[0]: Family name
        # splited_family_sequence_structure[1]: Sequence name
        ###
        family_plus_sequence = 'Family#%d#%s' % (family_list.index(splited_family_sequence_structure[0]), splited_family_sequence_structure[1])
        ###
        # splited_family_sequence_structure[2]: subopt + number
        ###
        sub_structure = splited_family_sequence_structure[2]
        if (i % num_of_structures == 0):
            dict_labeled[family_plus_sequence] = dict()
            dict_label_count = dict()
            '''
            11/17 Need to improve: too slow
            '''
            for i in range(len(family_list)):
                dict_label_count[i] = 0
        if (clustered_label[i] in dict_label_count):
            dict_label_count[clustered_label[i]] += 1
        else:
            dict_label_count[clustered_label[i]] = 1
        dict_labeled[family_plus_sequence][sub_structure] = clustered_label[i]
        if (i != 0 and i % (num_of_structures - 1) == 0):
            list_cluster_count.append(dict_label_count)

    df_labeled = pd.DataFrame.from_dict(dict_labeled, orient='index')
    print(df_labeled, '\n')

    new_df = pd.concat([pd.DataFrame(dict_labeled.keys()), pd.DataFrame(list_cluster_count)], axis=1)
    print(new_df)

    one_subopt_dict, input_for_SNM = find_which_subopt_to_use(dict_allFamilies, family_list, dict_labeled, rna_sequence)
    print(pd.DataFrame.from_dict(one_subopt_dict, orient='index'))

    print(input_for_SNM)
    runSNM(input_for_SNM)


    new_clusters = NewClusters(family_list, num_of_structures, num_of_sequences)
    new_x, new_y = new_clusters.getFeaturesAndTarget()
    return new_clusters, new_x, new_y


def find_which_subopt_to_use(dict_allFamilies, family_list, dict_labeled, rna_sequence):
    '''
    dict_allFamilies['RF02913']['URS0000D66AA1_12908/1-115']: dict_keys([
    'sequence', 
    '>RF02913#URS0000D66AA1_12908/1-115#subopt0', 
    '>RF02913#URS0000D66AA1_12908/1-115#subopt1', 
    '>RF02913#URS0000D66AA1_12908/1-115#subopt2', 
    '>RF02913#URS0000D66AA1_12908/1-115#subopt3',
    '>RF02913#URS0000D66AA1_12908/1-115#subopt4'
    ])
    '''
    dot_bracket_str = ""
    one_subopt_dict = dict()
    for family_plus_sequence, subopt_dict in dict_labeled.items():
        splited_family_sequence = family_plus_sequence.split('#')
        ###
        # splited_family_sequence: 'Family#0#URS0000D66AA1_12908/1-115'
        # splited_family_sequence[0]: 'Family'
        # splited_family_sequence[1]: Family number: '0' (index in 'family_list')
        # splited_family_sequence[2]: sequence: 'URS0000D66AA1_12908/1-115'
        ###
        subopt_list = list(subopt_dict.keys())
        # print(subopt_list[0])
        family_name = family_list[int(splited_family_sequence[1])]
        id = splited_family_sequence[2]
        find = False
        for subopt in subopt_list:
            # print("iterate: %s" % subopt)
            if (str(subopt_dict[subopt]) == splited_family_sequence[1]):
                find = True
                # print("match: %s" % subopt_dict[subopt])
                # print("match: %s" % splited_family_sequence[1])

                one_subopt_dict[family_plus_sequence] = subopt                
                id_subopt = ">%s#%s#%s" % (family_name, id, subopt)
                structure = dict_allFamilies[family_name][id]['%s' % id_subopt]
                dot_bracket_str = dot_bracket_str + (id_subopt + '\n' +
                                dict_allFamilies[family_name][id]['sequence'] + '\n' +
                                structure + '\n'
                                )
                # rna_sequence.append(dot_bracket_str)
                break
        if (find == False):   
            one_subopt_dict[family_plus_sequence] = 'subopt0'
            id_subopt = ">%s#%s#subopt%s" % (family_name, id, '0')
            structure = dict_allFamilies[family_name][id]['%s' % id_subopt]
            dot_bracket_str = dot_bracket_str + (id_subopt + '\n' +
                            dict_allFamilies[family_name][id]['sequence'] + '\n' +
                            structure + '\n'
                            )
    return one_subopt_dict, dot_bracket_str



if __name__ == "__main__":
    main()