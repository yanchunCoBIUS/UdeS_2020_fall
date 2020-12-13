from compareCluster import compare_clusters
from visualize_clusters import *
from singleRunCompareCluster import *
import pickle


def main():
    # 50 families
    all_family_list = ["RF02913", "RF02914", "RF02924", "RF03064", "RF02440", "RF01695", "RF03161", "RF01739", "RF03158", "RF03156",
                       "RF03106", "RF03157", "RF01704", "RF03144", "RF02929", "RF03072", "RF02987", "RF03169", "RF02958", "RF03130",
                       "RF03142", "RF02957", "RF01717", "RF01725", "RF00174", "RF03132", "RF02001", "RF03054", "RF03131", "RF03135",
                       "RF03127", "RF03000", "RF02035", "RF03069", "RF02969", "RF03012", "RF03086", "RF03032", "RF01419", "RF03087",
                       "RF01734", "RF03057", "RF02925", "RF01701", "RF03074", "RF00169", "RF03141", "RF03065", "RF02975", "RF02999"]
        
    #3 4 5 structures
    num_of_structures = 5

     #3 4 5 sequences
    num_of_sequences = 5
    
    # results = []
    # result = singleRunCompareCluster(all_family_list, num_of_structures, num_of_sequences)
    # results.append(result)
    # with open('result_50familiy_3structure_3sequence.data', 'wb') as filehandle:
    #     # store the data as binary data stream
    #     pickle.dump(results, filehandle)
    
    # result_ARI_dict = dict()

    # compare_clusters(all_family_list, num_of_structures, num_of_sequences)
    visulize()

if __name__ == "__main__":
    main()