from Clusters import *


def main():
    family_list = ["RF03064", "RF02913", "RF02914", "RF02924", "RF03064"]
    
    num_of_structures = 20
    num_of_sequences = 20

    clusters = Clusters(family_list, num_of_structures, num_of_sequences)
    
    dataFrameSNM = clusters.get_dataFrameSNM()
    print(dataFrameSNM)

    xTrain, xTest, yTrain, yTest = clusters.splitTrainAndTest(0.5, 0.5)

    clusters.k_means(xTrain, xTest, yTrain, yTest)
    clusters.miniBatchKMeans(xTrain, xTest, yTrain, yTest)
    clusters.agglomerativeClustering(xTrain, xTest, yTrain, yTest)
    clusters.birch(xTrain, xTest, yTrain, yTest)
    clusters.dbscan(xTrain, xTest, yTrain, yTest)
    # clusters.optics(xTrain, xTest, yTrain, yTest)
    clusters.meanShift(xTrain, xTest, yTrain, yTest)
    clusters.spectralClustering(xTrain, xTest, yTrain, yTest)
    clusters.gaussianMixture(xTrain, xTest, yTrain, yTest)

    # clusters.affinityPropagation(xTrain, xTest, yTrain, yTest)
    


if __name__ == "__main__":
    main()