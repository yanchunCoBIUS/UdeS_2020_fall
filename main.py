from Clusters import *
from sklearn.cluster import KMeans
from sklearn.metrics import accuracy_score


def main():
    family_list = ["RF03064", "RF02913", "RF02914", "RF02924", "RF03064"]
    
    num_of_structures = 20
    num_of_sequences = 20

    clusters = Clusters(family_list, num_of_structures, num_of_sequences)
    
    dataFrameSNM = clusters.get_dataFrameSNM()
    print(dataFrameSNM)

    xTrain, xTest, yTrain, yTest = clusters.splitTrainAndTest(0.5, 0.5)

    kmeans = KMeans(n_clusters=5, random_state=0).fit(xTrain)
    # yTrain_label = kmeans.labels_
    # print(yTrain_label)
    # print(np.array(yTrain['Family']))

    yTest_predict = kmeans.predict(xTest)
    print("Predict: \n", yTest_predict)
    print("Actual: \n", np.array(yTest['Family']))

    score = accuracy_score(np.array(yTest['Family']), yTest_predict)
    print(score)


if __name__ == "__main__":
    main()