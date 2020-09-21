from cluster import *


def main():
    family_name = "RF03064" # RF02913 RF02914 RF02924 RF03064
    
    num_of_clusters = 25
    num_of_sequences = 20

    generate_cluster(family_name, num_of_clusters, num_of_sequences)


if __name__ == "__main__":
    main()