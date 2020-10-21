from cluster import *


def main():
    family_list = ["RF03064", "RF02913", "RF02914", "RF02924", "RF03064"]
    
    num_of_structures = 2
    num_of_sequences = 2

    generate_cluster(family_list, num_of_structures, num_of_sequences)


if __name__ == "__main__":
    main()