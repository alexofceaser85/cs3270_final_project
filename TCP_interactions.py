map = {
    "A" : ("C"),
    "B" : ("C"),
    "C" : ("A", "B", "D", "E"),
    "D" : ("C", "F"),
    "E" : ("C", "F"),
    "F" : ("D", "E", "G", "H"),
    "H" : ("F"),
    "G" : ("F")
}

data_interactions = {}

def traverse_nodes(current_node, maximum_branches):

    paths = []
    immediate_neighbors = []
    for neighbour in map[current_node]:
        immediate_neighbors.append(neighbour)

    for node in immediate_neighbors:
        for secondary_neighbor in map[node]:
            if (current_node != node and current_node != secondary_neighbor and secondary_neighbor != node and secondary_neighbor != current_node):
                paths.append([current_node, node, secondary_neighbor])

    return paths

def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

if __name__ == '__main__':
    visited_nodes = []
    maximum_branches = 2
    node_paths = {}
    for node in map:
        node_paths[node] = traverse_nodes(node, 2)

    for node in node_paths:
        counter = 0
        list_length = len(node_paths[node])

        for node_path in node_paths[node]:
            current_path = node_path

            for jaccard_value in node_paths[node]:
                if jaccard_value != current_path:
                    print(str(current_path) + "   " + str(jaccard_value) + "   " + str(jaccard(current_path, jaccard_value)))
        print("_________________")
