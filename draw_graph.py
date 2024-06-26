import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import argparse

def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("-s", "--shortest-path", default="-1")
    parser.add_argument("-i", "--image-size", default="10")
    parser.add_argument("-p", "--path-only", action="store_true")
    args = parser.parse_args()
    size = int(args.image_size)
    node = 0
    chosen_path = []

    adj_list = np.loadtxt(args.filename)

    G = nx.Graph()

    if not bool(args.path_only):
        for i in range(len(adj_list)):
            G.add_node(i)
            for j in range(len(adj_list[i])):

                if adj_list[i][j] > 0:
                    G.add_edge(i, j, weight=adj_list[i][j], edge_color="blue")


        pos = nx.fruchterman_reingold_layout(G, k=12, scale=2)
        plt.figure(1, figsize=(size, size))
        nx.draw_networkx_nodes(G, pos)
        nx.draw_networkx_edges(G, pos, edge_color='black')
        
        if int(args.shortest_path) >= 0:
            paths = np.loadtxt("path.txt", dtype=str, delimiter="\n")
            colored = []

            if len(paths) - 1 < int(args.shortest_path):
                print("Invalid node number. Node number <= ", len(paths) - 1)
                return
            if int(paths[int(args.shortest_path)] == "-1"):
                print("Invalid node number. Node number cannot be equal to source node: ", int(args.shortest_path))
                return
            elif int(paths[int(args.shortest_path)] == "-2"):
                print("Node not connected. Node: ", int(args.shortest_path))
                return
            else:
                node = int(args.shortest_path)
                chosen_path = np.array(list(paths[node].split(" ")), dtype=int)

                for i in range(len(chosen_path) - 1):
                    x = chosen_path[i]
                    y = chosen_path[i+1]

                    colored.append((x, y))

                nx.draw_networkx_edges(G, pos, edgelist=colored, edge_color='red', width=3)
    else:
        paths = np.loadtxt("path.txt", dtype=str, delimiter="\n")
        node = int(args.shortest_path)
        chosen_path = np.array(list(paths[node].split(" ")), dtype=int)

        if(node >= 0):
            for i in range(len(chosen_path) - 1):
                x = chosen_path[i]
                y = chosen_path[i+1]

                G.add_node(x)
                G.add_node(y)
                G.add_edge(x, y, weight=adj_list[x][y])

            pos = nx.fruchterman_reingold_layout(G, k=12, scale=2)
            plt.figure(1, figsize=(size, size))
            nx.draw_networkx_nodes(G, pos)
            nx.draw_networkx_edges(G, pos, edge_color='black')

    nx.draw_networkx_labels(G, pos)
    weights = nx.get_edge_attributes(G, "weight")

    nx.draw_networkx_edge_labels(G, pos, weights)
            


    plt.savefig("graph.png")
    

if __name__ == "__main__":
    main()