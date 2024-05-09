#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <limits>
#include <mpi.h>

#define INFINITY std::numeric_limits<float>::infinity()
#define END MPI_Finalize(); return 0;

struct node_distance
{
    float distance;
    int node;
};

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size, n = 0, source = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int split_size, remainder, starting_node = 0;
    std::vector<std::vector<float>> adjacency_matrix;
    //std::vector<float> l;
    std::vector<int> recvcounts(size), displs(size);

    double start = MPI_Wtime();

    if (rank == 0)
    {
        if (argc < 2)
        {
            throw std::invalid_argument("File name not specified");
        }

        std::ifstream file(argv[1]);
        
        if(argc != 3)
        {
            std::cout << "No node specified, using default: " + source << std::endl;
        }
        else
        {
            source = std::stoi(argv[2]);
        }

        if (!file.is_open())
        {
            throw std::invalid_argument("Could not open the file");
        }

        std::cout << "File \"" << argv[1] << "\" opened, parsing...\n" ;

        std::string line;
        int row_iter = 0;
        while (std::getline(file, line))
        {
            adjacency_matrix.push_back(std::vector<float>());
            
            std::istringstream iss(line);
            std::string value;
            while(iss >> value)
            {
                adjacency_matrix[row_iter].push_back(std::stoi(value) > 0 ? std::stoi(value) : INFINITY);
            }
            adjacency_matrix[row_iter].shrink_to_fit();
            if (n == 0)
            {
                n = adjacency_matrix[0].size();
            }
            row_iter++;
        }

        std::cout << "Done\n";/*
        for (int i = 0; i < n; i++)
        {
            l.push_back(adjacency_matrix[0][i] != 0.f ? adjacency_matrix[0][i] : INFINITY);
        }
        l.shrink_to_fit();*/

        // determine how to split the matrix
        split_size = n / size;
        remainder = n - split_size * size;
        // 1. send size of matrix to all nodes
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        recvcounts[0] = split_size;
        displs[0] = 0;

        int split_size_send, split_sum = split_size;
        for (int process = 1; process < size; process++)
        {
            split_size_send = process + remainder < size ? split_size : split_size + 1;

            recvcounts[process] = split_size_send;
            displs[process] = split_sum;
            
            // 2. send split size for every node
            MPI_Send(&split_size_send, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            // 3. send starting node for every node
            MPI_Send(&split_sum, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            // 4. send part of l array to every node
            //MPI_Send(&l[split_sum], split_size_send, MPI_FLOAT, process, 0, MPI_COMM_WORLD);
            // 5. send part of adjacency matrix to every node
            for (int i = 0; i < split_size_send; i++)
            {
                MPI_Send(&adjacency_matrix[i + split_sum][0], n, MPI_FLOAT, process, 0, MPI_COMM_WORLD);
            }

            split_sum += split_size_send;
        }
    }
    else
    {
        // 1.
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // 2.
        MPI_Recv(&split_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // 3.
        MPI_Recv(&starting_node, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // 4.
        //l.resize(split_size);
        //MPI_Recv(&l[0], split_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // 5.
        for (int i = 0; i < split_size; i++)
        {
            adjacency_matrix.push_back(std::vector<float>(n));
            MPI_Recv(&adjacency_matrix[i][0], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    MPI_Bcast(&source, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<float> local_distance(split_size);
    for (int i = 0; i < split_size; i++)
    {
        local_distance[i] = adjacency_matrix[i][source];
    }
    
    // check if data has been send successfuly
    if (rank == size - 1)
    {
        std::cout << "Process " << rank << ": split size = " << split_size << "\n";
        std::cout << "l : ";
        for (int i = 0; i < split_size; i++)
        {
            std::cout << local_distance[i] << " ";
        }
        std::cout << "\nmatrix\n";
        for (int i = 0; i < split_size; i++)
        {
            for (int j = 0; j < n; j++)
            {
                std::cout << adjacency_matrix[i][j] << " ";
            }
            std::cout << "\n";
        }
    }
    
    // start calculation
    std::vector<int> local_predecessor_node(split_size, source);
    std::vector<bool> node_visited(split_size);
    node_distance local_min, global_min;

    if (source >= starting_node && source < starting_node + split_size)
    {
        node_visited[source - starting_node] = true;
    }
    
    int local_min_node;
    float shortest_distance, new_distance;
    for (int i = 0; i < n - 1; i++)
    {
        // find min local distance vertex
        local_min_node = -1;
        shortest_distance = INFINITY;
        
        for (int node = 0; node < split_size; node++) 
        {
            if (!node_visited[node] && (local_distance[node] < shortest_distance))
            {
                shortest_distance = local_distance[node];
                local_min_node = node;
            }
        }
        
        if (local_min_node != -1)
        {
            local_min.distance = local_distance[local_min_node];
            local_min.node = local_min_node + starting_node;
        }
        else
        {
            local_min.distance = INFINITY;
            local_min.node = -1;
        }
        
        MPI_Allreduce(&local_min, &global_min, 1, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);
        
        if (global_min.node == -1)
            break;

        if (global_min.node >= starting_node && global_min.node < starting_node + split_size)
        {
            node_visited[global_min.node - starting_node] = true;
        }
        
        for (int node = 0; node < split_size; node++)
        {
            if (!node_visited[node])
            {
                new_distance = global_min.distance + adjacency_matrix[node][global_min.node];
                if (new_distance < local_distance[node])
                {
                    local_distance[node] = new_distance;
                    local_predecessor_node[node] = global_min.node;
                }
            }
        }
    }

    std::vector<float> global_distance(n);
    std::vector<int> global_predecessor_node(n);
    MPI_Gatherv(&local_distance[0], split_size, MPI_FLOAT, &global_distance[0], &recvcounts[0], &displs[0], MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(&local_predecessor_node[0], split_size, MPI_INT, &global_predecessor_node[0], &recvcounts[0], &displs[0], MPI_FLOAT, 0, MPI_COMM_WORLD);
    std::cout << rank << "\n";
    
    if (rank == 0)
    {
        std::cout << "FINAL RESULTS\n";

        for (int i = 0; i < n; i++)
        {
            std::cout << global_distance[i] << " ";
        }
        std::cout << "\npredecessor ";

        for (int i = 0; i < n; i++)
        {
            std::cout << global_predecessor_node[i] << " ";
        }
        std::cout << "\n";
        
        std::ofstream Path("path.txt");

        for (int node = 0; node < n; node++)
        {
            if (node == source)
                continue;

            std::cout << "Path to node " << node << ": ";
            int current_node = node;

            while (current_node != source)
            {
                std::cout << current_node << " ";
                Path << current_node << " ";
                current_node = global_predecessor_node[current_node];
            } 
            std::cout << "\n";

            Path << source << "\n";
        }

        double end = MPI_Wtime();
        Path.close();
        std::cout << "Elapsed time: " << end - start << std::endl;


    }

    MPI_Finalize();

    return 0;
}