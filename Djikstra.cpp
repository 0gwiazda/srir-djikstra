#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <limits>
#include <chrono>
#include <mpi.h>

#define INFINITY std::numeric_limits<float>::infinity()

struct node_distance
{
    float distance;
    int node;
};

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size, n = 0, source = 0, sink;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int split_size, remainder, starting_node = 0;
    std::vector<std::vector<float>> adjacency_matrix;
    std::vector<int> recvcounts(size), displs(size);
    char file_type;
    if (rank == 0)
    {
        if (argc != 5)
        {
            throw std::logic_error("Not enough arguments! Usage: djikstra [file type(b/t)] [file name] [source] [sink]");
        }

        file_type = *argv[1];

        std::ifstream file;
        if (file_type == 'b')
        {
            file = std::ifstream(argv[2], std::ios::binary);
        }
        else
        {
            file = std::ifstream(argv[2]);
        }

        if (!file.is_open())
        {
            throw std::logic_error("Could not open the file!");
        }

        std::cout << "File \"" << argv[2] << "\" opened, parsing...\n";

        try
        {
            source = std::stoi(argv[3]);
            sink = std::stoi(argv[4]);
        }
        catch(...)
        {
            throw std::logic_error("Source and sink need to be integer numbers");
        }

        if (file_type == 'b')
        {
            int row_iter = 0;
            uint8_t value;
            file.read((char *) &n, sizeof(int));
            for (int row = 0; row < n; row++)
            {
                adjacency_matrix.push_back(std::vector<float>(n));
                for (int col = 0; col < n; col++)
                {
                    file.read((char *) &value, sizeof(uint8_t));
                    adjacency_matrix[row][col] = (value > 0) || (row == col) ? value : INFINITY;
                }
                //std::cout << '\r' << (static_cast<float>(row) / n * 100) << '%';
            }
        }
        else
        {
            std::string line;
            int row_iter = 0;
            while (std::getline(file, line))
            {
                adjacency_matrix.push_back(std::vector<float>());
            
                std::istringstream iss(line);
                std::string value;
                int col_iter = 0;
                while(iss >> value)
                {
                    adjacency_matrix[row_iter].push_back((std::stoi(value) > 0) || (row_iter == col_iter) ? std::stoi(value) : INFINITY);
                    col_iter++;
                }
                adjacency_matrix[row_iter].shrink_to_fit();
                if (n == 0)
                {
                    n = adjacency_matrix[0].size();
                }
                row_iter++;
            }
        }

        std::cout << "\rDone    \n\nCalculating the shortest path...\n";
        
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
            // 4. send part of adjacency matrix to every node
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
        for (int i = 0; i < split_size; i++)
        {
            adjacency_matrix.push_back(std::vector<float>(n));
            MPI_Recv(&adjacency_matrix[i][0], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    MPI_Bcast(&source, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sink, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<float> local_distance(split_size);
    for (int i = 0; i < split_size; i++)
    {
        local_distance[i] = adjacency_matrix[i][source];
    }
    
    // start calculation
    std::vector<int> local_predecessor_node(split_size, source);
    std::vector<bool> node_visited(split_size);
    node_distance local_min, global_min;

    if (source >= starting_node && source < starting_node + split_size)
    {
        node_visited[source - starting_node] = true;
        local_distance[source - starting_node] = 0;
        local_predecessor_node[source - starting_node] = source;
    }
    
    int local_min_node;
    float shortest_distance, new_distance;
    double start = MPI_Wtime();
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
    double end = MPI_Wtime();
    if (rank == 0)
    {
        if(sink != -1)
        {
            std::cout << "\nPath from node " << source << " to node " << sink << ":\n\n";
                int current_node = sink;
                int loop_detector = n;
                while (true)
                {
                    if (global_distance[sink] == INFINITY)
                    {
                        std::cout << "Nodes not connected\n\n";
                        break;
                    }
                    
                    std::cout << '(' << current_node << ")<-";
                    current_node = global_predecessor_node[current_node];
                    loop_detector--;

                    if (current_node == source)
                    {
                        std::cout << '(' << current_node << ")\n\n";
                        break;
                    }
                    
                    if (loop_detector < 0)
                    {
                        std::cout << "LOOP DETECTED";
                        break;
                    }
                }
                std::cout << "Path length: " << global_distance[sink] << "\n";
        }
        else
        {
            std::cout << "Saving to file \"paths.txt\"..." << std::endl;
            std::ofstream Path("path.txt");

            for (int node = 0; node < n; node++)
            {
                if (node == source)
                {
                    Path << -1 << "\n";
                    continue;
                }

                if (global_distance[node] == INFINITY)
                {
                    Path << -2 << std::endl;
                    break;
                }

                int current_node = node;
                int loop_detector = n;
                
                while (true)
                {
                    Path << current_node << " ";
                    current_node = global_predecessor_node[current_node];
                    loop_detector--;

                    if (current_node == source)
                    {
                        Path << current_node << '\n';
                        break;
                    }
                    
                    if (loop_detector < 0)
                    {
                        std::cout << "LOOP DETECTED";
                        break;
                    }
                }
            }

            Path.close();
        }
    }

    /*
    if (rank == 0)
        std::cout << "Time: " << (end - start) << "[s]\n";*/

    MPI_Finalize();

    return 0;
}