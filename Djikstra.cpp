#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <limits.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size, n = 0, source = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int split_size, remainder;
    std::vector<std::vector<float>> adjacency_matrix;
    std::vector<float> l;
    if (rank == 0)
    {
        if (argc != 2)
        {
            throw std::invalid_argument("File name not specified");
        }

        std::ifstream file(argv[1]);

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
                adjacency_matrix[row_iter].push_back(std::stoi(value));
            }
            adjacency_matrix[row_iter].shrink_to_fit();
            if (n == 0)
            {
                n = adjacency_matrix[0].size();
            }
            row_iter++;
        }

        std::cout << "Done\n";
        for (int i = 0; i < n; i++)
        {
            l.push_back(adjacency_matrix[0][i] != 0.f ? adjacency_matrix[0][i] : std::numeric_limits<float>::infinity());
        }
        l.shrink_to_fit();

        // determine how to split the matrix
        split_size = n / size;
        remainder = n - split_size * size;
        // 1. send size of matrix to all nodes
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int split_size_send, split_sum = split_size;
        for (int process = 1; process < size; process++)
        {
            split_size_send = process + remainder < size ? split_size : split_size + 1;
            // 2. send split size for every node
            MPI_Send(&split_size_send, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
            // 3. send part of l array to every node
            MPI_Send(&l[split_sum], split_size_send, MPI_FLOAT, process, 0, MPI_COMM_WORLD);
            // 4. send part of afjacency matrix to every node
            for (int i = 0; i < split_size_send; i++)
            {
                MPI_Send(&adjacency_matrix[i + split_size_send][0], n, MPI_FLOAT, process, 0, MPI_COMM_WORLD);
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
        for (int i = 0; i < split_size; i++)
        {
            l.push_back(0);
        }
        l.shrink_to_fit();
        MPI_Recv(&l[0], split_size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // 4.
        for (int i = 0; i < split_size; i++)
        {
            adjacency_matrix.push_back(std::vector<float>());
            for (int j = 0; j < n; j++)
            {
                adjacency_matrix[i].push_back(0.f);
            }
            MPI_Recv(&adjacency_matrix[i][0], n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    MPI_Bcast(&source, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // start calculation
    std::vector<int> V_t = { source };

    if (rank == size - 1)
    {
        std::cout << "Process " << rank << ": split size = " << split_size << "\n";
        std::cout << "l : ";
        for (int i = 0; i < split_size; i++)
        {
            std::cout << l[i] << " ";
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
    
    MPI_Finalize();
    return 0;
}