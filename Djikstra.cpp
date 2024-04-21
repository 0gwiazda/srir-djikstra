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
    int* l;
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
        std::vector<std::vector<int>> matrix;
        while (std::getline(file, line))
        {
            matrix.push_back(std::vector<int>());
            
            std::istringstream iss(line);
            std::string value;
            while(iss >> value)
            {
                matrix[row_iter].push_back(std::stoi(value));
            }
            matrix[row_iter].shrink_to_fit();
            if (n == 0)
            {
                n = matrix[0].size();
            }
            row_iter++;
        }

        std::cout << "Done\n";
        l = new int[n];
        for (int i = 0; i < n; i++)
        {
            l[i] = matrix[0][i] != 0 ? matrix[0][i] : INT_MAX;
        }

        // determine how to split the matrix
        split_size = n / size;
        remainder = n - split_size * size;

        int split_size_send;
        for (int process = 1; process < size; process++)
        {
            split_size_send = process + remainder < size ? split_size : split_size + 1;
            MPI_Send(&split_size_send, 1, MPI_INT, process, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&split_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        l = new int[n];
    }
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&source, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&l[0], n, MPI_INT, 0, MPI_COMM_WORLD);

    // start calculation
    std::vector<int> V_t({source});

    std::cout << "Process " << rank << ": split size = " << split_size << "\n";
    MPI_Finalize();
    delete[] l;
    return 0;
}