#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
    std::ofstream file(argv[1]);

    int N = std::stoi(argv[2]);

    int** matrix = new int*[N];
    for(int i = 0; i < N; ++i)
        matrix[i] = new int[N];
    
    std::srand(std::time(nullptr));
    int random_value;
    for (int i = 0; i < N; i++)
    {
        for (int j = i; j < N; j++)
        {
            if (i == j)
            {
                matrix[i][j] = 0;
                continue;
            }
            
            if (std::rand() % 2 < 1)
            {
                matrix[i][j] = 0;
                matrix[j][i] = 0;
            }
            else
            {
                random_value = std::rand() % 15 + 1;
                matrix[i][j] = random_value;
                matrix[j][i] = random_value;
            }
        }
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            file << matrix[i][j] << " ";
        }
        file << "\n";
    }

    file.close();

    for(int i = 0; i < N; ++i)
        delete[] matrix[i];
}