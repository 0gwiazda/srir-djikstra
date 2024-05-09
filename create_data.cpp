#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{
    std::cout << sizeof(uint8_t) << "\n";
    int N = std::stoi(argv[1]);
    
    char file_type = argv[2][0];
    std::ofstream file;
    if (file_type == 'b')
    {
        file = std::ofstream(argv[3], std::ios::binary);
    }
    else
    {
        file = std::ofstream(argv[3]);
    }

    uint8_t** matrix = new uint8_t*[N];
    for(int i = 0; i < N; ++i)
        matrix[i] = new uint8_t[N];
    
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
            
            if (std::rand() % 10 < 9)
            {
                matrix[i][j] = 0;
                matrix[j][i] = 0;
            }
            else
            {
                random_value = std::rand() % 255 + 1;
                matrix[i][j] = random_value;
                matrix[j][i] = random_value;
            }
        }
    }

    if (file_type == 'b')
    {
        file.write(reinterpret_cast<const char*>(&N), sizeof(int));
        
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                file.write(reinterpret_cast<const char*>(&matrix[i][j]), sizeof(uint8_t));
            }
        }
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                file << matrix[i][j] << " ";
            }
            file << "\n";
        }
    }

    file.close();

    for(int i = 0; i < N; ++i)
        delete[] matrix[i];
}