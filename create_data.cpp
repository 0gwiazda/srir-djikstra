#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstring>

int main(int argc, char* argv[])
{
    int N = std::stoi(argv[1]);
    float chance = std::atof(argv[2]);
    
    bool text_output = false, binary_output = false;
    if (strlen(argv[3]) > 1)
    {
        if (argv[3][0] == 't' || argv[3][1] == 't')
            text_output = true;
        if (argv[3][0] == 'b' || argv[3][1] == 'b')
            binary_output = true;
    }
    else
    { 
        if (argv[3][0] == 'b')
            binary_output = true;
        else
            text_output = true;
    }

    std::ofstream file_text, file_binary;
    
    if (binary_output)
    {
        file_binary = std::ofstream(std::string(argv[4]) + ".bin", std::ios::binary);
    }

    if (text_output)
    {
        file_text = std::ofstream(std::string(argv[4]) + ".txt");
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
            
            if (static_cast<float>(std::rand()) / RAND_MAX >= chance)
            {
                matrix[i][j] = 0;
                matrix[j][i] = 0;
            }
            else
            {
                random_value = std::rand() % 100 + 1;
                matrix[i][j] = random_value;
                matrix[j][i] = random_value;
            }
        }
    }

    if (binary_output)
    {
        file_binary.write(reinterpret_cast<const char*>(&N), sizeof(int));
        
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                file_binary.write(reinterpret_cast<const char*>(&matrix[i][j]), sizeof(uint8_t));
            }
        }
    }

    if (text_output)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                file_text << static_cast<int>(matrix[i][j]) << " ";
            }
            file_text << "\n";
        }
    }

    if (binary_output)
        file_binary.close();

    if (text_output)
        file_text.close();

    for(int i = 0; i < N; ++i)
        delete[] matrix[i];
    delete[] matrix;
}