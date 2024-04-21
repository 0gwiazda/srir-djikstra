#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <limits.h>


int main(int argc, char* argv[])
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

        

        
    return 0;
}
