#include <string>
#include <iostream>
#include <filesystem>
#include <sstream>
#include <fstream>

namespace fs = std::filesystem;

int main()
{
    std::string path = "fastq";
    for (const auto& entry : fs::directory_iterator(path)) {
        std::string fileName = fs::path(entry.path()).filename().string();
        if (fileName.rfind("J", 0) == 0) {
            //std::string pathF = entry.path().string();
            //std::cout << fs::path(entry.path()).filename() << std::endl;
            std::cout << pathF << std::endl;
            std::ifstream fileOpen(entry.path().string());
            if (!fileOpen.is_open()) {
                std::perror("File opening failed");
                return EXIT_FAILURE;
            }
            else
            {
                std::cout << "Otvoren file\n";
            }
            
        }
    }
}

