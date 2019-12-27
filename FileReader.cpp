#include <string>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

int main()
{
    std::string path = "fastq";
    for (const auto& entry : fs::directory_iterator(path)) {
        std::string fileName = fs::path(entry.path()).filename().string();
        if (fileName.rfind("J", 0) == 0) {
            std::cout << fs::path(entry.path()).filename() << std::endl;
        }
    }
}

