#include <string>
#include <iostream>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

int readFile(std::string* path, bool read_collect, int *lengthValue); //if false then just read, else collect valid chains
int findMostCommonLength(std::vector<unsigned int>* allLengthsVector);
void collectChains(std::string* path, int* length);
int findMinimumDistance(std::string chain1, std::string chain2, int chainLength);

std::vector<unsigned int> allLengthsVector;
std::vector<std::string> validChains;

int main()
{
    int lengthValue = 0;
    std::string path = "fastq";
  
    std::string gene1 = "AGGGCT";
    std::string gene2 = "AGGCAT";
    findMinimumDistance(gene1, gene2, 6);

    readFile(&path, false, &lengthValue);
    
    //std::cout << allLengthsVector.size() << std::endl;

    lengthValue = findMostCommonLength(&allLengthsVector);
    readFile(&path, true, &lengthValue);

    std::cout << "Validnih lanaca ima " << validChains.size() <<std::endl;
     
    return 0;
}

int readFile(std::string *path, bool read_collect, int* lengthValue) {
    //std::string path = "fastq";
    for (const auto& entry : fs::directory_iterator(*path)) {
        std::string fileName = fs::path(entry.path()).filename().string();
        if (fileName.rfind("J", 0) == 0) {
            std::ifstream fileOpen(entry.path().string());
            if (!fileOpen.is_open()) {
                std::perror("File opening failed");
                return EXIT_FAILURE;
            }
            else
            {
                std::cout << "Reading: " << fileName<<"\n";
                std::string line;
                bool flag = false;
                while (fileOpen.good())
                {
                    getline(fileOpen, line);
                    if (flag == true) {
                        unsigned int lineLength = line.length();
           
                        if (read_collect && lineLength == *lengthValue) {
                            validChains.push_back(line);
                        }
                        else {
                            allLengthsVector.push_back(line.length());
                        }
                        flag = false;
                    }
                    if (line.rfind("@16WBS", 0) == 0) {
                        flag = true;        //sljedeća linija je lanac
                    }
                }
                fileOpen.close();
                if (!read_collect) {
                    std::cout << allLengthsVector.size() << std::endl;
                }
            }

        }
    }
    return 0;
}

int findMostCommonLength(std::vector<unsigned int> *allLengthsVector) {
    std::vector <unsigned int> lengthVector; //izvojene duljine lanaca 

    int maxRepetition = 0;
    int lengthValue;   //na kojem mjestu unutar lengthVector-a se nalazi duljina s najvise ponavljanja 
    int flag = 1;

    for (int i = 0; i < (*allLengthsVector).size();i++) {   //izdvoji duljine lanaca bez ponavljanja
        flag = 1;
        for (int j = 0; j < lengthVector.size();j++) {
            if (lengthVector.at(j) == (*allLengthsVector).at(i))
                flag = 0;
        }
        if (flag == 1) {
            lengthVector.push_back((*allLengthsVector).at(i));
        }
    }
    std::cout << "Ukupno ima: " <<lengthVector.size() <<" razlicitih duljina lanaca"<< std::endl;

    //prebroji koliko ponavljanja ima za svaku duljinu
    for (int i = 0; i < lengthVector.size();i++) {      
        int mycount = std::count((*allLengthsVector).begin(), (*allLengthsVector).end(), lengthVector.at(i));
        if (mycount > maxRepetition) {     //odredi maksimum
            maxRepetition = mycount;
            lengthValue = lengthVector.at(i);
        }
    }
    std::cout << "\nNajveci broj pojavljivanja je za: " << lengthValue << std::endl;
    return lengthValue;
}

void collectChains(std::string* path, int* length) {

}

int findMinimumDistance(std::string chain1, std::string chain2, int chainLength) {
    int misMatchPenalty = 3;
    int gapPenalty = 2;

    int i, j; // intialising variables  

    // table for storing optimal substructure answers 
    int** dp = new int* [2 * chainLength + 1];
    for (int i = 0; i < 2 * chainLength + 1; ++i)
        dp[i] = new int[2 * chainLength + 1];

    for (int i = 0; i < 2 * chainLength + 1; i++) {
        for (int j = 0; j < 2 * chainLength + 1; j++) {
            dp[i][j] = 0;
        }
    }
    //std::cout << "kraj";
    // intialising the table 
    for (i = 0; i <= (chainLength + chainLength); i++)
    {
        dp[i][0] = i * gapPenalty;
        dp[0][i] = i * gapPenalty;
    }

    // calcuting the minimum penalty 
    for (i = 1; i <= chainLength; i++)
    {
        for (j = 1; j <= chainLength; j++)
        {
            if (chain1[i - 1] == chain2[j - 1])
            {
                dp[i][j] = dp[i - 1][j - 1];
            }
            else
            {
                dp[i][j] = std::min({ dp[i - 1][j - 1] + misMatchPenalty ,
                                dp[i - 1][j] + gapPenalty ,
                                dp[i][j - 1] + gapPenalty });
            }
        }
    }

    // Reconstructing the solution 
    int l = 2 * chainLength; // maximum possible length 

    i = chainLength; j = chainLength;

    int xpos = l;
    int ypos = l;

    //const int dulj = 11; //6+5

    // Final answers for the respective strings 
    int* yans = new int[l + 1];
    int* xans = new int[l + 1];

    while (!(i == 0 || j == 0))
    {
        if (chain1[i - 1] == chain2[j - 1])
        {
            xans[xpos--] = (int)chain1[i - 1];
            yans[ypos--] = (int)chain2[j - 1];
            i--; j--;
        }
        else if (dp[i - 1][j - 1] + misMatchPenalty == dp[i][j])
        {
            xans[xpos--] = (int)chain1[i - 1];
            yans[ypos--] = (int)chain2[j - 1];
            i--; j--;
        }
        else if (dp[i - 1][j] + gapPenalty == dp[i][j])
        {
            xans[xpos--] = (int)chain1[i - 1];
            yans[ypos--] = (int)'_';
            i--;
        }
        else if (dp[i][j - 1] + gapPenalty == dp[i][j])
        {
            xans[xpos--] = (int)'_';
            yans[ypos--] = (int)chain2[j - 1];
            j--;
        }
    }
    while (xpos > 0)
    {
        if (i > 0) xans[xpos--] = (int)chain1[--i];
        else xans[xpos--] = (int)'_';
    }
    while (ypos > 0)
    {
        if (j > 0) yans[ypos--] = (int)chain2[--j];
        else yans[ypos--] = (int)'_';
    }

    // Since we have assumed the answer to be n+m long, 
    // we need to remove the extra gaps in the starting 
    // id represents the index from which the arrays 
    // xans, yans are useful 
    int id = 1;
    for (i = l; i >= 1; i--)
    {
        if ((char)yans[i] == '_' && (char)xans[i] == '_')
        {
            id = i + 1;
            break;
        }
    }

    // Printing the final answer 
    std::cout << "Minimum Penalty in aligning the genes = ";
    std::cout << dp[chainLength][chainLength] << "\n";
    std::cout << "The aligned genes are :\n";
    for (i = id; i <= l; i++)
    {
        std::cout << (char)xans[i];
    }
    std::cout << "\n";
    for (i = id; i <= l; i++)
    {
        std::cout << (char)yans[i];
    }

    int rezultat = dp[chainLength][chainLength];

    for (int i = 0; i < 2 * chainLength; ++i)
        delete[] dp[i];
    delete[] dp;
    delete[] xans;
    delete[] yans;
    return rezultat;
}

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>
#include <string>

using namespace std;
using DataFrame = vector<string>;

double square(double value) {
    return value * value;
}

double distance(string first, string second) {
    //return square(first.x - second.x) + square(first.y - second.y);
    //GLOBAL SEQUENCE ALLIGNMENT
}

DataFrame k_means(const DataFrame& data,
    size_t k,
    size_t number_of_iterations) {
    static std::random_device seed;
    static std::mt19937 random_number_generator(seed());
    std::uniform_int_distribution<size_t> indices(0, data.size() - 1);

    // Pick centroids as random points from the dataset.
    DataFrame means(k);
    for (auto& cluster : means) {
        cluster = data[indices(random_number_generator)];
    }

    std::vector<size_t> assignments(data.size());
    for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
        // Find assignments.
        for (size_t chain = 0; chain < data.size(); ++chain) {
            double best_distance = std::numeric_limits<double>::max();
            size_t best_cluster = 0;
            for (size_t cluster = 0; cluster < k; ++cluster) {
                const double distance =
                    squared_l2_distance(data[chain], means[cluster]);
                if (distance < best_distance) {
                    best_distance = distance;
                    best_cluster = cluster;
                }
            }
            assignments[chain] = best_cluster;
        }

        // Divide sums by counts to get new centroids.
        //TU SPOA
        for (size_t cluster = 0; cluster < k; ++cluster) {
            // Turn 0/0 into 0/1 to avoid zero division.
            const auto count = std::max<size_t>(1, counts[cluster]);
            means[cluster].x = new_means[cluster].x / count;
            means[cluster].y = new_means[cluster].y / count;
        }
    }

    return means;
}


