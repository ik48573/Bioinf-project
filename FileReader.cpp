#include <string>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>
#include <map>
#include "spoa/include/spoa/spoa.hpp"

using namespace std;
namespace fs = std::filesystem;
using DataFrame = vector<string>;

int readFileNames(std::string* path);
int findMostCommonLength(std::vector<unsigned int>* allLengthsVector);
std::vector<std::string> collectChains(std::string path, std::string fileName);
int findMinimumDistance(std::string chain1, std::string chain2, int chain1Length, int chain2Length);
DataFrame k_means(const DataFrame& data, size_t k, size_t number_of_iterations);

std::vector<std::string> fileNames;

int main(int argc, char *argv[])
{
    int lengthValue = 0;
    std::string path = "fastq";
  
    std::string gene1 = "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCTGGCGAGGTGACGCGAA";
    std::string gene2 = "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA";
    int dist = findMinimumDistance(gene1, gene2, gene1.length(), gene2.length());
    cout << "udaljenost: " <<dist << endl;
    cout << argc << endl;

    if (argc == 1) {
        readFileNames(&path);
        //lengthValue = findMostCommonLength(&allLengthsVector);
    //readFile(&path, true, &lengthValue);

        for (int i = 0; i < fileNames.size(); i++) {
            //std::cout << fileNames.at(i) << std::endl;
            std::vector<std::string> chainsFromFile = collectChains(path, fileNames.at(i));
            std::cout << "Dohvaćeni lanci" << std::endl;
            vector<string> consensus = k_means(chainsFromFile, 4, 1);
            cout << consensus.at(0);
        }
    }
    else if (argc == 2) {
        string fileName = argv[1];
        vector<string> chainsFromFile = collectChains(path, fileName);
        cout << "Dohvaćeni lanci" << endl;

        for (int i = 0; i < chainsFromFile.size(); i++) {
            //cout << chainsFromFile.at(i) << endl;
        }
        vector<string> consensus = k_means(chainsFromFile, 4, 1);
        for (int i = 0; i < consensus.size(); i++) {
            cout << i << " klaster: " << consensus.at(i) << endl;
        }
    }

    
     
    return 0;
}

int readFileNames(std::string *path) {
    for (const auto& entry : fs::directory_iterator(*path)) {
        std::string fileName = fs::path(entry.path()).filename().string();
        if (fileName.rfind("J", 0) == 0) {
            if (std::find(fileNames.begin(), fileNames.end(), fileName) != fileNames.end()) {}
            else {
                fileNames.push_back(fileName);
            }
        }
    }
    return 0;
}

int findMostCommonLength(std::vector<unsigned int> *allLengthsVector) {
    std::vector <unsigned int> lengthVector; //izvojene duljine lanaca 

    int maxRepetition = 0;
    int lengthValueIndex;   //na kojem mjestu unutar lengthVector-a se nalazi duljina s najvise ponavljanja 
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
            lengthValueIndex = lengthVector.at(i);
        }
    }
    std::cout << "Najveci broj pojavljivanja je za: " << lengthValueIndex << std::endl;
    return lengthValueIndex;
}

std::vector<std::string> collectChains(std::string path, std::string fileName) {

    std::string fullPath = path + "/" + fileName;
    std::vector<unsigned int> allLengthsVector;
    std::vector<std::string> validChains;

    std::ifstream fileOpen(fullPath);
    if (!fileOpen.is_open()) {
        std::perror("File opening failed");
    }
    else
    {
        std::cout << "Reading: " << fullPath << std::endl;
        std::string line;
        bool flag = false;
        while (fileOpen.good())
        {
            getline(fileOpen, line);
            if (flag == true) {
                unsigned int lineLength = line.length();
                allLengthsVector.push_back(line.length());
                flag = false;
            }
            if (line.rfind("@16WBS", 0) == 0) {
                flag = true;        //sljedeća linija je lanac
            }
        }
        fileOpen.close();
        std::cout << allLengthsVector.size() << std::endl;
    }

    int lengthValue = findMostCommonLength(&allLengthsVector);
    std::cout << "Duljina je: " << lengthValue << std::endl;

    std::ifstream fileOpen2(fullPath);
    if (!fileOpen2.is_open()) {
        std::perror("File opening failed");
    }
    else
    {
        std::cout << "Reading: " << fullPath << std::endl;
        std::string line;
        bool flag = false;
        while (fileOpen2.good())
        {
            getline(fileOpen2, line);
            if (flag == true) {
                unsigned int lineLength = line.length();
                //allLengthsVector.push_back(line.length());
                if (line.length() == lengthValue) {
                    //cout << line << endl;
                    validChains.push_back(line);
                }
                flag = false;
            }
            if (line.rfind("@16WBS", 0) == 0) {
                flag = true;        //sljedeća linija je lanac
            }
        }
        fileOpen2.close();
    }
    return validChains;
}

int findMinimumDistance(std::string chain1, std::string chain2, int chain1Length, int chain2Length) {
    int misMatchPenalty = 7;
    int gapPenalty = 8;

    int i, j; // intialising variables  

    const int n = chain1Length;
    const int m = chain2Length;

    // table for storing optimal substructure answers 
    int** dp = new int* [n+m+1];
    for (int i = 0; i < (n+m+1); ++i)
        dp[i] = new int[n+m+1];

    for (int i = 0; i < (n + m + 1); i++) {
        for (int j = 0; j < (n + m + 1); j++) {
            dp[i][j] = 0;
        }
    }
    //std::cout << "kraj";
    // intialising the table 
    for (i = 0; i <= (n + m); i++)
    {
        dp[i][0] = i * gapPenalty;
        dp[0][i] = i * gapPenalty;
    }

    // calcuting the minimum penalty 
    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= m; j++)
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
    int l = n+m; // maximum possible length 

    i = n; j = m;

    int xpos = l;
    int ypos = l;

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
    /*std::cout << "Minimum Penalty in aligning the genes = ";
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
    }*/

    int rezultat = dp[m][n];

    for (int i = 0; i < (n+m+1); ++i)
        delete[] dp[i];
    delete[] dp;
    delete[] xans;
    delete[] yans;
    return rezultat;
}


double square(double value) {
    return value * value;
}


DataFrame k_means(const DataFrame& data,
    size_t k,
    size_t number_of_iterations) {
    static std::random_device seed;
    static std::mt19937 random_number_generator(seed());
    std::uniform_int_distribution<size_t> indices(0, data.size() - 1);

    //map<int, vector<string> > clusterChainMap;
    vector<vector<string>> clusterChainMap(k);
    for (int i = 0; i < k; i++) {
        clusterChainMap.at(i).push_back("");
    }
    cout << "odaberi random centroid" << endl;
    // Pick centroids as random points from the dataset.
    DataFrame means(k);
    for (auto& cluster : means) {
        cluster = data[indices(random_number_generator)];
        cout << cluster << endl;
    }
    cout << "Odabran centroid\n";

    std::vector<size_t> assignments(data.size());
    for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
        // Find assignments.
        cout << "Pocetak iteracije " + (iteration+1)<< endl;
        for (size_t chain = 0; chain < data.size(); ++chain) {
            int best_distance = std::numeric_limits<int>::max();
            size_t best_cluster = 7;
            for (size_t cluster = 0; cluster < k; ++cluster) {
                cout << "trazi klaster" << endl;
                cout << "racunaj distance" << endl;
                //const int distance =
                    //distanceBetweenTwoSequences(data[chain], means[cluster]);
                const int distance = findMinimumDistance(data[chain], means[cluster], data[chain].length(), means[cluster].length());
                cout << data[chain]<< endl;
                cout << means[cluster] << endl;
                cout << "Izracunao distance" << endl;
                cout << "Najbolja distanca: " << best_distance << " i distanca: " << distance<<endl;
                cout << "Najbolji klaster: " << best_cluster << "i trenutni klaster: " << cluster<<endl;
                if (distance <= best_distance) {
                    best_distance = distance;
                    best_cluster = cluster;
                }
                cout << "Nasao klaster" << endl;
            }
            //clusterChainMap[best_cluster].push_back(data[chain]);
            cout << "Spremi lanac u klaster" << endl;
            clusterChainMap.at(best_cluster).push_back(data[chain]);
            cout << "Spremio lanac u klaster" << endl;
            cout << "zavrsen lanac: " << chain<< "/" <<data.size()<< endl;

        }


        // Divide sums by counts to get new centroids.
        //TU SPOA


        auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(1),
            5, -4, -8, -6);

        auto graph = spoa::createGraph();
        for (int i = 0; i < k; i++) {
            for (const auto& it : clusterChainMap[i]) {
                auto alignment = alignment_engine->align(it, graph);
                graph->add_alignment(alignment, it);
            }
            std::string consensus = graph->generate_consensus();
            means[i] = consensus.c_str();
        }

       /* std::string consensus = graph->generate_consensus();

        fprintf(stderr, "Consensus (%zu)\n", consensus.size());
        fprintf(stderr, "%s\n", consensus.c_str());
        for (size_t cluster = 0; cluster < k; ++cluster) {
            // Turn 0/0 into 0/1 to avoid zero division.
            const auto count = std::max<size_t>(1, counts[cluster]);
            means[cluster] = consensus.c_str();
        }*/
        cout << "kraj iteracije" << endl;
    }

    return means;
}


