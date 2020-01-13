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
#include "include/spoa/spoa.hpp"
using namespace std;
namespace fs = filesystem;

int readFileNames(string* path);
int findMostCommonLength(vector<unsigned int>* allLengthsVector);
vector<string> collectChains(string path, string fileName);
int findMinimumDistance(string chain1, string chain2, int chain1Length, int chain2Length);
vector<string> k_means(const vector<string>& data, int k, int number_of_iterations);

vector<string> fileNames;

int main(int argc, char *argv[])
{
    //int lengthValue = 0;
    string path = "fastq";
  
    /*string gene1 = "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCTGGCGAGGTGACGCGAA";
    string gene2 = "GATCCTCTCTCTGCAGCACATTTCCTGCTGTATGCTAAGAGCGAGTGTCATTTCTCCAACGGGACGCAGCGGGTGGGGTTCCTGGACAGATACTTCTATAACGGAGAAGAGTTCGTGCGCTTCGACAGCGACTGGGGCGAGTACCGGGCGGTGACAGAGCTGGGGCGGCCGGTGGCCGAGTACCTGAACAGCCAGAAGGAGTACATGGAGCAGACGCGGGCCGAGGTGGACACGTACTGCAGACACAACTACGGCGGCGTTGAGAGTTTCACTGTGCAGCGGCGAGGTGACGCGAA";
    int dist = findMinimumDistance(gene1, gene2, gene1.length(), gene2.length());*/

    if (argc == 1) {
        readFileNames(&path);
        //lengthValue = findMostCommonLength(&allLengthsVector);
    //readFile(&path, true, &lengthValue);

        for (int i = 0; i < fileNames.size(); i++) {
            //cout << fileNames.at(i) << endl;
            vector<string> chainsFromFile = collectChains(path, fileNames.at(i));
            cout << "Dohvaćeni lanci" << endl;
            vector<string> consensus = k_means(chainsFromFile, 4, 1);
        }
    }
    else if (argc == 2) {
        string fileName = argv[1];
        vector<string> chainsFromFile = collectChains(path, fileName);

        for (int i = 0; i < chainsFromFile.size(); i++) {
            //cout << chainsFromFile.at(i) << endl;
        }
        vector<string> consensus = k_means(chainsFromFile, 4, 1);
        for (int i = 0; i < consensus.size(); i++) {
            cout << consensus.at(i) << endl;
        }
    }
    return 0;
}

int readFileNames(string *path) {
    for (const auto& entry : fs::directory_iterator(*path)) {
        string fileName = fs::path(entry.path()).filename().string();
        if (fileName.rfind("J", 0) == 0) {
            if (find(fileNames.begin(), fileNames.end(), fileName) != fileNames.end()) {}
            else {
                fileNames.push_back(fileName);
            }
        }
    }
    return 0;
}

int findMostCommonLength(vector<unsigned int> *allLengthsVector) {
    vector <unsigned int> lengthVector; //izvojene duljine lanaca 

    int maxRepetition = 0;
    int lengthValueIndex=0;   //na kojem mjestu unutar lengthVector-a se nalazi duljina s najvise ponavljanja 
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
    cout << "Ukupno ima: " <<lengthVector.size() <<" razlicitih duljina lanaca"<< endl;

    //prebroji koliko ponavljanja ima za svaku duljinu
    for (int i = 0; i < lengthVector.size();i++) {      
        int mycount = count((*allLengthsVector).begin(), (*allLengthsVector).end(), lengthVector.at(i));
        if (mycount > maxRepetition) {     //odredi maksimum
            maxRepetition = mycount;
            lengthValueIndex = lengthVector.at(i);
        }
    }
    cout << "Najveci broj pojavljivanja je za: " << lengthValueIndex << endl;
    return lengthValueIndex;
}

vector<string> collectChains(string path, string fileName) {

    string fullPath = path + "/" + fileName;
    vector<unsigned int> allLengthsVector;
    vector<string> validChains;

    ifstream fileOpen(fullPath);
    if (!fileOpen.is_open()) {
        perror("File opening failed");
    }
    else
    {
        cout << "Reading: " << fullPath << endl;
        string line;
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
        cout << allLengthsVector.size() << endl;
    }

    int lengthValue = findMostCommonLength(&allLengthsVector);
    cout << "Duljina je: " << lengthValue << endl;

    ifstream fileOpen2(fullPath);
    if (!fileOpen2.is_open()) {
        perror("File opening failed");
    }
    else
    {
        cout << "Reading: " << fullPath << endl;
        string line;
        bool flag = false;
        while (fileOpen2.good())
        {
            getline(fileOpen2, line);
            if (flag == true) {
                unsigned int lineLength = line.length();
                //allLengthsVector.push_back(line.length());
                //if (line.length() <= lengthValue + 1 || line.length()>=lengthValue - 1) {
                if (line.length()== lengthValue) {

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

int findMinimumDistance(string chain1, string chain2, int chain1Length, int chain2Length) {
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
    //cout << "kraj";
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
                dp[i][j] = min({ dp[i - 1][j - 1] + misMatchPenalty ,
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


    int result = dp[m][n];

    for (int i = 0; i < (n+m+1); ++i)
        delete[] dp[i];
    delete[] dp;
    delete[] xans;
    delete[] yans;
    return result;
}

vector<string> k_means(const vector<string>& data,
    int k,
    int number_of_iterations) {

    //map<int, vector<string> > clusterChainMap;
    vector<vector<string>> clusterChainMap(k);


    // Pick centroids as random points from the dataset.
    vector<string> means;
    for (int i = 0; i < k; ++i) {
        clusterChainMap.at(i).push_back("");
        means.push_back(data[rand() % data.size()]);
    }

    vector<size_t> assignments(data.size());
    for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
        // Find assignments.
        for (size_t chain = 0; chain < data.size(); ++chain) {
            int best_distance = numeric_limits<int>::max();
            size_t best_cluster = 0;
            for (size_t cluster = 0; cluster < k; ++cluster) {
                //const int distance =
                    //distanceBetweenTwoSequences(data[chain], means[cluster]);
                const int distance = findMinimumDistance(data[chain], means[cluster], data[chain].length(), means[cluster].length());
                if (distance < best_distance) {
                    best_distance = distance;
                    best_cluster = cluster;
                }
            }
            //clusterChainMap[best_cluster].push_back(data[chain]);
            clusterChainMap.at(best_cluster).push_back(data[chain]);

        }
        //TU SPOA
        for (int i = 0; i < k; ++i) {
            auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(1),
                0, -4, -8, -6);
            auto graph = spoa::createGraph();
            for (const auto& it : clusterChainMap[i]) {
                auto alignment = alignment_engine->align(it, graph);
                graph->add_alignment(alignment, it);
            }

            //GREŠKA SEGMENTATION FAULT (CORE DUMPED) KOD POZIVA SLJEDEĆE LINIJE:
            string consensus = graph->generate_consensus();
            means[i] = consensus;
        }
        for (int i = 0; i < k; ++i) {
            cout << clusterChainMap.at(i).size() << endl;
        }
    }

    return means;
}
