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
vector<string> merge_clusters(vector<string> first_cluster, vector<string> second_cluster);

vector<string> file_names;

int main(int argc, char *argv[])
{
    int deer_k = 4;
    int chamois_k = 6;
    int number_of_iterations = 5;
    string path = "fastq";
    fs::create_directory("fasta");

    if (argc == 1) {
        readFileNames(&path);

        for (int i = 0; i < file_names.size(); i++) {
            vector<string> chains_from_file = collectChains(path, file_names.at(i));
            vector<string> consensus = k_means(chains_from_file, deer_k, number_of_iterations);
        }
    }
    else if (argc == 2) {
        string file_name = argv[1];
        vector<string> chains_from_file = collectChains(path, file_name);
        vector<string> consensus;
        if (file_name.rfind("J", 0) == 0) {
            consensus = k_means(chains_from_file, deer_k, number_of_iterations);
        }
        else if (file_name.rfind("L", 0) == 0) {
            consensus = k_means(chains_from_file, chamois_k, number_of_iterations);
        }
        else {
            cout << "Please put valid file name!" << endl;
        }
        string path_save = "fasta/";
        string new_file_name = file_name.substr(0, file_name.size() - path.size()-1); //remove .fastq
        
        string new_file;
        time_t t = time(0);   // get time now
        struct tm* now = localtime(&t);

        char buffer[80];
        strftime(buffer, 80, "%d%m%Y_%H%M%S", now);
        new_file.append(path_save);
        new_file.append(new_file_name);
        new_file.append("_");
        new_file.append(buffer);
        new_file.append(".fasta");
        cout << new_file << endl;
        ofstream file_save;
        file_save.open(new_file);
        for (int i = 0; i < consensus.size(); i++) {
            if (consensus.at(i) != "") {
                file_save << ">Varijanta gena: " << (i + 1) << endl;
                file_save << consensus.at(i) << endl;
                cout << consensus.at(i) << endl;
            }
        }
        file_save.close();
    }
    return 0;
}

int readFileNames(string *path) {
    for (const auto& entry : fs::directory_iterator(*path)) {
        string file_name = fs::path(entry.path()).filename().string();
        if (file_name.rfind("J", 0) == 0) {
            if (find(file_names.begin(), file_names.end(), file_name) != file_names.end()) {}
            else {
                file_names.push_back(file_name);
            }
        }
    }
    return 0;
}

int findMostCommonLength(vector<unsigned int> *allLengthsVector) {
    vector <unsigned int> length_vector; //izvojene duljine lanaca 

    int max_repetition = 0;
    int length_value_index = 0;   //na kojem mjestu unutar lengthVector-a se nalazi duljina s najvise ponavljanja 
    int flag;

    for (int i = 0; i < (*allLengthsVector).size();i++) {   //izdvoji duljine lanaca bez ponavljanja
        flag = 1;
        for (int j = 0; j < length_vector.size();j++) {
            if (length_vector.at(j) == (*allLengthsVector).at(i))
                flag = 0;
        }
        if (flag == 1) {
            length_vector.push_back((*allLengthsVector).at(i));
        }
    }

    //cout << "Ukupno ima: " <<lengthVector.size() <<" razlicitih duljina lanaca"<< endl;

    //prebroji koliko ponavljanja ima za svaku duljinu
    for (int i = 0; i < length_vector.size();i++) {      
        int my_count = count((*allLengthsVector).begin(), (*allLengthsVector).end(), length_vector.at(i));
        if (my_count > max_repetition) {     //odredi maksimum
            max_repetition = my_count;
            length_value_index = length_vector.at(i);
        }
    }
    return length_value_index;
}

vector<string> collectChains(string path, string fileName) {

    string full_path = path + "/" + fileName;
    vector<unsigned int> all_lengths_vector;
    vector<string> valid_chains;

    ifstream fileOpen(full_path);
    if (!fileOpen.is_open()) {
        perror("File opening failed");
    }
    else
    {
        cout << "Reading: " << full_path << endl;
        string line;
        bool flag = false;
        while (fileOpen.good())
        {
            getline(fileOpen, line);
            if (flag == true) {
                unsigned int lineLength = line.length();
                all_lengths_vector.push_back(line.length());
                flag = false;
            }
            if (line.rfind("@16WBS", 0) == 0) {
                flag = true;        //sljedeća linija je lanac
            }
        }
        fileOpen.close();
    }

    int length_value = findMostCommonLength(&all_lengths_vector); //najčešća duljina prisutna u datoteci

    ifstream fileOpen2(full_path);
    if (!fileOpen2.is_open()) {
        perror("File opening failed");
    }
    else
    {
        string line;
        bool flag = false;
        while (fileOpen2.good())
        {
            getline(fileOpen2, line);
            if (flag == true) {
                unsigned int line_length = line.length();
                if (line.length()== length_value) {
                    valid_chains.push_back(line);
                }
                flag = false;
            }
            if (line.rfind("@16WBS", 0) == 0) {
                flag = true;        //sljedeća linija je lanac
            }
        }
        fileOpen2.close();
    }
    return valid_chains;
}

int findMinimumDistance(string chain1, string chain2, int chain1Length, int chain2Length) {
    int mis_match_penalty = 7;
    int gap_penalty = 8;

    int i, j; // Intialising variables.  

    const int n = chain1Length;
    const int m = chain2Length;

    // Table for storing optimal substructure answers. 
    int** dp = new int* [n+m+1];
    for (int i = 0; i < (n+m+1); ++i)
        dp[i] = new int[n+m+1];

    for (int i = 0; i < (n + m + 1); i++) {
        for (int j = 0; j < (n + m + 1); j++) {
            dp[i][j] = 0;
        }
    }

    for (i = 0; i <= (n + m); i++)
    {
        dp[i][0] = i * gap_penalty;
        dp[0][i] = i * gap_penalty;
    }

    // Calcuting the minimum penalty. 
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
                dp[i][j] = min({ dp[i - 1][j - 1] + mis_match_penalty ,
                                dp[i - 1][j] + gap_penalty ,
                                dp[i][j - 1] + gap_penalty });
            }
        }
    }

    int result = dp[n][m];

    for (int i = 0; i < (n+m+1); ++i)
        delete[] dp[i];
    delete[] dp;
    return result;
}

vector<string> k_means(const vector<string>& data,
    int k,
    int number_of_iterations) {

    vector<vector<string>> cluster_chain_map(k);
    int cluster_merge_threshold = 20;

    // Pick centroids as random points from the dataset.
    vector<string> means;
    vector<string> check_means;

    for (int i = 0; i < k; ++i) {
        cluster_chain_map.at(i).push_back("");
        means.push_back(data[rand() % data.size()]);
        check_means.push_back("");
    }

    for (int iteration = 0; iteration < number_of_iterations; ++iteration) {

        // Find assignments.
        for (int chain = 0; chain < data.size(); ++chain) {
            int best_distance = numeric_limits<int>::max();
            int best_cluster = 0;
            for (int cluster = 0; cluster < k; ++cluster) {
                if (iteration > 0) {
                    check_means[cluster] = means[cluster];
                }
                const int distance = findMinimumDistance(data[chain], means[cluster], data[chain].length(), means[cluster].length());
                if (distance < best_distance) {
                    best_distance = distance;
                    best_cluster = cluster;
                }
            }
            cluster_chain_map.at(best_cluster).push_back(data[chain]);

        }
 
        for (int i = 0; i < k; ++i) {
            auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(1), 0, -4, -8, -6);
            auto graph = spoa::createGraph();
            if (cluster_chain_map[i].size() > 1) {
                for (const auto& it : cluster_chain_map[i]) {
                    auto alignment = alignment_engine->align(it, graph);
                    graph->add_alignment(alignment, it);
                }
                string consensus = graph->generate_consensus();
                means[i] = consensus;
            }
        }

        int flag = 0;
        if (iteration != 0) {
            cout << "Check last session." << endl;
            for (int i = 0; i < means.size(); i++) {
                if (check_means[i] == means[i]) {
                    flag++;
                }
            }
            if (flag == means.size()) {
                cout << "Found stationary state." << endl;
                break;
            }
        }
    }
    cout << "Merging threshold: " << cluster_merge_threshold << "." << endl;
    while (1) {
        int a = -1;
        int b = -1;
        int flag = -1;

        for (int i = k - 1; i > 0; --i) {
            flag = 1;
            for (int j = 0; j < i; ++j) {
                a = i;
                b = j;
                if (means[i] != "" && means[j] != "") {
                    int distance = findMinimumDistance(means[i], means[j], means[i].length(), means[j].length());
                    if (distance < cluster_merge_threshold) {
                        cout << "Merging clusters " << i << " and " << j << " because distance between them is: " << distance << endl;
                        means[i] = "";
                        cluster_chain_map[j] = merge_clusters(cluster_chain_map[i], cluster_chain_map[j]);
                        cluster_chain_map[i].clear();
                        auto alignment_engine = spoa::createAlignmentEngine(static_cast<spoa::AlignmentType>(1), 0, -4, -8, -6);
                        auto graph = spoa::createGraph();
                        for (const auto& it : cluster_chain_map[j]) {
                            auto alignment = alignment_engine->align(it, graph);
                            graph->add_alignment(alignment, it);
                        }
                        means[j] = graph->generate_consensus();
                        flag = 0;
                        break;
                    }
                }
            }
            if (flag == 0) {
                break;
            }
        }
        if (a == 1 && b == 0) {
            break;
        }
        if (cluster_chain_map.size() == 1 || flag == 1) {
            break;
        }
    }
    return means;
}

vector<string> merge_clusters(vector<string> first_cluster, vector<string> second_cluster) {
    if (first_cluster.size() >= second_cluster.size()) {
        for (int i = 0; i < second_cluster.size(); ++i) {
            first_cluster.push_back(second_cluster[i]);
        }
        return first_cluster;
    }
    else
    {
        for (int i = 0; i < first_cluster.size(); ++i) {
            second_cluster.push_back(first_cluster[i]);
        }
        return second_cluster;
    }
 }