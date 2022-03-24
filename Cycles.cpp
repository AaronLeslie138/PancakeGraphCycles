
#include <iostream>
#include <unordered_map>
#include <string>
#include <vector> 
#include <functional>
#include <unordered_set>
#include <algorithm> 
#include <chrono>
#include <math.h>   
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/graph/hawick_circuits.hpp>
#include <boost/graph/visitors.hpp>
#include <deque>
#include <fstream>

std::string get_child(std::string vertex, int instance, int side) {
    std::string original_name = vertex;
    std::string end_name = "";
    int count = 0;
    for (int j = 0; j < instance; j++) {
        for (int k = 1; k < original_name.length(); k++) {
            if (original_name.at(k) == 'b') {
                count = count + 1;
            }
            else {
                break;
            }
        }
        count = count + 1;
        std::string temp_name = original_name.substr(0, count);
        original_name = original_name.substr(count, INT_MAX);
        if (temp_name.length() == side) {
            temp_name = temp_name[0];
        }
        else {
            temp_name = temp_name + 'b';
        }
        end_name = temp_name + end_name;
        count = 0;
    }
    end_name = end_name + original_name;
    return end_name;
}


void generate_graph(std::string vertex, std::unordered_map<std::string, std::unordered_set<std::string>> &gdict, int s, int p) {
    if (gdict.find(vertex) != gdict.end()) {
        return;
    }
    int instances = 1;
    std::unordered_set<std::string> children;
    for (int flip = 0; flip < p; flip++) {
        children.insert(get_child(vertex, instances, s));
        instances = instances + 1;
        gdict[vertex] = children;
        for(std::string vrt : gdict[vertex]){
            generate_graph(vrt, gdict, s, p);
        }
    }

}


static std::unordered_map<std::string, int> pairings;

int stringToInt(std::string str) {
    static int biggestNum = 0;
    if (pairings.find(str) != pairings.end()) {
        return pairings[str];
    }
    pairings[str] = biggestNum++;
    return pairings[str];
}

std::string intToString(int k) {


    const int prevToFind = k;
    auto findResult = std::find_if(std::begin(pairings), std::end(pairings), [&](const std::pair<std::string, int>& pair)
    {
        return pair.second == prevToFind;
    });


    if (findResult != std::end(pairings))
    {
        return findResult->first;
    }

    return "oops";
}



#include <string>
#include <sys/stat.h> // stat
#include <errno.h>    // errno, ENOENT, EEXIST
#if defined(_WIN32)
#include <direct.h>   // _mkdir
#endif

bool isDirExist(const std::string& path)
{
#if defined(_WIN32)
    struct _stat info;
    if (_stat(path.c_str(), &info) != 0)
    {
        return false;
    }
    return (info.st_mode & _S_IFDIR) != 0;
#else 
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
    {
        return false;
    }
    return (info.st_mode & S_IFDIR) != 0;
#endif
}

bool makePath(const std::string& path)
{
#if defined(_WIN32)
    int ret = _mkdir(path.c_str());
#else
    mode_t mode = 0755;
    int ret = mkdir(path.c_str(), mode);
#endif
    if (ret == 0)
        return true;

    switch (errno)
    {
    case ENOENT:
        // parent didn't exist, try to create it
    {
        int pos = path.find_last_of('/');
        if (pos == std::string::npos)
#if defined(_WIN32)
            pos = path.find_last_of('\\');
        if (pos == std::string::npos)
#endif
            return false;
        if (!makePath(path.substr(0, pos)))
            return false;
    }
    // now, try to create again
#if defined(_WIN32)
    return 0 == _mkdir(path.c_str());
#else 
    return 0 == mkdir(path.c_str(), mode);
#endif

    case EEXIST:
        // done!
        return isDirExist(path);

    default:
        return false;
    }
}

int fact(int n) {

    return (n == 0) || (n == 1) ? 1 : n * fact(n - 1);
}

typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS > Graph;
typedef boost::reverse_graph<Graph> Rgraph;
typedef Graph::vertex_descriptor Vertex;
typedef std::vector<Vertex> Stack;

class{
    public:
        
        void cycle(const Stack& s, const Graph& g) {
            if (s.size() == 2)
                return;
            static bool lengths[400];
            static std::string examples[400];
            if (s.size() == 0) {
                std::string output = "";
                for (int i = 0; i < 400; i++) {
                    if (lengths[i]) output += std::to_string(i) + ",";
                }
                std::cout << output << "\n";
                output = "";
                for (int i = 0; i < 400; i++) {
                    if(lengths[i]) output += "[" + std::to_string(i) + "]: " + examples[i] + " ";
                }
                std::cout << output << "\n";
                std::ofstream outfile("finished.txt");
                outfile << "done!" << std::endl;
                outfile.close();
                return;
            }
            if (lengths[s.size()]) return;

            std::string example = "";
            for (int v : s) {
                example += intToString(v) + "-";
            }
            examples[s.size()] = example;
            std::ofstream outfile(std::to_string(s.size()) + ".txt", std::ios_base::trunc);
            outfile << example << std::endl;
            outfile.close();

            lengths[s.size()] = true;
            std::string output = "";
            for (int i = 0; i < 400; i++) {
                if (lengths[i]) output += std::to_string(i) + ",";
            }
            std::cout << output << "\n";
        }
        void print_lengths() {
            //for (int length : lengths) std::cout << length + ",";
        }
    
} visitor;

#ifdef _WIN32
#include <direct.h>
// MSDN recommends against using getcwd & chdir names
#define cwd _getcwd
#define cd _chdir
#else
#include "unistd.h"
#define cwd getcwd
#define cd chdir
#endif

int main()
{
    using namespace boost;
    std::cout << "Enter the number of Sides per Permutation: ";
    int sides = 0;
    std::cin >> sides;
    std::cout << "\n";

    std::cout << "Enter the number of  Permutation: ";
    int pancakes = 0;
    std::cin >> pancakes;
    std::cout << "\n";

    makePath("" + std::to_string(sides) + "-" + std::to_string(pancakes) + "-U");
    cd((std::to_string(sides) + "-" + std::to_string(pancakes) + "-U").c_str());

    auto start = std::chrono::high_resolution_clock::now();

    std::string root = "";
    for (int pancake = 1; pancake < pancakes + 1; pancake++) {
        root = root + std::to_string(pancake);
    }

    std::unordered_map<std::string, std::unordered_set<std::string>> graphdict;
    generate_graph(root, graphdict, sides, pancakes);
   
    
    Graph graph;
    for (const auto& iter : graphdict) {
        for (const auto& end : iter.second) {
            boost::add_edge(stringToInt(iter.first), stringToInt(end), graph);
        }
    }
    boost::hawick_unique_circuits(graph, visitor);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Elapsed time: " << (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) << " ms";
    visitor.print_lengths();
}
