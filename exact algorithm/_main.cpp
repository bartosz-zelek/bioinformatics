#include <iostream>
#include <map>
#include <vector>

using namespace std;

map<string, string> nucleotide_to_weak_strong = {
    {"A", "W"},
    {"T", "W"},
    {"C", "S"},
    {"G", "S"},
};

map<string, string> nucleotide_to_purine_pyrimidine = {
    {"A", "R"},
    {"T", "Y"},
    {"C", "Y"},
    {"G", "R"},
};

class WSRY
{
    // public:
private:
    string start_converted;
    map<string, string> map_convertion;
    map<string, bool> map_cells;
    vector<string> path;
    vector<int> depth;

    string convert_oligo(string oligo)
    {
        string converted_oligo = "";
        for (int i = 0; i < oligo.size() - 1; i++)
        {
            converted_oligo += 'TODO';
        }
        return converted_oligo;
    }
};