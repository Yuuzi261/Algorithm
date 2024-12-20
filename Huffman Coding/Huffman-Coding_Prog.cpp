#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <string>
#include <sstream>
#include <bitset>
#include <unordered_map>
#include <algorithm>
#include <filesystem>

using namespace std;

class Node {
public:
    string val;
    int freq;
    Node* lchild;
    Node* rchild;

    Node(string _val, int _freq, Node* _lchild, Node* _rchild) {
        val = _val;
        freq = _freq;
        lchild = _lchild;
        rchild = _rchild;
    }
};

struct cmp {
    bool operator() (Node* a, Node* b) {
        return a->freq > b->freq;
    }
};

unordered_map<string, int> freq;
unordered_map<string, string> HuffmanMap;    
priority_queue<Node*, vector<Node*>, cmp> forest;

void ReadPGMFile(string);
void BuildHuffmanTree();
void Compression(vector<string>);
void Decompression(string, string);
void CompressionTree(Node*, string&);
Node* DecompressionTree(istringstream&);
void MakeMap(Node*, string);
void OutputCSV();

int main(int argc, char *argv[]) {

    if(argc < 3) {
        cout << "Error mode or input" << endl;
        return 0;
    }
    else {
        string mode = string(argv[1]);
        if(mode == "-c") {
            vector<string> input_files;
            for(int i = 2;i < argc;i++) {
                input_files.push_back(string(argv[i]));
                cout << "analyzing " << filesystem::path(input_files.back()).filename().string() << "..." << endl;
                ReadPGMFile(input_files.back());
            }
            cout << endl << "Frequencies:" << endl << endl;
            OutputCSV();
            BuildHuffmanTree();
            vector<pair<string, int>> sorted_freq;
            for (const auto &item : freq) {
                sorted_freq.emplace_back(item);
            }
            sort(sorted_freq.begin(), sorted_freq.end(), [](const auto &x, const auto &y) { return x.second > y.second; });
            for(const auto &[key, value] : sorted_freq) cout << key << " (freq: " << value << "): " << HuffmanMap[key] << endl;
            cout << endl << "Files size:" << endl << endl;
            Compression(input_files);
        }
        else if(mode == "-d") {
            Decompression(string(argv[2]), string(argv[3]));
            cout << "Completed" << endl;
        }
        else cout << "Error mode" << endl;
    }

    return 0;
}

void ReadPGMFile(string input_file) {

    ifstream ifs;

    ifs.open(input_file);

    string line, pixel;
    freq["_s_"] = 0;
    while(getline(ifs, line)) {
        istringstream iss(line);
        while(iss >> pixel) {
            if(freq.count(pixel)) freq[pixel]++;
            else freq[pixel] = 1;
            freq["_s_"]++;
        }
        freq["\\n"]++;
    }
    ifs.close();

}

void BuildHuffmanTree() {

    for(auto it = freq.begin(); it != freq.end(); it++)
        forest.push(new Node((*it).first, (*it).second, nullptr, nullptr));

    for(int i = 0; i < freq.size()-1; i++) {
        Node* ptr1 = forest.top(); forest.pop();
        Node* ptr2 = forest.top(); forest.pop();
        Node* parentNode = new Node("`", ptr1->freq + ptr2->freq, ptr1, ptr2);
        forest.push(parentNode);
    }
    MakeMap(forest.top(), string(""));

}

void Compression(vector<string> input_files) {

    ifstream ifs;
    ofstream ofs;

    string header = "";
    CompressionTree(forest.top(), header);

    ofs.open("output.hc", std::ios::binary);
    int length = header.length();
    ofs.write((char*)&(length), sizeof(int));
    for (int i = 0; i < length; i++) 
        ofs.write(&header[i], 1);

    string content = "";
    vector<int> content_info;
    for(string input_file : input_files) {

        cout << filesystem::path(input_file).filename().string() << ": " << filesystem::file_size(input_file) / 1024.0 << " KB" << endl;
        ifs.open(input_file);

        string tmp = "";
        string line, pixel;
        while(getline(ifs, line)) {
            istringstream iss(line);
            while(iss >> pixel) {
                tmp += HuffmanMap[pixel];
                tmp += HuffmanMap["_s_"];
            }
            tmp += HuffmanMap["\\n"];
        }
        ifs.close();
        for(;tmp.length() % 8 != 0;tmp += "0");
        content += tmp;
        content_info.push_back(tmp.length() / 8);

    }

    int content_nums = content_info.size();
    ofs.write((char*)&content_nums, sizeof(int));
    for(int i = 0;i < content_nums;i++)
        ofs.write((char*)&content_info[i], sizeof(int));

    for(int i = 0;i < content.length();i+=8) {
        bitset<8> btmp(content.substr(i, 8));
        ofs.write((char*)&btmp, 1);
    }
    ofs.close();
    cout << "compressed file: " << filesystem::file_size("output.hc") / 1024.0 << " KB" << endl;

}

void Decompression(string zip_file, string output_file) {

    int tree_data_len, file_nums;
    vector<int> files_len;
    string tree_data = "";
    vector<string> huffman_codes;
    ifstream ifs;
    ofstream ofs;

    ifs.open(zip_file, std::ios::binary);
    ifs.read((char*)&tree_data_len, sizeof(int));

    char ctmp;
    for(int i=0;i<tree_data_len;i++) {
        ifs.read((char*)&ctmp, 1);
        tree_data += ctmp;
    }

    istringstream iss(tree_data);
    Node* root = DecompressionTree(iss);

    int itmp;
    ifs.read((char*)&file_nums, sizeof(int));
    for(int i = 0;i < file_nums;i++) {
        ifs.read((char*)&itmp, sizeof(int));
        files_len.push_back(itmp);
    }

    bitset<8> btmp;
    for(int i = 0;i < file_nums;i++) {
        string stmp = "";
        for(int j = 0;j < files_len[i];j++) {
            ifs.read((char*)&btmp, 1);
            stmp += btmp.to_string();
        }
        huffman_codes.push_back(stmp);
    }

    ifs.close();

    for(int i = 0;i < huffman_codes.size();i++) {
        Node* now = root;
        string output = "";
        cout << "unzip file " << i+1 << "..." << endl;
        for(char dir : huffman_codes[i]) {
            if(dir == '0') now = now->lchild;
            else now = now->rchild;

            if(now->lchild == nullptr && now->rchild == nullptr) {
                if(now->val == "_s_") output += " ";
                else if(now->val == "\\n") output += "\n";
                else output += now->val;
                now = root;
            }
        }

        string ofile_name = output_file.substr(0, output_file.find(".pgm"));
        ofs.open(ofile_name + "_" + to_string(i+1) + ".pgm", std::ios::binary);
        ofs << output;
        ofs.close();

    }

}

void CompressionTree(Node* ptr, string& output) { 
    if(ptr == nullptr) { output += "_n_ "; return; }

    output += ptr->val + ' ';
    CompressionTree(ptr->lchild, output);
    CompressionTree(ptr->rchild, output);
}

Node* DecompressionTree(istringstream& iss) {
    string next;
    iss >> next;
    if (next == "_n_") {
        return nullptr;
    }
    else {
        Node* node = new Node(next, 0, nullptr, nullptr);
        node->lchild = DecompressionTree(iss);
        node->rchild = DecompressionTree(iss);
        return node;
    }
}

void MakeMap(Node* ptr, string s) {
    if(ptr->lchild == nullptr || ptr->rchild == nullptr) {
        HuffmanMap[ptr->val] = s;
        return;
    }
    if(ptr->lchild) MakeMap(ptr->lchild, string(s + "0"));
    if(ptr->rchild) MakeMap(ptr->rchild, string(s + "1"));
}

void OutputCSV() {
    ofstream ofs;

    ofs.open("freq.csv");
    for(int i = 0;i <= 255;i++) 
        ofs << i << ", " << freq[to_string(i)] << endl;

    ofs.close();
}