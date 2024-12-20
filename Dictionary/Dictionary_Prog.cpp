#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <chrono>

using namespace std;

class Node {
    public:
        Node();
        Node(Node*, Node*);
        Node(Node*, Node*, Node*);
        Node(string);
        Node(Node*, Node*, string);
        Node(Node*, Node*, Node*, string);
        Node* getLeft();
        Node* getRight();
        Node* getParent();
        string getWord();
        void setLeft(Node*);
        void setRight(Node*);
        void setParent(Node*);
        void setWord(string);
        void setWord(Node*);
        void setWord(Node);
    private:
        string removePunctuationAndLower(string);
        Node *left, *right, *parent;
        string word;
};

class TreapNode : public Node {
    public:
        TreapNode();
        TreapNode(TreapNode*, TreapNode*);
        TreapNode(TreapNode*, TreapNode*, TreapNode*);
        TreapNode(string);
        TreapNode(TreapNode*, TreapNode*, string);
        TreapNode(TreapNode*, TreapNode*, TreapNode*, string);
        TreapNode* getLeft();
        TreapNode* getRight();
        TreapNode* getParent();
        int getPriority();
        void setLeft(TreapNode*);
        void setRight(TreapNode*);
        void setParent(TreapNode*);
        void setPriority(int);
    private:
        int priority;
};

class BinarySearchTree {
    public:
        BinarySearchTree();
        BinarySearchTree(Node*);
        void Insert(Node*);
        void Insert(Node);
        void Insert(string);
        Node* Search(Node*);
        Node* Search(Node);
        Node* Search(string);
        void Delete(Node*);
        void Delete(Node);
        void Delete(string);
    protected:
        Node* getRoot();
        void setRoot(Node*);
        Node* minValueNode(Node*);
    private:
        void insert(Node*, Node*);
        Node* search(Node*, string);
        Node* del(Node*, string);
        Node* root;
};

class SplayTree : public BinarySearchTree {
    public:
        SplayTree();
        SplayTree(Node*);
        void Insert(Node*);
        void Insert(Node);
        void Insert(string);
        Node* Search(Node*);
        Node* Search(Node);
        Node* Search(string);
        void Delete(Node*);
        void Delete(Node);
        void Delete(string);
    private:
        void splay(Node*);
        void leftRotation(Node*);
        void rightRotation(Node*);
        void del();
};

class Treap {
    public:
        Treap();
        Treap(TreapNode*);
        void Insert(TreapNode*);
        void Insert(TreapNode);
        void Insert(string);
        TreapNode* Search(TreapNode*);
        TreapNode* Search(TreapNode);
        TreapNode* Search(string);
        void Delete(TreapNode*);
        void Delete(TreapNode);
        void Delete(string);
    private:
        TreapNode* insert(TreapNode*, TreapNode*);
        TreapNode* search(TreapNode*, string);
        TreapNode* del(TreapNode*, string);
        TreapNode* rotateLeft(TreapNode*);
        TreapNode* rotateRight(TreapNode*);
        TreapNode* root;
};

void TestBinarySearchTree(string&, string&, vector<string>&);
void TestSplayTree(string&, string&, vector<string>&);
void TestTreap(string&, string&, vector<string>&);

int main(void) {

    srand(time(nullptr));

    string input_file1 = "TestFile1.txt", input_file2 = "TestFile2.txt";
    string delete_file1 = "TestFile14.txt", delete_file2 = "TestFile24.txt";
    vector<string> search_files1 = {"TestFile11.txt", "TestFile12.txt", "TestFile13.txt"}, search_files2 = {"TestFile21.txt", "TestFile22.txt", "TestFile23.txt"};

    // Binary Search Tree
    cout << endl << "--------Binary Search Tree--------" << endl;
    TestBinarySearchTree(input_file1, delete_file1, search_files1);
    TestBinarySearchTree(input_file2, delete_file2, search_files2);

    // Splay Tree
    cout << "\n--------Splay Tree--------" << endl;
    TestSplayTree(input_file1, delete_file1, search_files1);
    TestSplayTree(input_file2, delete_file2, search_files2);

    // Treap
    cout << "\n--------Treap--------" << endl;
    TestTreap(input_file1, delete_file1, search_files1);
    TestTreap(input_file2, delete_file2, search_files2);

    return 0;

}

Node::Node() { left = nullptr; right = nullptr; parent = nullptr; }
Node::Node(Node* left, Node* right) { this->left = left; this->right = right; parent = nullptr; }
Node::Node(Node* left, Node* right, Node* parent) { this->left = left; this->right = right; this->parent = parent; }
Node::Node(string word) : Node() { this->word = removePunctuationAndLower(word); }
Node::Node(Node* left, Node* right, string word) : Node(left, right) { this->word = removePunctuationAndLower(word); }
Node::Node(Node* left, Node* right, Node* parent, string word) : Node(left, right, parent) { this->word = removePunctuationAndLower(word); }
Node* Node::getLeft() { return left; }
Node* Node::getRight() { return right; }
Node* Node::getParent() { return parent; }
string Node::getWord() { return word; }
void Node::setLeft(Node* left) { this->left = left; }
void Node::setRight(Node* right) { this->right = right; }
void Node::setParent(Node* parent) { this->parent = parent; }
void Node::setWord(string word) { this->word = removePunctuationAndLower(word); }
void Node::setWord(Node* node) { this->word = node->getWord(); }
void Node::setWord(Node node) { this->word = node.getWord(); }

string Node::removePunctuationAndLower(string word) {
    word.erase(remove_if(word.begin(), word.end(), static_cast<int(*)(int)>(&ispunct)), word.end());
    transform(word.begin(), word.end(), word.begin(), static_cast<int(*)(int)>(&tolower));
    return word;
}

TreapNode::TreapNode() : Node() { priority = rand(); }
TreapNode::TreapNode(TreapNode* left, TreapNode* right) : Node(left, right) { priority = rand(); }
TreapNode::TreapNode(TreapNode* left, TreapNode* right, TreapNode* parent) : Node(left, right, parent) { priority = rand(); }
TreapNode::TreapNode(string word) : Node(word) { priority = rand(); }
TreapNode::TreapNode(TreapNode* left, TreapNode* right, string word) : Node(left, right, word) { priority = rand(); }
TreapNode::TreapNode(TreapNode* left, TreapNode* right, TreapNode* parent, string word) : Node(left, right, parent, word) { priority = rand(); }
TreapNode* TreapNode::getLeft() { return static_cast<TreapNode*>(Node::getLeft()); }
TreapNode* TreapNode::getRight() { return static_cast<TreapNode*>(Node::getRight()); }
TreapNode* TreapNode::getParent() { return static_cast<TreapNode*>(Node::getParent()); }
int TreapNode::getPriority() { return priority; }
void TreapNode::setLeft(TreapNode* left) { Node::setLeft(static_cast<Node*>(left)); }
void TreapNode::setRight(TreapNode* right) { Node::setRight(static_cast<Node*>(right)); }
void TreapNode::setParent(TreapNode* parent) { Node::setRight(static_cast<Node*>(parent)); }
void TreapNode::setPriority(int priority) { this->priority = priority; }

BinarySearchTree::BinarySearchTree() { root = nullptr; }
BinarySearchTree::BinarySearchTree(Node* root) { this->root = root; }
void BinarySearchTree::Insert(Node* node) {
    if(root == nullptr) root = node;
    else insert(root, node);
}
void BinarySearchTree::Insert(Node node) { Insert(&node); }
void BinarySearchTree::Insert(string word) { Insert(new Node(word)); }
Node* BinarySearchTree::Search(Node* node) { if(root) return search(root, node->getWord()); else return root; }
Node*  BinarySearchTree::Search(Node node) { if(root) return search(root, node.getWord()); else return root; }
Node*  BinarySearchTree::Search(string word) { if(root) return search(root, Node(word).getWord()); else return root; }
void BinarySearchTree::Delete(Node* node) { root = del(root, node->getWord()); }
void BinarySearchTree::Delete(Node node) { root = del(root, node.getWord()); }
void BinarySearchTree::Delete(string word) { root = del(root, word); }

Node* BinarySearchTree::getRoot() { return root; }
void BinarySearchTree::setRoot(Node* root) { this->root = root; }
Node* BinarySearchTree::minValueNode(Node* from) {
    if(from->getLeft()) return minValueNode(from->getLeft());
    else return from;
}

void BinarySearchTree::insert(Node* now, Node* node) {

    if(node->getWord() < now->getWord()) {
        if(now->getLeft()) insert(now->getLeft(), node);
        else { now->setLeft(node); node->setParent(now); }
    }
    else if(node->getWord() > now->getWord()) {
        if(now->getRight()) insert(now->getRight(), node);
        else { now->setRight(node); node->setParent(now); }
    }
    else return;

}
Node* BinarySearchTree::search(Node* now, string word) {

    if(word < now->getWord()) {
        if(now->getLeft()) return search(now->getLeft(), word);
        else return nullptr;
    }
    else if(word > now->getWord()) {
        if(now->getRight()) return search(now->getRight(), word);
        else return nullptr;
    }

    return now;

}
Node* BinarySearchTree::del(Node* now, string word) {

    if(now == nullptr) return now;

    if(word < now->getWord()) now->setLeft(del(now->getLeft(), word));
    else if(word > now->getWord()) now->setRight(del(now->getRight(), word));
    else {

        if(now->getLeft() == nullptr) {
            Node* temp = now->getRight();
            free(now);
            return temp;
        }
        else if(now->getRight() == nullptr) {
            Node* temp = now->getLeft();
            free(now);
            return temp;
        }
        
        Node* temp = minValueNode(now->getRight());
        now->setWord(temp);
        now->setRight(del(now->getRight(), temp->getWord()));

    }

    return now;

}

SplayTree::SplayTree() : BinarySearchTree() {};
SplayTree::SplayTree(Node* root) : BinarySearchTree(root) {};
void SplayTree::Insert(Node* node) { BinarySearchTree::Insert(node); splay(node); }
void SplayTree::Insert(Node node) { Insert(&node); }
void SplayTree::Insert(string word) { Insert(new Node(word)); }
Node* SplayTree::Search(Node* node) { Node* target = BinarySearchTree::Search(node); if(target) splay(target); return target; }
Node* SplayTree::Search(Node node) { Node* target = BinarySearchTree::Search(&node); if(target) splay(target); return target; }
Node* SplayTree::Search(string word) { Node* target = BinarySearchTree::Search(word); if(target) splay(target); return target; }
void SplayTree::Delete(Node* node) { Node* target = Search(node); if(target) del(); }
void SplayTree::Delete(Node node) { Node* target = Search(node); if(target) del(); }
void SplayTree::Delete(string word) { Node* target = Search(word); if(target) del(); }

void SplayTree::splay(Node* node) {

    while(node->getParent()) {
        if(node->getParent()->getParent() == nullptr) {
            if(node->getParent()->getLeft() == node) rightRotation(node->getParent());
            else leftRotation(node->getParent());
        }
        else if(node->getParent()->getParent()->getLeft() == node->getParent() && node->getParent()->getLeft() == node) {
            rightRotation(node->getParent()->getParent());
            rightRotation(node->getParent());
        }
        else if(node->getParent()->getParent()->getRight() == node->getParent() && node->getParent()->getRight() == node) {
            leftRotation(node->getParent()->getParent());
            leftRotation(node->getParent());
        }
        else if(node->getParent()->getParent()->getLeft() == node->getParent() && node->getParent()->getRight() == node) {
            leftRotation(node->getParent());
            rightRotation(node->getParent());
        }
        else {
            rightRotation(node->getParent());
            leftRotation(node->getParent());
        }
    }

}
void SplayTree::leftRotation(Node* rotatedNode) {

    Node* rightChild = rotatedNode->getRight();

    rotatedNode->setRight(rightChild->getLeft());
    if(rightChild->getLeft()) rightChild->getLeft()->setParent(rotatedNode);

    rightChild->setParent(rotatedNode->getParent());
    if(rotatedNode->getParent() == nullptr) setRoot(rightChild);
    else if(rotatedNode == rotatedNode->getParent()->getLeft()) rotatedNode->getParent()->setLeft(rightChild);
    else rotatedNode->getParent()->setRight(rightChild);

    rightChild->setLeft(rotatedNode);
    rotatedNode->setParent(rightChild);

}
void SplayTree::rightRotation(Node* rotatedNode) {

    Node* leftChild = rotatedNode->getLeft();

    rotatedNode->setLeft(leftChild->getRight());
    if(leftChild->getRight()) leftChild->getRight()->setParent(rotatedNode);

    leftChild->setParent(rotatedNode->getParent());
    if(rotatedNode->getParent() == nullptr) setRoot(leftChild);
    else if(rotatedNode == rotatedNode->getParent()->getLeft()) rotatedNode->getParent()->setLeft(leftChild);
    else rotatedNode->getParent()->setRight(leftChild);

    leftChild->setRight(rotatedNode);
    rotatedNode->setParent(leftChild);

}
void SplayTree::del() {

    Node* root = getRoot();

    if(root->getLeft() && root->getRight()) {
        Node* leftChild = root->getLeft();
        setRoot(root->getRight());
        Search(minValueNode(getRoot()));
        getRoot()->setLeft(leftChild);
        leftChild->setParent(getRoot());
    }
    else {
        if(root->getLeft()) setRoot(root->getLeft());
        else setRoot(root->getRight());
    }

    if(getRoot()) getRoot()->setParent(nullptr);
    free(root);

}

Treap::Treap() { root = nullptr; }
Treap::Treap(TreapNode* root) { this->root = root; }
void Treap::Insert(TreapNode* node) { root = insert(root, node); }
void Treap::Insert(TreapNode node) { Insert(&node); }
void Treap::Insert(string word) { Insert(new TreapNode(word)); }
TreapNode* Treap::Search(TreapNode* node) { if(root) return search(root, node->getWord()); else return root; }
TreapNode* Treap::Search(TreapNode node) { if(root) return search(root, node.getWord()); else return root; }
TreapNode* Treap::Search(string word) { if(root) return search(root, Node(word).getWord()); else return root; }
void Treap::Delete(TreapNode* node) { del(root, node->getWord()); }
void Treap::Delete(TreapNode node) { del(root, node.getWord()); }
void Treap::Delete(string word) { del(root, Node(word).getWord()); }

TreapNode* Treap::insert(TreapNode* now, TreapNode* node) {

    if(now == nullptr) return node;

    if(node->getWord() < now->getWord()) {
        now->setLeft(insert(now->getLeft(), node));
        if(now->getLeft()->getPriority() < now->getPriority()) now = rotateRight(now);
    }
    else if(node->getWord() > now->getWord()) {
        now->setRight(insert(now->getRight(), node));
        if(now->getRight()->getPriority() < now->getPriority()) now = rotateLeft(now);
    }

    return now;

}
TreapNode* Treap::search(TreapNode* now, string word) {

    if(word < now->getWord()) {
        if(now->getLeft()) return search(now->getLeft(), word);
        else return nullptr;
    }
    else if(word > now->getWord()) {
        if(now->getRight()) return search(now->getRight(), word);
        else return nullptr;
    }

    return now;

}
TreapNode* Treap::del(TreapNode* node, string word) {

    if (node == nullptr) return nullptr;

    if (word < node->getWord()) node->setLeft(del(node->getLeft(), word));
    else if (word > node->getWord()) node->setRight(del(node->getRight(), word));
    else {
        if (node->getLeft() == nullptr) {
            TreapNode* temp = node->getRight();
            delete node;
            return temp;
        }
        else if (node->getRight() == nullptr) {
            TreapNode* temp = node->getLeft();
            delete node;
            return temp;
        }

        if (node->getLeft()->getPriority() > node->getRight()->getPriority()) {
            node = rotateRight(node);
            node->setRight(del(node->getRight(), word));
        } else {
            node = rotateLeft(node);
            node->setLeft(del(node->getLeft(), word));
        }
    }

    return node;
}
TreapNode* Treap::rotateLeft(TreapNode* node) {
    TreapNode *rightChild = node->getRight(), *T2 = rightChild->getLeft();
    rightChild->setLeft(node);
    node->setRight(T2);
    return rightChild;
}
TreapNode* Treap::rotateRight(TreapNode* node) {
    TreapNode *leftChild = node->getLeft(), *T2 = leftChild->getRight();
    leftChild->setRight(node);
    node->setLeft(T2);
    return leftChild;
}

void TestBinarySearchTree(string& input_file, string& delete_file, vector<string>& search_file_vector) {

    BinarySearchTree dict;
    
    ifstream ifs;
    ifs.open(input_file);
    cout << "Input file: " << input_file << endl << endl;

    string line;
    while (getline(ifs, line)) {
        stringstream ss(line);
        string word;
        while (ss >> word) dict.Insert(word);
    }
    ifs.close();

    for(string search_file : search_file_vector) {

        int success = 0, fail = 0;
        ifs.open(search_file);
        auto start = chrono::high_resolution_clock::now();

        while (getline(ifs, line)) {
            stringstream ss(line);
            string word;
            while (ss >> word) {
                if(dict.Search(word)) success++;
                else fail++;
            }
        }
        ifs.close();

        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        cout << "Search file: " << search_file << endl;
        cout << "Success rate: " << 100.0 * success / (success + fail) << '%' << endl;
        cout << "Failure rate: " << 100.0 * fail / (success + fail) << '%' << endl;
        cout << "Execution time: " << duration.count() << " mus" << endl << endl;

    }

    ifs.open(delete_file);
    while (getline(ifs, line)) {
        stringstream ss(line);
        string word;
        while (ss >> word) {
            cout << "Delete: " << word << endl;
            dict.Delete(word);
            cout << "Search " << word << ": " << dict.Search(word) << endl << endl;
        }
    }
    ifs.close();

}

void TestSplayTree(string& input_file, string& delete_file, vector<string>& search_file_vector) {

    SplayTree dict;
    
    ifstream ifs;
    ifs.open(input_file);
    cout << "Input file: " << input_file << endl << endl;

    string line;
    while (getline(ifs, line)) {
        stringstream ss(line);
        string word;
        while (ss >> word) dict.Insert(word);
    }
    ifs.close();

    for(string search_file : search_file_vector) {

        int success = 0, fail = 0;
        ifs.open(search_file);
        auto start = chrono::high_resolution_clock::now();

        while (getline(ifs, line)) {
            stringstream ss(line);
            string word;
            while (ss >> word) {
                if(dict.Search(word)) success++;
                else fail++;
            }
        }
        ifs.close();

        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        cout << "Search file: " << search_file << endl;
        cout << "Success rate: " << 100.0 * success / (success + fail) << '%' << endl;
        cout << "Failure rate: " << 100.0 * fail / (success + fail) << '%' << endl;
        cout << "Execution time: " << duration.count() << " mus" << endl << endl;

    }

    ifs.open(delete_file);
    while (getline(ifs, line)) {
        stringstream ss(line);
        string word;
        while (ss >> word) {
            cout << "Delete: " << word << endl;
            dict.Delete(word);
            cout << "Search " << word << ": " << dict.Search(word) << endl << endl;
        }
    }
    ifs.close();

}

void TestTreap(string& input_file, string& delete_file, vector<string>& search_file_vector) {

    Treap dict;
    
    ifstream ifs;
    ifs.open(input_file);
    cout << "Input file: " << input_file << endl << endl;

    string line;
    while (getline(ifs, line)) {
        stringstream ss(line);
        string word;
        while (ss >> word) dict.Insert(word);
    }
    ifs.close();

    for(string search_file : search_file_vector) {

        int success = 0, fail = 0;
        ifs.open(search_file);
        auto start = chrono::high_resolution_clock::now();

        while (getline(ifs, line)) {
            stringstream ss(line);
            string word;
            while (ss >> word) {
                if(dict.Search(word)) success++;
                else fail++;
            }
        }
        ifs.close();

        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        cout << "Search file: " << search_file << endl;
        cout << "Success rate: " << 100.0 * success / (success + fail) << '%' << endl;
        cout << "Failure rate: " << 100.0 * fail / (success + fail) << '%' << endl;
        cout << "Execution time: " << duration.count() << " mus" << endl << endl;

    }

    ifs.open(delete_file);
    while (getline(ifs, line)) {
        stringstream ss(line);
        string word;
        while (ss >> word) {
            cout << "Delete: " << word << endl;
            dict.Delete(word);
            cout << "Search " << word << ": " << dict.Search(word) << endl << endl;
        }
    }
    ifs.close();

}