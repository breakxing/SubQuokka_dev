#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <iterator>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

class Node {
public:
    //set<int> dependency = {};
    set<int> myQubit = {};
    set<int> otherDependencies = {};
    set<string> dependenciesName = {};
    vector<Node *> preNodes;
    vector<Node *> nextNodes;
    string name;
    Node *next;
    Node()=default;
};

void readCircuitToVec(char *(&), vector<string> &);
void setGateInfo(const string &, vector<string> &);
void setDAG(vector<string> &, vector<Node*> &, Node *(&));
void printData(Node *(&));
void findChunkbyMaxNodeOri(int, Node *(&), set<int> &);
void findChunkbyMaxNode(int, Node *(&), set<int> &);
void popGate(Node *(&), set<int> &, vector<string> &, vector<string> &);
void freeNodes(Node *(&));

void readCircuitToVec(char *(& filename), vector<string> &gateVec) {
    ifstream cirfile;
    cirfile.open(filename);
    string line;
    int iter = 0;
    while (getline(cirfile, line)) {
        stringstream ss(line);
        string gateType;
        ss >> gateType;
        line += " " + gateType + to_string(iter);
        gateVec.push_back(line);
        iter++;
    }
    cirfile.close();
}

void setGateInfo(const string &curGate, vector<string> &gateInfo, int iter) {
    size_t pos = 0;
    string token;
    stringstream ss(curGate);
    istream_iterator<string> begin(ss);
    istream_iterator<string> end;
    vector<string> tokens(begin, end);
    gateInfo.assign(tokens.begin(), tokens.end());
    //for (auto &s: gateInfo)
    //    cout << s << " ";
    //cout << endl;
}

void printObj(const Node &ptr) {
    for (const auto &d: ptr.nextNodes) {
        cout << d->name << ", [dep]: "; // << d->dependenciesName << endl;
        for (const string &name : d->dependenciesName) {
            cout << name << " ";
        }
        if(d->nextNodes.size() > 0)
            cout << "n[0]: " << d->nextNodes[0]->name;
        if(d->nextNodes.size() > 1)
            cout << "====>n[1]: " << d->nextNodes[1]->name;
        cout << endl;
        printObj(*d);
    }
}

void setDAG(vector<string> &gateVec, vector<Node*> &tails, Node *(&nList)) {
    Node *pre = nList;
    for(int i = 0; i < gateVec.size(); i++) {
        vector<string> gateInfo;
        setGateInfo(gateVec[i], gateInfo, i);
        Node *node = new Node();
        pre->next = node;
        pre = node;
    }
    pre = nList;

    for(int i = 0; i < gateVec.size(); i++) {
        vector<string> gateInfo;
        setGateInfo(gateVec[i], gateInfo, i);
        Node *node = pre->next; //Node *node = new Node();
        if(gateInfo[0] == "H") {
            int targ = stoi(gateInfo[1]);
            // setup current node
            node->myQubit.insert(targ); // set<int> myQubit = {};
            node->dependenciesName.insert(gateInfo[gateInfo.size() - 1]); // set<string> dependenciesName = {};
            node->name = gateInfo[gateInfo.size() - 1];  // string name;
            // get the pre info
            Node *pre = tails[targ];
            set_union(
                pre->myQubit.begin(), pre->myQubit.end(),
                pre->otherDependencies.begin(), pre->otherDependencies.end(),
                inserter(node->otherDependencies, node->otherDependencies.begin()));
            for (const int &n : node->myQubit) {
                node->otherDependencies.erase(n);
            }
            node->dependenciesName.insert(pre->dependenciesName.begin(), pre->dependenciesName.end());
            node->preNodes.push_back(pre);
            pre->nextNodes.push_back(node);
            tails[targ] = node;
        }
        else if(gateInfo[0] == "CZ") {
            int targ1 = stoi(gateInfo[1]);
            int targ2 = stoi(gateInfo[2]);
            // setup current node
            node->myQubit.insert(targ1); // set<int> myQubit = {};
            node->myQubit.insert(targ2);
            node->dependenciesName.insert(gateInfo[gateInfo.size() - 1]); // set<string> dependenciesName = {};
            node->name = gateInfo[gateInfo.size() - 1];  // string name;

            // get the pre info
            Node *pre1 = tails[targ1];
            Node *pre2 = tails[targ2];

            set_union(
                pre1->myQubit.begin(), pre1->myQubit.end(),
                pre1->otherDependencies.begin(), pre1->otherDependencies.end(),
                inserter(node->otherDependencies, node->otherDependencies.begin()));
            set_union(
                pre2->myQubit.begin(), pre2->myQubit.end(),
                pre2->otherDependencies.begin(), pre2->otherDependencies.end(),
                inserter(node->otherDependencies, node->otherDependencies.begin()));
            for (const int &n : node->myQubit) {
                node->otherDependencies.erase(n);
            }
            node->dependenciesName.insert(pre1->dependenciesName.begin(), pre1->dependenciesName.end());
            node->dependenciesName.insert(pre2->dependenciesName.begin(), pre2->dependenciesName.end());

            node->preNodes.push_back(pre1);
            node->preNodes.push_back(pre2);
            pre1->nextNodes.push_back(node);
            pre2->nextNodes.push_back(node);
            tails[targ1] = node;
            tails[targ2] = node;
        }
        else if(gateInfo[0] == "TOF") {
            int targ1 = stoi(gateInfo[1]);
            int targ2 = stoi(gateInfo[2]);
            int targ3 = stoi(gateInfo[3]);
            // setup current node
            node->myQubit.insert(targ1); // set<int> myQubit = {};
            node->myQubit.insert(targ2);
            node->myQubit.insert(targ3);
            node->dependenciesName.insert(gateInfo[gateInfo.size() - 1]); // set<string> dependenciesName = {};
            node->name = gateInfo[gateInfo.size() - 1];  // string name;

            // get the pre info
            Node *pre1 = tails[targ1];
            Node *pre2 = tails[targ2];
            Node *pre3 = tails[targ3];

            set_union(
                pre1->myQubit.begin(), pre1->myQubit.end(),
                pre1->otherDependencies.begin(), pre1->otherDependencies.end(),
                inserter(node->otherDependencies, node->otherDependencies.begin()));
            set_union(
                pre2->myQubit.begin(), pre2->myQubit.end(),
                pre2->otherDependencies.begin(), pre2->otherDependencies.end(),
                inserter(node->otherDependencies, node->otherDependencies.begin()));
            set_union(
                pre3->myQubit.begin(), pre3->myQubit.end(),
                pre3->otherDependencies.begin(), pre3->otherDependencies.end(),
                inserter(node->otherDependencies, node->otherDependencies.begin()));
            for (const int &n : node->myQubit) {
                node->otherDependencies.erase(n);
            }
            node->dependenciesName.insert(pre1->dependenciesName.begin(), pre1->dependenciesName.end());
            node->dependenciesName.insert(pre2->dependenciesName.begin(), pre2->dependenciesName.end());
            node->dependenciesName.insert(pre3->dependenciesName.begin(), pre3->dependenciesName.end());

            node->preNodes.push_back(pre1);
            node->preNodes.push_back(pre2);
            node->preNodes.push_back(pre3);
            pre1->nextNodes.push_back(node);
            pre2->nextNodes.push_back(node);
            pre3->nextNodes.push_back(node);
            tails[targ1] = node;
            tails[targ2] = node;
            tails[targ3] = node;
        }
        pre = pre->next;
    }
}

void printData(Node *(&nList)) {
    cout << "========" << endl;
    Node *head = nList;
    while(nList->next) {
        cout << nList->next->name << ": ";
        nList = nList->next;
        set<int> tmpSet = {};
        set_union(
            nList->myQubit.begin(), nList->myQubit.end(),
            nList->otherDependencies.begin(), nList->otherDependencies.end(),
            inserter(tmpSet, tmpSet.begin()));
        for (const int &n : tmpSet) {
            cout << n << " ";
        }
        // cout << endl;
        cout << "[Dep]: ";
        for (const string &n : nList->dependenciesName) {
            cout << n << " ";
        }
        cout << "[MaxNodes]: " << nList->dependenciesName.size();
        cout << endl;
    }
    nList = head; // go back to head;
}

// TODO: 沒有合併的接腳時，有bug
void findChunkbyMaxNodeOri(int chunkSize, Node *(&nList), set<int> &chunkSet) {
    Node *head = nList;
    int maxDepNameSize = 0;
    while(nList->next) {
        nList = nList->next;
        set<int> tmpSet = {};
        int curNodeSize = 0;
        set_union(
            nList->myQubit.begin(), nList->myQubit.end(),
            nList->otherDependencies.begin(), nList->otherDependencies.end(),
            inserter(tmpSet, tmpSet.begin()));
        curNodeSize = tmpSet.size();
        if(curNodeSize <= chunkSize) {
            int curDependenciesNameSize = nList->dependenciesName.size();
            if (maxDepNameSize <= curDependenciesNameSize) {
                maxDepNameSize = curDependenciesNameSize;
                chunkSet = tmpSet;
            }
        }
    }
    nList = head; // go back to head;
}

void findChunkbyMaxNode(int chunkSize, Node *(&nList), set<int> &chunkSet) {
    Node *head = nList;
    int maxNodeSize = 0;
    chunkSet = {};
    set<string> names = {};
    int changeNames = 1;
    while(maxNodeSize < chunkSize && changeNames == 1) {
        changeNames = 0;
        int tmpDepNameSize = 0;
        set<int> tmpChunkSet = {};
        string tmpName = {};
        while(nList->next) {
            nList = nList->next;
            set<int> tmpSet = {};
            int tmpNodeSize = 0;
            set_union(
                nList->myQubit.begin(), nList->myQubit.end(),
                nList->otherDependencies.begin(), nList->otherDependencies.end(),
                inserter(tmpSet, tmpSet.begin()));
            tmpNodeSize = tmpSet.size();
            if(tmpNodeSize <= chunkSize - maxNodeSize) {
                int curDependenciesNameSize = nList->dependenciesName.size();
                if (tmpDepNameSize <= curDependenciesNameSize &&
                    names.find(nList->name) == names.end())
                {
                    tmpDepNameSize = curDependenciesNameSize;
                    tmpChunkSet = tmpSet;
                    tmpName = nList->name;
                    changeNames = 1;

                    
                }
            }
        }
        nList = head; // go back to head;
        chunkSet.insert(tmpChunkSet.begin(), tmpChunkSet.end());
        names.insert(tmpName);
        maxNodeSize = chunkSet.size();
    }
}

void popGate(Node *(&nList), set<int> &chunkSet, vector<string> &gateVec, vector<string> &popGateVec) {
    Node *head = nList;
    while(nList->next) {
        Node *tmpPtr = nList->next;
        set<int> tmpSet = {};
        set_union(
            tmpPtr->myQubit.begin(), tmpPtr->myQubit.end(),
            tmpPtr->otherDependencies.begin(), tmpPtr->otherDependencies.end(),
            inserter(tmpSet, tmpSet.begin()));
        if (includes(chunkSet.begin(), chunkSet.end(),
            tmpSet.begin(), tmpSet.end())) {
            for(int i = 0; i < gateVec.size(); i++) {
                vector<string> gateInfo;
                setGateInfo(gateVec[i], gateInfo, i);
                string tmpName = gateInfo[gateInfo.size() - 1];
                if(tmpName == tmpPtr->name) {
                    gateVec.erase(gateVec.begin() + i);
                    popGateVec.push_back(tmpName);
                    break;
                }
            }
            nList->next = nList->next->next;
            delete tmpPtr; // delete the Node
        }
        else{
            nList = nList->next;
        }
    }
    nList = head; // go back to head;
}

void freeNodes(Node *(&nList)){
    while(nList->next) {
        Node *tmpPtr = nList->next;
        nList->next = nList->next->next;
        delete tmpPtr; // delete the Node
    }
}

int main(int argc, char *argv[]) {
    vector<string> gateVec;
    readCircuitToVec(argv[1], gateVec);
    size_t qubits = 40;
    vector<Node> dags (qubits);
    vector<Node*> tails (qubits);
    for(int i = 0; i < qubits; i++) {
        tails[i] = &(dags[i]);
    }

    Node *nList = new Node();
    Node *head = nList; // backup the head;
    setDAG(gateVec, tails, nList);
    //set<int> chunkSet {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; // set<int> chunkSet {0, 1, 2, 3};
    set<int> chunkSet {0, 1}; // small set
    int chunkSize = chunkSet.size();
    int iter = 0;
    while(nList->next) {
        //printData(nList);
        vector<string> popGateVec;
        if(iter != 0)
            findChunkbyMaxNode(chunkSize, nList, chunkSet);
   
        // Displaying set elements
        //for (set<int>::iterator itr = chunkSet.begin(); itr != chunkSet.end(); itr++) {
        //    cout << *itr << " ";
        //}
        //cout << endl;

        popGate(nList, chunkSet, gateVec, popGateVec);

        // reset all Nodes
        for(int i = 0; i < qubits; i++) { // reset tail to head of array
            tails[i] = &(dags[i]);
        }
        freeNodes(nList);
        nList = new Node();
        head = nList; // backup the head;

        setDAG(gateVec, tails, nList);
        iter++;

        // print
        cout << "[iter " << iter << "]" << " ";
        for(int i = 0; i < popGateVec.size(); i++)
            cout << popGateVec[i] << " ";
        if(iter != 1) {
            cout << ", Swap in chunk: ";
            for(auto &c: chunkSet)
                cout << c << " ";
        }
        cout << endl;
        nList = head;
    }
}