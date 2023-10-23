// backup the source code
#include <iostream>
#include <iomanip>
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
    string type;
    string name;
    Node *next;
    Node()=default;
    Node(const Node& other) {
        myQubit = other.myQubit;
        otherDependencies = other.otherDependencies;
        dependenciesName = other.dependenciesName;
        preNodes = other.preNodes;
        nextNodes = other.nextNodes;
        type = other.type;
        name = other.name;
        next = other.next;
    }
};

void readCircuitToVec(char* &, vector<string> &);
void setGateInfo(const string &, vector<string> &);
void setDagFromVec(vector<string> &, vector<Node*> &, Node* &);

void copyNodes(Node* &, Node* &);
void freeNodes(Node* &);
void printNodes(Node* &);
void popNodes(Node* &, set<int> &, vector<string> &);
void setDagFromNodes(vector<Node> &, vector<Node*> &, Node* &);
void findChunkbyMaxNodePrinciple(const int, Node* &, set<int> &, int &);
void findChunkbyMaxNodeWithNextPrinciple(const int, Node* &, set<int> &, vector<Node> &, vector<Node*> &);

void readCircuitToVec(char* &filename, vector<string> &gateVec) {
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

void setDagFromVec(vector<string> &gateVec, vector<Node*> &tails, Node* &nList) {
    Node *pre = nList;
    for(int i = 0; i < gateVec.size(); i++) {
        vector<string> gateInfo;
        setGateInfo(gateVec[i], gateInfo, i);
        Node *node = new Node();

        if(gateInfo[0] == "H") {
            int targ = stoi(gateInfo[1]);
            // setup current node
            node->myQubit.insert(targ); // set<int> myQubit = {};
            node->dependenciesName.insert(gateInfo[gateInfo.size() - 1]); // set<string> dependenciesName = {};
            node->type = gateInfo[0];
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
        else if(gateInfo[0] == "CZ" || gateInfo[0] == "U2" || gateInfo[0] == "SWAP") {
            int targ1 = stoi(gateInfo[1]);
            int targ2 = stoi(gateInfo[2]);
            // setup current node
            node->myQubit.insert(targ1); // set<int> myQubit = {};
            node->myQubit.insert(targ2);
            node->dependenciesName.insert(gateInfo[gateInfo.size() - 1]); // set<string> dependenciesName = {};
            node->type = gateInfo[0];            
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
            node->type = gateInfo[0];
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
        pre->next = node;
        pre = node;
    }
    pre = nList;
}

void setDagFromNodes(vector<Node> &dags, vector<Node*> &tails, Node* &nList) {
    for(int i = 0; i < dags.size(); i++) // reset tail to head of array
       tails[i] = &(dags[i]);
    Node *ptr = nList;
    while(ptr->next) {
        Node* tmpPtr = ptr->next;
        if(tmpPtr->type == "H") {
            vector<int> targs;
            for(auto &s: tmpPtr->myQubit)
                targs.push_back(s);
            tmpPtr->otherDependencies = {};
            tmpPtr->dependenciesName = {};
            tmpPtr->dependenciesName.insert(tmpPtr->name);
            Node *pre = tails[targs[0]];
            set_union(
                pre->myQubit.begin(), pre->myQubit.end(),
                pre->otherDependencies.begin(), pre->otherDependencies.end(),
                inserter(tmpPtr->otherDependencies, tmpPtr->otherDependencies.begin()));
            for (auto &n: tmpPtr->myQubit) {
                tmpPtr->otherDependencies.erase(n);
            }
            tmpPtr->dependenciesName.insert(pre->dependenciesName.begin(), pre->dependenciesName.end());
            tmpPtr->preNodes.push_back(pre);
            pre->nextNodes.push_back(tmpPtr);
            tails[targs[0]] = tmpPtr;
            //node->myQubit.insert(targ); // set<int> myQubit = {};
            //node->dependenciesName.insert(gateInfo[gateInfo.size() - 1]); // set<string> dependenciesName = {};
            //node->name = gateInfo[gateInfo.size() - 1];  // string name;
            // get the pre info
            //Node *pre = tails[targ];
        }
        else if(tmpPtr->type == "CZ" || tmpPtr->type == "U2" || tmpPtr->type == "SWAP") {
            vector<int> targs;
            for(auto &s: tmpPtr->myQubit)
                targs.push_back(s);
            tmpPtr->otherDependencies = {};
            tmpPtr->dependenciesName = {};
            tmpPtr->dependenciesName.insert(tmpPtr->name);
            Node *pre1 = tails[targs[0]];
            Node *pre2 = tails[targs[1]];

            set_union(
                pre1->myQubit.begin(), pre1->myQubit.end(),
                pre1->otherDependencies.begin(), pre1->otherDependencies.end(),
                inserter(tmpPtr->otherDependencies, tmpPtr->otherDependencies.begin()));
            set_union(
                pre2->myQubit.begin(), pre2->myQubit.end(),
                pre2->otherDependencies.begin(), pre2->otherDependencies.end(),
                inserter(tmpPtr->otherDependencies, tmpPtr->otherDependencies.begin()));

            for (auto &n: tmpPtr->myQubit) {
                tmpPtr->otherDependencies.erase(n);
            }
            tmpPtr->dependenciesName.insert(pre1->dependenciesName.begin(), pre1->dependenciesName.end());
            tmpPtr->dependenciesName.insert(pre2->dependenciesName.begin(), pre2->dependenciesName.end());

            tmpPtr->preNodes.push_back(pre1);
            tmpPtr->preNodes.push_back(pre2);
            pre1->nextNodes.push_back(tmpPtr);
            pre2->nextNodes.push_back(tmpPtr);
            tails[targs[0]] = tmpPtr;
            tails[targs[1]] = tmpPtr;
        }
        else if(tmpPtr->type ==  "TOF") {
            vector<int> targs;
            for(auto &s: tmpPtr->myQubit)
                targs.push_back(s);
            tmpPtr->otherDependencies = {};
            tmpPtr->dependenciesName = {};
            tmpPtr->dependenciesName.insert(tmpPtr->name);
            Node *pre1 = tails[targs[0]];
            Node *pre2 = tails[targs[1]];
            Node *pre3 = tails[targs[2]];

            set_union(
                pre1->myQubit.begin(), pre1->myQubit.end(),
                pre1->otherDependencies.begin(), pre1->otherDependencies.end(),
                inserter(tmpPtr->otherDependencies, tmpPtr->otherDependencies.begin()));
            set_union(
                pre2->myQubit.begin(), pre2->myQubit.end(),
                pre2->otherDependencies.begin(), pre2->otherDependencies.end(),
                inserter(tmpPtr->otherDependencies, tmpPtr->otherDependencies.begin()));
            set_union(
                pre3->myQubit.begin(), pre3->myQubit.end(),
                pre3->otherDependencies.begin(), pre3->otherDependencies.end(),
                inserter(tmpPtr->otherDependencies, tmpPtr->otherDependencies.begin()));

            for (auto &n: tmpPtr->myQubit) {
                tmpPtr->otherDependencies.erase(n);
            }
            tmpPtr->dependenciesName.insert(pre1->dependenciesName.begin(), pre1->dependenciesName.end());
            tmpPtr->dependenciesName.insert(pre2->dependenciesName.begin(), pre2->dependenciesName.end());
            tmpPtr->dependenciesName.insert(pre3->dependenciesName.begin(), pre3->dependenciesName.end());
            tmpPtr->preNodes.push_back(pre1);
            tmpPtr->preNodes.push_back(pre2);
            tmpPtr->preNodes.push_back(pre3);
            pre1->nextNodes.push_back(tmpPtr);
            pre2->nextNodes.push_back(tmpPtr);
            pre3->nextNodes.push_back(tmpPtr);
            tails[targs[0]] = tmpPtr;
            tails[targs[1]] = tmpPtr;
            tails[targs[2]] = tmpPtr;
        }    
        ptr = ptr->next;
    }
}

void copyNodes(Node* &nList, Node* &nListSecond) {
    Node *ptr = nList;
    Node *ptrCopy = nListSecond;
    while(ptr->next) {
        ptrCopy->next = new Node(*(ptr->next));
        ptrCopy = ptrCopy->next;
        ptr = ptr->next;
    }
}

void freeNodes(Node* &nList) {
    while(nList->next) {
        Node *tmpPtr = nList->next;
        nList->next = nList->next->next;
        delete tmpPtr; // delete the Node
    }
}

void printNodes(Node* &nList) {
    cout << "========" << endl;
    Node *ptr = nList;
    while(ptr->next) {
        cout << ptr->next->name << ": ";
        ptr = ptr->next;
        set<int> tmpSet = {};
        set_union(
            ptr->myQubit.begin(), ptr->myQubit.end(),
            ptr->otherDependencies.begin(), ptr->otherDependencies.end(),
            inserter(tmpSet, tmpSet.begin()));
        for (const int &n : tmpSet) {
            cout << n << " ";
        }
        // cout << endl;
        cout << "[Dep]: ";
        for (const string &n : ptr->dependenciesName) {
            cout << n << " ";
        }
        cout << "[MaxNodes]: " << ptr->dependenciesName.size();
        cout << endl;
    }
}

void popNodes(Node* &nList, set<int> &chunkSet, vector<string> &popGateVec) {
    Node *ptr = nList;
    while(ptr->next) {
        Node *tmpPtr = ptr->next;
        set<int> tmpSet = {};
        set_union(
            tmpPtr->myQubit.begin(), tmpPtr->myQubit.end(),
            tmpPtr->otherDependencies.begin(), tmpPtr->otherDependencies.end(),
            inserter(tmpSet, tmpSet.begin())
        );
        if (includes(chunkSet.begin(), chunkSet.end(), tmpSet.begin(), tmpSet.end())) 
        {
            string tmpName = tmpPtr->name;
            popGateVec.push_back(tmpName);
            ptr->next = ptr->next->next;
            delete tmpPtr; // delete the Node
        }
        else
            ptr = ptr->next;
    }
}

void findChunkbyMaxNodePrinciple(const int chunkSize, Node *(&nList), set<int> &chunkSet, int &maxNameSize) {
    Node *nListSecond = new Node();
    copyNodes(nList, nListSecond);
    Node *ptr = nListSecond;
    int maxNodeSize = 0; // total node size in this iteration
    maxNameSize = 0; // total dep name in this iteration
    chunkSet = {};
    set<string> nameSet = {};
    int changeNames = 1;
    while(maxNodeSize < chunkSize && changeNames == 1) {
        changeNames = 0;
        int tmpNameSize = 0; // 目前的 name
        set<int> tmpChunkSet = {};
        set<string> tmpNameSet = {};
        string tmpName = "";
        while(ptr->next) {
            Node *tmpPtr = ptr->next;
            set<int> tmpSet = {};
            set_union(
                tmpPtr->myQubit.begin(), tmpPtr->myQubit.end(),
                tmpPtr->otherDependencies.begin(), tmpPtr->otherDependencies.end(),
                inserter(tmpSet, tmpSet.begin()));
            tmpSet.insert(chunkSet.begin(), chunkSet.end());
            auto tmpNodeSize = tmpSet.size();
            if(tmpNodeSize <= chunkSize) {
            //if(tmpNodeSize + maxNodeSize <= chunkSize) { // 若能裝進 chunk 中
                int curNameSize = tmpPtr->dependenciesName.size();
                if (tmpNameSize <= curNameSize)
                {
                    tmpNameSize = curNameSize;
                    tmpChunkSet = tmpSet;
                    tmpName = tmpPtr->name;
                    changeNames = 1;
                    tmpNameSet = tmpPtr->dependenciesName;
                }
            }
            ptr = ptr->next;
        }
        vector<string> popGateVec;
        popNodes(nListSecond, tmpChunkSet, popGateVec);
        chunkSet.insert(tmpChunkSet.begin(), tmpChunkSet.end());
        nameSet.insert(tmpNameSet.begin(), tmpNameSet.end());
        maxNodeSize = chunkSet.size();
        maxNameSize = nameSet.size();
        ptr = nListSecond; // back to head
    }
    freeNodes(nListSecond);
}

void findChunkbyMaxNodeWithNextPrinciple(const int chunkSize, Node* &nList, set<int> &chunkSet, vector<Node> &dags, vector<Node*> &tails) {
    Node *nListSecond = new Node();
    copyNodes(nList, nListSecond);
    Node *ptr = nListSecond;
    int maxNodeSize = 0; // total node size in this iteration
    int maxNameSize = 0; // total dep name in this iteration
    chunkSet = {};
    int changeNames = 1;
    while(maxNodeSize < chunkSize && changeNames == 1) {
        changeNames = 0;
        int tmpNameSize = 0; // 目前的 name 的上限
        set<int> tmpChunkSet = {};
        string tmpName = "";
        while(ptr->next) {
            Node *tmpPtr = ptr->next;
            set<int> tmpSet = {};
            set_union(
                tmpPtr->myQubit.begin(), tmpPtr->myQubit.end(),
                tmpPtr->otherDependencies.begin(), tmpPtr->otherDependencies.end(),
                inserter(tmpSet, tmpSet.begin()));
            // auto tmpNodeSize = tmpSet.size();
            tmpSet.insert(chunkSet.begin(), chunkSet.end());
            auto tmpNodeSize = tmpSet.size();
            if(tmpNodeSize <= chunkSize) {
            //if(tmpNodeSize + maxNodeSize <= chunkSize) { // 若能裝進 chunk 中
                int curNameSize = tmpPtr->dependenciesName.size();

                Node *nListNextIter = new Node();
                copyNodes(nListSecond, nListNextIter);
                vector<string> popGateVecNextIter;
                popNodes(nListNextIter, tmpSet, popGateVecNextIter);
                setDagFromNodes(dags, tails, nListNextIter);

                //cout << "\ncurNode: " << tmpPtr->name << endl;
                //printNodes(nListNextIter);
                set <int> chunkSetNextIter = {};
                int nameNextIterSize = 0;

                findChunkbyMaxNodePrinciple(chunkSize, nListNextIter, chunkSetNextIter, nameNextIterSize);
                //cout << "(curName, nextName):" << curNameSize << ", " << nameNextIterSize << endl;
                curNameSize += nameNextIterSize;
                if (tmpNameSize <= curNameSize)
                {
                    tmpNameSize = curNameSize;
                    tmpChunkSet = tmpSet;
                    tmpName = tmpPtr->name;
                    //cout << "(tmpPtr->name):" << tmpPtr->name << ", tmpNameSize: " << tmpNameSize << endl;
                    changeNames = 1;
                }
                freeNodes(nListNextIter);
            }
            ptr = ptr->next;
        }
        vector<string> popGateVec;
        popNodes(nListSecond, tmpChunkSet, popGateVec);
        chunkSet.insert(tmpChunkSet.begin(), tmpChunkSet.end());
        maxNodeSize = chunkSet.size();
        ptr = nListSecond; // back to head
    }
    freeNodes(nListSecond);
}

int main(int argc, char *argv[]) {
    vector<string> gateVec;
    readCircuitToVec(argv[1], gateVec); // read from file and gen the uni name
    int mode = 0;
    if(argc >= 3)
        mode = stoi(argv[2]);
    size_t qubits = 40;
    vector<Node> dags (qubits);
    vector<Node*> tails (qubits);
    for(int i = 0; i < qubits; i++)
        tails[i] = &(dags[i]);
    Node *nList = new Node();
    Node *head = nList; // backup the head;
    setDagFromVec(gateVec, tails, nList);
    //printNodes(nList);
    set<int> chunkSet {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    //set<int> chunkSet {0, 1, 2};
    set<int> chunkSetSwapOut;

    int chunkSize = chunkSet.size();
    int iter = 0;
    int nameSize = 0;
    while(nList->next) {
        chunkSetSwapOut = chunkSet;
        if(iter != 0) {
            if(mode == 0)
                findChunkbyMaxNodePrinciple(chunkSize, nList, chunkSet, nameSize);
            else
                findChunkbyMaxNodeWithNextPrinciple(chunkSize, nList, chunkSet, dags, tails);
        }
        vector<string> popGateVec;
        popNodes(nList, chunkSet, popGateVec);
        setDagFromNodes(dags, tails, nList);
        //printNodes(nList);
        iter++;

        set<int> intersect, uChunkSet;
        set_intersection(chunkSet.begin(), chunkSet.end(), chunkSetSwapOut.begin(), chunkSetSwapOut.end(),
            std::inserter(intersect, intersect.begin()));
        set_union(chunkSet.begin(), chunkSet.end(), chunkSetSwapOut.begin(), chunkSetSwapOut.end(),
            std::inserter(uChunkSet, uChunkSet.begin()));
        if(iter != 1) {
            vector<int> sIn, sOut;
            set_difference(chunkSet.begin(), chunkSet.end(), intersect.begin(), intersect.end(),
                std::inserter(sIn, sIn.begin()));
            set_difference(chunkSetSwapOut.begin(), chunkSetSwapOut.end(), intersect.begin(), intersect.end(),
                std::inserter(sOut, sOut.begin()));
            
            int cntOut = 0;
            for(auto &c: sOut) {
                uChunkSet.erase(c);
                if (cntOut >= sIn.size() - 1)
                    break;
                cntOut++;
            }

            cout << ", 當前存在 in chunk: ";
            for(auto &c: chunkSetSwapOut)
                cout << setw(2) << c << ", ";
            cout << endl;

            cout << ", 下輪存在 in chunk: ";
            for(auto &c: uChunkSet)
                cout << setw(2) << c << ", ";
            cout << endl;


           
            /*
            cout << ", 下輪必有 in chunk: ";
            for(auto &c: chunkSet)
                cout << c << ", ";
            cout << endl;
            cout << ", 交集 in chunk: ";
            for(auto &c: intersect)
                cout << c << ", ";
            cout << endl;
 
            cout << ", Swap in chunk: ";
            for(auto &c: sIn) {
                cout << c << ", ";
            }
            cout << endl;
            
            cout << ", Swap out chunk: ";
            int cntOut = 0;
            for(auto &c: sOut) {
                uChunkSet.erase(c);
                cout << c << ", ";
                if (cntOut >= sIn.size() - 1)
                    break;
                cntOut++;
                //cout << "cntOut: " << cntOut << ", len: " << sIn.size() << endl;
            }
            cout << endl;
            */
            chunkSet = uChunkSet;
            /*
            cout << ", Swap in chunk: ";
            for(auto &c: chunkSet)
                cout << c << " ";

            cout << ", Swap out chunk: ";
            for(auto &c: chunkSetSwapOut)
                cout << c << " ";
            */
        }
        //cout << endl;
        cout << "[iter " << iter << "]" << " ";
        for(int i = 0; i < popGateVec.size(); i++)
            cout << popGateVec[i] << " ";
        nList = head;
        cout << endl;
    }

    //cout << "[ ori_gate ]" << endl;
    //for(int i = 0; i < gateVec.size(); i++)
    //    cout << gateVec[i] << " " << endl;
}