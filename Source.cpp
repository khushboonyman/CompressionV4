#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <set>
#include <math.h>
#include <chrono>
#include <vector>
#include<algorithm>

using namespace std;
//GLOBAL VARIABLES THAT NEED TO BE CHANGED ACCORDINGLY
int version = 4;
int limit = 10;
int recursiveLimit = 1000;
int runLimit = 10000;
//int runLimit = 1000000;
string location_main = "C:\\Users\\Bruger\\Desktop\\books\\THESIS start aug 3\\datasets\\";
//file name here
//string fileName = "genome.fa";
//string fileName = "Gen178.fa";
//string fileName = "embl50.h178.fa";
//FILE TO BE LOGGED
//string fileName = "my_complete_genome.txt";
//string fileName = "genome_selfmade_100.txt";
//string fileName = "main_test_file_1.fsa";
string fileName = "dna";
string* dnaArray;
unordered_map<string, vector<int>> fingerPrints;
unordered_map<char, int> singleChar;
int memory = 0;
string relativeString = "";
int relativeSize;
int numberOfStrings = 0;
int memoryVar = 0;
int memoryOld = 0;
int base = 0;
int mb = 1024 * 1024;
int kb = 1024;

vector<int>* indexRelative;
vector<int>* indexCString;
vector<int>* originalIndex = new vector<int>;
vector<int>* pointerIndex = new vector<int>;
int recursive = 0;
string extra = "";

//FIND HOW MANY STRINGS ARE THERE
void findSize(string &location) {
    ifstream myfile;
    myfile.open(location);
    cout << "dna file opened" << endl;
    string x;

    while (myfile >> x) {
        if (x[0] == '>')
            numberOfStrings += 1;
    }

    myfile.close();
}

//CREATE ARRAY OF DNA STRINGS FROM THE FILE
void readDna(string &location) {
    ifstream myfile;
    myfile.open(location);
    int size = -1;
    string temp, x;

    while (myfile >> x) {
        if (x[0] == '>') {
            cout << "index "<<size + 1 << endl;
            if (size > -1) {
                dnaArray[size] = temp;
            }
            temp = "";
            size ++;
        }
        else {
            temp += x;
        }
            
    }
    dnaArray[size] = temp;
    myfile.close();
    cout << endl << "dna file closed" << endl;
}

//FIND HOW MANY STRINGS ARE THERE
void findSizePizzaChilli(string& location) {
    ifstream myfile;
    myfile.open(location);
    cout << "dna file opened" << endl;
    string x;

    while (myfile >> x) {
        numberOfStrings += 1;
    }

    myfile.close();
}

//CREATE ARRAY OF DNA STRINGS FROM THE FILE
void readDnaPizzaChilli(string& location) {
    ifstream myfile;
    myfile.open(location);
    int size = 0;
    string x;
    while (myfile >> x) {
        dnaArray[size] = x;
        //TEST ONLY!!!
        relativeString += x;
        size++;
    }
    myfile.close();
    cout << endl << "dna file closed and size " << size << endl;
}


void findBase() {
    set<char> setOfChars;
    
    for (int i = 0; i < numberOfStrings; i++) {
        for (int j = 0; j < dnaArray[i].size(); j++) {
            if (setOfChars.find(dnaArray[i][j]) == setOfChars.end()) {
                cout << "char not found " << dnaArray[i][j] << endl;
                setOfChars.insert(dnaArray[i][j]);
            }
        }
    }

    base = setOfChars.size();
    cout << " base is " << base << endl;
}

//WRITE A LOG TO THE FILE, WHICH WILL BE USED TO PLOT TIME AND SPACE USED BY COMPRESSION TECHNIQUE
void writeLog(string &location, string &fileName, int &version, int &memoryVar, int &time) {
    fstream myfile;
    myfile.open(location, fstream::app);
    cout << "log file opened" << endl;
    myfile << "\n";
    myfile << fileName << ";" << version << ";" << memoryVar << ";" << time;
    myfile.close();
    cout << "log file closed" << endl;
}

//CREATING MAPS <FINGERPRINT,INDEXOFSTRING> and <FINGERPRINT,NUMBEROFOCCURRENCES>
void updateFingerPrints(unordered_map<string, vector<bool>>& fingerPrintsIndex, string &currentString, int &i, int &j) {
    string fingerPrint;
    fingerPrint.reserve(limit);
    fingerPrint.append(currentString, j, limit);
    unordered_map<string, vector<bool>>::iterator fpiCheck = fingerPrintsIndex.find(fingerPrint);
    if (fpiCheck != fingerPrintsIndex.end()) {
        vector<bool> stringIndex = fpiCheck->second;
        if (!stringIndex[i]) {
            stringIndex[i] = true;
            stringIndex[numberOfStrings] = true;
            fingerPrintsIndex[fingerPrint] = stringIndex;
        }
    }
    else {
        vector<bool> stringIndex(numberOfStrings+1);
        for (int x = 0; x < numberOfStrings+1; x++) {
            stringIndex[x] = false;
        }
        stringIndex[i] = true;
        fingerPrintsIndex[fingerPrint] = stringIndex;
    }
}

int findNextRelativeString(vector<int>& common , int& relativeFirstStr, int& minCommon) {
    cout << " finding next relative ";
    int toReturn = 0;

    for (int x = 0; x < numberOfStrings; x++) {
        if (x != relativeFirstStr && common[x] < minCommon) {
            minCommon = common[x];
            toReturn = x;
        }
    }

    cout << "minimum common fingerprint " << minCommon << endl;

    return toReturn;
}

int findFirstRelativeString(unordered_map<string, vector<bool>>& fingerPrintsIndex, vector<vector<int> >& common, int& maxCommon) {

    cout << "finding first part of relative string, size of all fingerPrints  " <<fingerPrintsIndex.size() << endl;
    int relativeFirstStr = 0, i = 0;
    
    unordered_map<string, vector<bool>>::iterator fpiLoop = fingerPrintsIndex.begin();

    cout << "creating two dimensional array " << endl;
    while (fpiLoop != fingerPrintsIndex.end()) {
        if (i % 1000 == 0) {
            cout <<"processing "<< i << endl;
        }
        vector<bool> stringIndex = fpiLoop->second;
        for (int x = 0; x < numberOfStrings; x++) {
            if (stringIndex[x]) {
                for (int y = 0; y < numberOfStrings; y++) {
                    if (x != y && stringIndex[y]) {
                        common[x][y]++;
                    }
                }
            }
        }
        i++;
        fpiLoop++;
    }

    cout << "finding maximum common " << endl;

    for (int x = 0; x < numberOfStrings; x++) {
        int maxX = 0;
        for (int y = 0; y < numberOfStrings; y++) {
            if (x != y) {
                maxX += common[x][y];
            }
        }
        if (maxX > maxCommon) {
            maxCommon = maxX;
            relativeFirstStr = x;
        }
    }

    cout << "max common " << maxCommon << endl;

    return relativeFirstStr;
}

//FIND A GOOD RELATIVE STRING
string findRelativeString() {
    string toReturn = "";
    cout << "finding relative string" << endl;

    int totalSize=0, relativeFirstStr=-1, relativeSecondStr=-1, maxFingerPrint=0, maxCommon=0, chosenStringSize=0;
    unordered_map<string, vector<bool>> fingerPrintsIndex ;
    unordered_map<string, vector<bool>> fingerPrintsRemoveFirst;

    cout << "creating vector of fingerprints" << endl;

    for (int i = 0; i < numberOfStrings; i++) {
        cout << "first loop " << i << " size of string " << dnaArray[i].size() ;
        totalSize += dnaArray[i].size();
        for (int j = 0; j <= dnaArray[i].size() - limit; j++) {
            updateFingerPrints(fingerPrintsIndex, dnaArray[i], i, j);
        }        
        cout << " size of fingerprints " << fingerPrintsIndex.size() << endl;
    }

    cout << "intializing number of occurrences" << endl;

    vector<int> numberOfOccurrences(numberOfStrings);
    for (int x = 0; x < numberOfStrings; x++) {
        numberOfOccurrences[x] = 0;
    }

    cout << "creating array numberOfOccurrences for each string" << endl;

    int z = 0;
    //GO THROUGH EACH OF THE FINGERPRINTS AND THEIR CORRESPONDING VECTORS AND MARK THE NUMBER OF OCCURRENCES (HOW MANY FINGERPRINTS THE STRING HAS, THAT OTHERS ALSO HAVE)
    unordered_map<string, vector<bool>>::iterator fpiLoop = fingerPrintsIndex.begin();
    while (fpiLoop != fingerPrintsIndex.end()) {
        if (z % 10000 == 0) {
            cout << "in the loop : value of z " << z << endl;
        }
        vector<bool> stringIndex = fpiLoop->second;
        for (int x = 0; x < numberOfStrings; x++) {
            if (stringIndex[x] && stringIndex[numberOfStrings]) {
                numberOfOccurrences[x] += 1;
            }
        }
        fpiLoop++;
        z++;
    }

    cout << "checking which string has most number of occurrences with others" << endl;

    //GETTING THE FIRST PART OF RELATIVE STRING
    for (int i = 0; i < numberOfStrings; i++) {
        if (numberOfOccurrences[i] > maxFingerPrint) {
            maxFingerPrint = numberOfOccurrences[i];
            relativeFirstStr = i;
        }
    }

    cout << "moving fingerprints that are not in first part of string" << endl;
    
    z = 0;
    fpiLoop = fingerPrintsIndex.begin();
    while (fpiLoop != fingerPrintsIndex.end()) {
        if (z % 10000 == 0) {
            cout << "in the loop : size of fingerprints" << fingerPrintsRemoveFirst.size() <<endl ;
        }
        vector<bool> stringIndex = fpiLoop->second;
        if (!(stringIndex[relativeFirstStr]) && stringIndex[numberOfStrings]) {
            stringIndex[numberOfStrings] = false;
            bool count = false;
            for (int x = 0; x < numberOfStrings; x++) {
                if (stringIndex[x]) {
                    if (!count) {
                        count = true;
                    }
                    else {
                        stringIndex[numberOfStrings] = true;
                        break;
                    }
                }
            }
            if(stringIndex[numberOfStrings])
                fingerPrintsRemoveFirst[fpiLoop->first] = stringIndex; 
        }
        z++;
        fpiLoop++;
    }

    cout << " size of fingerprints " << fingerPrintsRemoveFirst.size() << endl;

    cout << "creating array numberOfOccurrences for each string" << endl;
    
    for (int x = 0; x < numberOfStrings; x++) {
        numberOfOccurrences[x] = 0;
    }

    //GO THROUGH EACH OF THE FINGERPRINTS AND THEIR CORRESPONDING VECTORS AND MARK THE NUMBER OF OCCURRENCES (EXCLUDING FIRST RELATIVE STRING)
    fpiLoop = fingerPrintsRemoveFirst.begin();
    z = 0;
    while (fpiLoop != fingerPrintsRemoveFirst.end()) {
        if (z % 10000 == 0) {
            cout << "in the loop : z is " << z << endl;
        }
        vector<bool> stringIndex = fpiLoop->second;
        for (int x = 0; x < numberOfStrings; x++) {
            if (stringIndex[x] && stringIndex[numberOfStrings]) {  //KHUSH : THIS CHECK NEEDS TO BE CORRECTED ABOVE
                numberOfOccurrences[x] += 1;
            }
        }
        fpiLoop++;
        z++;
    }

    cout << "checking which string has most number of occurrences excluding the first one" << endl;

    maxFingerPrint = 0;
    //GETTING THE SECOND PART OF RELATIVE STRING
    for (int i = 0; i < numberOfStrings; i++) {
        if (numberOfOccurrences[i] > maxFingerPrint) {
            maxFingerPrint = numberOfOccurrences[i];
            relativeSecondStr = i;
        }
    }

    toReturn = dnaArray[relativeSecondStr] + dnaArray[relativeFirstStr];
    
    cout << " size of string array " << totalSize;
    cout << " first relative string " << relativeFirstStr << " second relative string " << relativeSecondStr <<endl;

    return toReturn;
}

//FIND IF THE FINGERPRINT EXISTS AND WHICH INDICES ARE THERE
vector<int> findFingerPrint(string& checkFingerPrint) {
    unordered_map<string, vector<int>>::iterator itCheck = fingerPrints.find(checkFingerPrint);
    vector<int> returnVector;
    if (itCheck != fingerPrints.end())
        returnVector = itCheck->second;
    return returnVector;
}

//PRINT LIST OF SINGLE CHARACTERS AND THEIR FIRST POSITION IN THE STRING
void printSingleChar() {
    unordered_map<char, int>::iterator itChar = singleChar.begin();
    while (itChar != singleChar.end()) {
        cout << "char : " << itChar->first << " index : " << itChar->second << endl;
        itChar++;
    }
}

//FIND IF THE SINGLE CHARACTER EXISTS AND WHAT IS THE INDEX
int findSingleChar(char& checkChar) {
    unordered_map<char, int>::iterator itCheck = singleChar.find(checkChar);
    if (itCheck != singleChar.end())
        return itCheck->second;
    return -1;
}

//SET UP THE DATASTRUCTURE CONTAINING FINGERPRINTS AND SINGLE CHARACTERS ALONG WITH THE INDICES
void setFingerPrintSingleChar() {
    int i;
    for (i = 0; i <= relativeSize - limit; i++) {
        string fingerPrint;
        fingerPrint.reserve(limit);
        fingerPrint.append(relativeString, i, limit);
        char single = relativeString[i];
        unordered_map<string, vector<int>>::const_iterator it = fingerPrints.find(fingerPrint);
        unordered_map<char, int>::const_iterator itC = singleChar.find(single);
        if (itC == singleChar.end())
            singleChar[single] = i;
        if (it == fingerPrints.end()) {
            vector<int> newVector;
            newVector.push_back(i);
            fingerPrints[fingerPrint] = newVector;
        }
        else {
            vector<int> existingVector = it->second;
            existingVector.push_back(i);
            fingerPrints[fingerPrint] = existingVector;
        }
    }

    while (i < relativeSize) {
        char single = relativeString[i];
        unordered_map<char, int>::const_iterator itC = singleChar.find(single);
        if (itC == singleChar.end())
            singleChar[single] = i;
        i++;
    }
}

//ADD A NEW CHARACTER TO THE RELATIVE STRING, IF IT WASN'T EXISTING INITIALLY WHEN THE RELATIVE STRING WAS JUST THE FIRST STRING FROM THE FILE 
int expandRelative(char &charToAdd) {
    relativeString += charToAdd;
    relativeSize += 1;
    singleChar[charToAdd] = relativeSize - 1;
    return relativeSize - 1;
}

//process compression based on if fingerPrint found or not 
void findFingerPrint(int &start, string &toCompress, vector<int> &indexRelativeElement, vector<int> &indexCStringElement) {
    vector<int> indices;
    string checkFingerPrint;
    checkFingerPrint.reserve(limit);
    checkFingerPrint.append(toCompress, start, limit);
    indices = findFingerPrint(checkFingerPrint);
    //limit size substring fingerprint not found
    if (indices.size() == 0) {
        int index = findSingleChar(toCompress[start]);
        if (index == -1) {
            cout << "char not found : " << toCompress[start] << endl;
            index = expandRelative(toCompress[start]);
        }
        indexRelativeElement.push_back(index);
        indexCStringElement.push_back(start);
        start++;
    }
    //fingerprint(s) found, so now we find the longest substring that matches through the fingerprint(s)
    else {
        int max_length = limit;
        int max_index = -1;
        for (vector<int>::iterator itV = indices.begin(); itV != indices.end(); itV++) {
            if (max_index == -1) {
                max_index = *itV;
            }
            int currIndexString = start + limit;
            int currIndexRelativeString = *itV + limit;
            int length = limit;
            while (true) {
                if (currIndexString >= toCompress.size() || currIndexRelativeString >= relativeSize) {
                    break;
                }
                char stringChar = toCompress[currIndexString];
                char relativeChar = relativeString[currIndexRelativeString];
                if (stringChar == relativeChar) {
                    currIndexString++;
                    currIndexRelativeString++;
                    length++;
                }
                else {
                    break;
                }
            }
            if (length > max_length) {
                max_length = length;
                max_index = *itV;
            }
        }
        indexRelativeElement.push_back(max_index);
        indexCStringElement.push_back(start);
        start += max_length;
    }
}

//COMPRESS ONE STRING AT A TIME
void compress(string& toCompress, vector<int>& indexRelativeElement, vector<int>& indexCStringElement) {
    int start = 0;
    while (start <= toCompress.size() - limit) {
        findFingerPrint(start, toCompress, indexRelativeElement, indexCStringElement);
    }

    while (start < toCompress.size()) {
        int index = findSingleChar(toCompress[start]);
        if (index == -1) {
            cout << "char not found : " << toCompress[start] << endl;
            index = expandRelative(toCompress[start]);
        }
        indexRelativeElement.push_back(index);
        indexCStringElement.push_back(start);
        start++;
    }
}

//FIND LOCATION OF THE INDEX OF SEARCHED STRING IN THE RELATIVE STRING
int findLocation(vector<int>& indexCStringElement, int& charIndex) {
    int count = 0;
    int first = 0;
    int last = ((charIndex < (int)indexCStringElement.size()) ? charIndex : (indexCStringElement.size() - 1));
    //int last = indexCStringElement.size() - 1;
    int mid = last / 2;
    int indexCStringCurrent, indexCStringNext, indexCStringNextNext;
    while (true) {
        count++;
        if (count > recursiveLimit) {
            cout << "first " << first << " mid " << mid << " last " << last << endl;
        }
        if (count > 1010) {
            cout << "limit reached in findLocation" << endl;
            break;
        }
        indexCStringCurrent = indexCStringElement[mid];

        if (mid < (int)(indexCStringElement.size() - 1)) {
            indexCStringNext = indexCStringElement[mid + 1];
            if (indexCStringCurrent <= charIndex && indexCStringNext > charIndex) {
                break;
            }
            if (mid < (int)(indexCStringElement.size() - 2)) {
                indexCStringNextNext = indexCStringElement[mid + 2];
                if (indexCStringNext <= charIndex && indexCStringNextNext > charIndex) {
                    mid++;
                    break;
                }
            }
        }
        else {
            if (indexCStringCurrent <= charIndex) {
                break;
            }
        }

        if (mid == first) {
            mid = last;
        }
        else {
            if (indexCStringCurrent > charIndex) {
                last = mid;
            }
            else {
                first = mid;
            }

            mid = (first + last) / 2;
        }
    }
    return mid;
}

//FIND CHARACTER FROM COMPRESSED REFERENCE STRING
char findCharIndexRefString(int& index) {
    recursive++;
    if (recursive > recursiveLimit) {
        cout << "something went wrong in recursion " << endl;
        return '@';
    }

    int first = 0;
    int last = (index <= extra.size() ? index : extra.size() - 1);
    int mid = last / 2;
    int count = 0;
    while (true) {
        if (originalIndex->at(mid) == index) {
            if (originalIndex->at(mid) == pointerIndex->at(mid)) {
                return extra[mid];
            }
            else {
                return findCharIndexRefString(pointerIndex->at(mid));
            }
        }

        if (mid == first) {
            if (originalIndex->at(last) == pointerIndex->at(last) && pointerIndex->at(last) == index) {
                return extra[last];
            }
        }
        if (mid != extra.size() - 1) {
            int next = mid + 1;
            if (originalIndex->at(mid) < index && originalIndex->at(next) > index) {
                if (originalIndex->at(next) == index + 1)
                    return extra[mid];
                else {
                    int callIndex = pointerIndex->at(mid) + index - originalIndex->at(mid);
                    return findCharIndexRefString(callIndex);
                }
            }
            else {
                if (mid == first) {
                    mid = last;
                }
                else {
                    if (originalIndex->at(mid) > index)
                        last = mid;
                    else
                        first = mid;
                    mid = (first + last) / 2;
                }
            }
        }
        else {
            if (index == extra.size() - 1)
                return extra[mid];
            else {
                int callIndex = pointerIndex->at(mid) + index - originalIndex->at(mid);
                return findCharIndexRefString(callIndex);
            }
        }
        count++;
        if (count > 50) {
            cout << "something went wrong again " << endl;
            break;
        }
    }
    return '$';
}

//FIND WHICH CHARACTER IS PRESENT ON THE STRING INDEX
char findCharacter(vector<int>& indexRelativeElement, vector<int>& indexCStringElement, int& charIndex) {
    int indexFound = findLocation(indexCStringElement, charIndex);
    int distance = charIndex - indexCStringElement[indexFound];
    int finalIndex = indexRelativeElement[indexFound] + distance;
    recursive = 0;
    return findCharIndexRefString(finalIndex);
}

//PROCESS MANY RANDOM REQUESTS OF FINDING A CHARACTER ON AN INDEX FROM A STRING
auto processRandomRequests(int* sizes) {
    cout << " processing " << runLimit << " requests " << endl;
    int counter = 0;
    int stringIndex, charIndex;
    //measuring time start
    auto start = chrono::high_resolution_clock::now();

    while (counter < runLimit) {
        if (counter % 100 == 0) {
            cout << counter << endl;
        }
        stringIndex = rand() % (numberOfStrings - 1);
        charIndex = rand() % (sizes[stringIndex] - 1);

        vector<int> indexRelativeElement = indexRelative[stringIndex];
        vector<int> indexCStringElement = indexCString[stringIndex];
        //this function finds the character from the compressed datastructure
        char charFound = findCharacter(indexRelativeElement, indexCStringElement, charIndex);
        /* TEST */
        /*if (charFound != dnaArray[stringIndex][charIndex]) {
            cout << "some issue fetching character " << stringIndex << " " << charIndex << endl;
            break;
        }*/
        counter++;
    }
    //measuring time end
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    return duration;
}

//FIND INDEX OF THE SUBSTRING FROM HASHMAP. IF NOT FOUND, THEN RETURN -1
int findIndexSubString(string& current, unordered_map<string, int>& subStrings) {
    unordered_map<string, int>::iterator itCheck = subStrings.find(current);
    if (itCheck != subStrings.end()) {
        return itCheck->second;
    }
    return -1;
}

//UPDATE INT, INT, CHAR VECTORS WITH ORIGINAL INDEX, POINTER INDEX AND EXTRA CHARACTER AND ALSO UPDATE THE HASHMAP TO STORE SUBSTRING TO INDEX MAPPING
void updateVector(int& i, char& additional, int& prevIndex, string& currentString, unordered_map<string, int>& subStrings, int& memory) {
    originalIndex->push_back(i);
    extra += additional;
    if (prevIndex == -1) {
        pointerIndex->push_back(i);
    }
    else {
        pointerIndex->push_back(prevIndex);
    }
    memory += 9;
    subStrings[currentString] = i;
}

//MAIN LOGIC TO COMPRESS THE REFERENCE STRING
void compressReference() {
    unordered_map<string, int> subStrings;
    string currentString = "";
    memory = sizeof(originalIndex) + sizeof(pointerIndex) + sizeof(extra) ;
    cout << "initial memory when compressing reference string " << memory << endl;

    int currentIndex = -1, prevIndex = -1, i;
    int baseIndex = 0;

    for (i = 0; i < relativeString.size(); i++) {
        prevIndex = currentIndex;
        currentString += relativeString[i];
        currentIndex = findIndexSubString(currentString, subStrings);

        if (currentIndex == -1) {
            updateVector(baseIndex, relativeString[i], prevIndex, currentString, subStrings, memory);
            currentString = "";
            baseIndex = i + 1;
        }
    }

    char additional = '$';
    currentString += '$';
    updateVector(baseIndex, additional, currentIndex, currentString, subStrings, memory);
    
    cout << "memory used by compression of main string : " << memory << "original size " << relativeString.size() << endl;
}


int main() {
    cout << "PROGRAM STARTING WITH LIMIT " << limit <<" FILE NAME "<< fileName << " !!!!!"<< endl;

    int i;
    string location = location_main + fileName;

    if (fileName.substr(0, 3) == "dna" || fileName.substr(0, 3) == "pro") {
        findSizePizzaChilli(location);
    }
    else {
        findSize(location);
    }
    
    dnaArray = new string[numberOfStrings];
    cout << "NUMBER OF STRINGS " << numberOfStrings << endl;
    
    if (fileName.substr(0, 3) == "dna" || fileName.substr(0, 3) == "pro") {
        readDnaPizzaChilli(location);
    }
    else {
        readDna(location);
    }
    
    //findBase();    //in case we are making our own hash function

    int* sizes = new int[numberOfStrings];

    memoryVar += (sizeof(sizes) + (4 * numberOfStrings));  //adding size of pointer for size array and each element of size array

    for (i = 0; i < numberOfStrings; i++) {
        sizes[i] = dnaArray[i].size();
        memoryOld += sizes[i]; //adding length of each string
    }

    cout << "DNA ARRAY READ !!!" << endl;
    
    relativeString = findRelativeString();
    //relativeString = dnaArray[0];
    //ONLY TEST!!!
    //relativeString = dnaArray[0] + dnaArray[14];
    relativeSize = relativeString.size();

    cout << "relative string found and its size is " << relativeSize<< endl;

    setFingerPrintSingleChar();
    printSingleChar();

    indexRelative = new vector<int>[numberOfStrings];
    indexCString = new vector<int>[numberOfStrings];

    memoryVar += (sizeof(indexRelative) + sizeof(indexCString)); //adding pointer size for both index arrays

    //to check how long it takes to compress all the data
    auto start = chrono::high_resolution_clock::now();

    cout << "BEFORE COMPRESSION !!!" << endl;

    for (i = 0; i < numberOfStrings; i++) {
        cout << "compressing " << i << " of length " << dnaArray[i].size() << endl;

        compress(dnaArray[i], indexRelative[i], indexCString[i]);

        memoryVar += 8 * indexRelative[i].size(); //adding space for encoding
        
    }

    //to check how long it takes to compress all the data
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << " duration in millisec " << duration.count() << endl;

    delete[] dnaArray;

    cout << "AFTER COMPRESSION!!!" << endl;

    compressReference();

    auto durationRandomRequests = processRandomRequests(sizes);
    //delete[] dnaArray;
    delete[] sizes;
    delete[] indexRelative;
    delete[] indexCString;
    delete originalIndex;
    delete pointerIndex;

    memoryVar += memory;
    cout << "PROGRAM ENDING!!! " << endl;

    cout << "old memory : " << memoryOld << " compressed memory : " << memoryVar << " compression " << ((float)memoryVar/(float)memoryOld)*100.0 << endl;
    /*//string headers = "FILE_NAME;VERSION;MEMORY;TIME";
    location = location_main + "LOGS.csv";
    writeLog(location, fileName, version, memoryVar, (int)durationRandomRequests.count());*/
}
