#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include "LinearAlgebra.h"
#include <iterator>

using namespace std;

template<typename T>
void Print_matrix(const vector<vector<T>>& matrix) {
    int m = matrix.size();
    cout << "Printing Matrix ...  " << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

vector<vector<double>> Data_to_Perms(const string& filename){
    int NumOfParticles = 0;
    vector<vector<double>> pmatrix;
    ifstream inputFile(filename);
    
    if (!inputFile.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        return pmatrix;
    }

    // Finding the total number of particles:
    string line;
    while(getline(inputFile, line)){
        int number;
        istringstream iss(line);
        vector<string> tokens(istream_iterator<string>{iss},istream_iterator<string> ());
        for(int i=1; i < tokens.size()-1; i++){
            if(i == 1){
                number = stoi(tokens[i].substr(1,tokens[0].size()));
                if( number > NumOfParticles){
                    NumOfParticles = number;
                }
            }
            else{
                number = stoi(tokens[i]);
                if(stoi(tokens[i]) > NumOfParticles){
                    NumOfParticles = stoi(tokens[i]);
                }
            }
        }
    }
    inputFile.close();


    // Now building the permutation matrices
    inputFile.open(filename);
    while(getline(inputFile, line)){
        vector<double> perm(NumOfParticles+1 , 0);
        int number;
        istringstream iss(line);
        vector<string> tokens(istream_iterator<string>{iss},istream_iterator<string> ());
        for(int i=1; i < tokens.size()-1; i++){
            // Determining whether dagger or not:
            double creation = 1.0;
            size_t pos = -1;
            pos = tokens[i].rfind('^');
            if(pos != -1){
                creation = -1.0;
            }
            if(i == 1){
                number = stoi(tokens[i].substr(1,tokens[0].size()));
                perm[number] = creation;
            }
            else{
                number = stoi(tokens[i]);
                perm[number] = creation;
            }
        }
        pmatrix.push_back(perm);
    }
    inputFile.close();

    return Transpose(pmatrix);
}

int main(int argc , char* argv[]) {
    vector<vector<double>> PermMatrix;

    string filename(argv[1]);
    PermMatrix = Data_to_Perms(filename);

    vector<vector<int>> Cycles = Nullspace(PermMatrix);

    Print_matrix(Cycles);

    return 0;
}