#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <algorithm>
#include "LinearAlgebra.h"
#include <iterator>

using namespace std;

typedef vector<vector<double>> CrAnOps;
typedef vector<vector<int>> NumberOps; 

struct DPdataElec{
    vector<vector<complex<double>>> Coeffs;
    vector<NumberOps> Diagonals;
    CrAnOps Permutations;
    vector<complex<double>> Coeffs0;
    vector<vector<int>> D0;  
};

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

double Sort_permute(vector<int> Vec){
    vector<int> VecSorted = Vec;
    int n = VecSorted.size();
    double permfactor = 1.0;
    sort(VecSorted.begin() , VecSorted.end());

    for(int i = 0 ; i < n-1; ++i){
        for(int j=i+1; j < n; ++j){
            if(Vec[i] > VecSorted[j]){
                permfactor = permfactor*(-1.0);
            }
        }
    }
    return permfactor;
}

bool Diag_compare(vector<int> d1 , vector<int> d2){
    int s1 = d1.size() , s2 = d2.size();
    vector<int> D2 = d2;
    bool DiagMatched;
    if(s1 != s2){
        return false;
    }
    else{
        for(int i = 0; i < s1; i++){
            DiagMatched = false;
            for(int j = 0; j < D2.size(); j++){
                if(abs(d1[i]- D2[j]) < 1E-6){
                    D2.erase(D2.begin() + j);
                    DiagMatched = true;
                    break;
                }
            }
            if(!DiagMatched){
                return false;
            }
        }
        return true;
    }
}

vector<double> Perm_sub(vector<double> p1 , vector<double> p2 ){
    vector<double> output = p1;
    for(int i = 0; i < p2.size(); ++i){
        output[i] -= p2[i];
    }
    return output;
}

bool All_zeros(vector<double> v){
    for(int i = 0; i < v.size() ; ++i){
        if(abs(v[i]) > 1E-6){
            return false;
        }
    }
    return true;
}

bool Perm_compare(vector<double> p1 , vector<double> p2){
    int s1 = p1.size() , s2 = p2.size();
    vector<double> diff = Perm_sub(p1, p2);
    if(All_zeros(diff)){
        return true;
    }
    else{
        return false;
    }
}

pair<bool , int> Diag_found(vector<int> v , NumberOps Vs){
    for(int i = 0; i < Vs.size(); ++i){
        if(Diag_compare(v , Vs[i])){
            return make_pair(true , i);
        }
    }
    return make_pair(false , -1);    
}

pair<bool , int> Perm_found(vector<double> v , CrAnOps Vs){
    for(int i = 0; i < Vs.size(); ++i){
        if(Perm_compare(v , Vs[i])){
            return make_pair(true , i);
        }
    }
    return make_pair(false , -1); 
}

complex<double> Read_coeff(const string complexStr){
    int plusPos = -1 , minPos = -1 , iPos = -1;
    plusPos = complexStr.find('+');
    minPos = complexStr.find('-');
    iPos = complexStr.find('j');

    string realPartStr , imagPartStr;

    if(iPos != -1){
        if(plusPos != -1){
            realPartStr = complexStr.substr(0, plusPos);
            imagPartStr = complexStr.substr(plusPos + 1, iPos - plusPos);
        }
        else if(minPos > 0){
            realPartStr = complexStr.substr(0, minPos);
            imagPartStr = complexStr.substr(minPos, iPos - minPos);
        }
        else{
            realPartStr = '0';
            imagPartStr = complexStr.substr(0 , iPos);
        }
    }
    else{
        realPartStr = complexStr;
        imagPartStr = '0';
    }

    // Convert the substrings to double
    double realPart = stod(realPartStr);
    double imagPart = stod(imagPartStr);

    // Create a complex<double> object
    complex<double> complexNumber(realPart, imagPart);
    return complexNumber;
}   

DPdataElec Data_to_Perms(const string& filename){
    int NumOfParticles = 0;
    DPdataElec Data;
    CrAnOps pmatrix;
    ifstream inputFile(filename);
    
    if (!inputFile.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        return Data;
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
        vector<double> perm(NumOfParticles+1 , 0.0);
        vector<int> diagonal;
        int number;
        istringstream iss(line);
        vector<int> CrString = {}, AnString = {};
        vector<string> tokens(istream_iterator<string>{iss},istream_iterator<string> ());
        complex<double> coeff = Read_coeff(tokens[0]);
        for(int i=1; i < tokens.size()-1; i++){
            // Determining whether dagger or not:
            double creation = 1.0;
            size_t pos = -1;
            pos = tokens[i].rfind('^');
            
            if(i == 1){
                number = stoi(tokens[i].substr(1,tokens[0].size()));
            }
            else{
                number = stoi(tokens[i]);
            }
            if(pos == -1){
                creation = -1.0;
                AnString.push_back(number);
            }
            else{
                coeff = coeff * pow(-1.0 , AnString.size());
                CrString.push_back(number);
            }
            perm[number] += creation;
            if(abs(perm[number]) < 1E-6){
                diagonal.push_back(number);
            }
        }
        coeff = coeff * Sort_permute(CrString) * Sort_permute(AnString);

        // Determine if the permutation is found and if it is, check if the new diagonal is a part of the D_i: if not, add it to the list with the corresponding coefficient
        // if it already is in the list of diagonals just add the coefficients!

        if(All_zeros(perm)){
            Data.Coeffs0 ;
            Data.D0 ;
            pair<bool , int> Dfound = Diag_found(diagonal , Data.D0);
            if(Dfound.first){
                Data.Coeffs0[Dfound.second] += coeff;
            }
            else{
                Data.Coeffs0.push_back(coeff);
                Data.D0.push_back(diagonal);
            }
        }
        else{
            pair<bool , int> Pfound = Perm_found(perm , Data.Permutations);
            if(Pfound.first){
                pair<bool , int> Dfound = Diag_found(diagonal , Data.Diagonals[Pfound.second]);
                if(Dfound.first){
                    Data.Coeffs[Pfound.second][Dfound.second] += coeff;
                }
                else{
                    Data.Coeffs[Pfound.second].push_back(coeff);
                    Data.Diagonals[Pfound.second].push_back(diagonal);
                }
            }
            else{
                Data.Coeffs.push_back({coeff});
                Data.Diagonals.push_back({diagonal});
                Data.Permutations.push_back(perm);
            } 
        }
    }
    inputFile.close();
    return Data;
}

int main(int argc , char* argv[]) {
    vector<vector<double>> PermMatrix;
    vector<vector<complex<double>>> Cs;
    DPdataElec ElecData;

    string filename(argv[1]);
    ElecData = Data_to_Perms(filename);

    PermMatrix = Transpose(ElecData.Permutations);
    Cs = ElecData.Coeffs;
    vector<vector<int>> Cycles = Nullspace(PermMatrix);

    cout << endl;
    cout << "The permutations are: "  << endl;
    Print_matrix(PermMatrix);
    cout << endl;

    cout << "The coefficients are: " << endl;
    Print_matrix(Cs);
    cout << endl;


    cout << "The cycles are: " << endl;
    Print_matrix(Cycles);

    return 0;
}