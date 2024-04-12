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
    vector<vector<complex<double>>> Coeffs , CoeffsConj;
    vector<NumberOps> Diagonals;
    CrAnOps Permutations;
    vector<complex<double>> Coeffs0;
    vector<vector<int>> D0;  
};

template<typename T>
void Print_matrix(const vector<vector<T>>& matrix) {
    int m = matrix.size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            cout << matrix[i][j] << "  ";
        }
        cout << endl;
    }
}

pair<bool , int> Find_number(int Num , vector<int> NumVec){
    for(int i = 0; i < NumVec.size() ; i++){
        if(NumVec[i] == Num){
            return {true , i};
        }
    }
    return {false , 0};
}

double Sort_permute(vector<int> Vec){
    vector<int> VecSorted = Vec;
    int n = VecSorted.size();
    double permfactor = 1.0;
    sort(VecSorted.begin() , VecSorted.end());

    for(int i = 0 ; i < n-1; ++i){
        for(int j=i; j < n; ++j){
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

vector<double> Perm_add(vector<double> p1 , vector<double> p2 ){
    vector<double> output = p1;
    for(int i = 0; i < p2.size(); ++i){
        output[i] += p2[i];
    }
    return output;
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

vector<double> Perm_conjugate(vector<double> v){
    vector<double> vconj;
    for(int i=0; i < v.size(); i++){
        vconj.push_back(-1.0*v[i]);
    }
    return vconj;
}


complex<double> Read_coeff(const string complexStr){
    int plusPos = -1 , minPos = -1 , iPos = -1 ;
    plusPos = complexStr.find('+');
    minPos = complexStr.find('-');
    iPos = complexStr.find('j');

    // Determine whether the string only starts with a minus sign or there's a minus in the middle of the string
    if(minPos == 0){
        minPos = -1;
        minPos = complexStr.substr(1 , complexStr.size()).find('-');
        minPos += 1;
    }

    string realPartStr , imagPartStr;

    if(iPos != -1){
        if(plusPos != -1){
            realPartStr = complexStr.substr(0, plusPos);
            imagPartStr = complexStr.substr(plusPos + 1, iPos - plusPos);
        }
        else if(minPos == 0){
            realPartStr = '0';
            imagPartStr = complexStr.substr(0 , iPos);
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

pair<int, int> Perm_lengths(vector<double> perm){
    int l = 0 , k = 0;
    pair<int , int> landk;
    for(int i = 0 ; i < perm.size() ; ++i){
        if(abs(perm[i] - 1.0) < 1E-6){
            l++;
        }
        else if(abs(perm[i] + 1.0) < 1E-6){
            k++;
        }
    }
    landk.first = l;
    landk.second = k;
    return landk;
}

vector<vector<complex<double>>> Conjugate_coeffs(vector<vector<complex<double>>> Cs , vector<vector<double>> Perms){
    vector<vector<complex<double>>> CsConj;
    for(int i=0; i < Perms.size() ; i++){
        pair<int, int> cran_lengths = Perm_lengths(Perms[i]);
        double sgn_factor = pow(-1.0 ,cran_lengths.first + cran_lengths.second);
        vector<complex<double>> CsConj_i;
        for(int j=0; j < Cs[i].size() ; ++j){
            // Adding/subtracting Cs[i][j] with its conjugate if the l + k is even/odd
            // l is the number of creation operators and k is the number of annihilation operators
            CsConj_i.push_back(conj(Cs[i][j])*sgn_factor);
        }
        CsConj.push_back(CsConj_i);
    }
    return CsConj;
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
                if(number > NumOfParticles){
                    NumOfParticles = stoi(tokens[i]);
                }
            }
        }
    }
    inputFile.close();

    // Now building the permutation matrices
    inputFile.open(filename);
    while(getline(inputFile, line)){
        vector<double> perm(NumOfParticles+1 , 0.0) , permconj(NumOfParticles+1 , 0.0);
        vector<int> diagonal = {};
        int number;
        istringstream iss(line);
        vector<int> CrString = {}, AnString = {};
        vector<string> tokens(istream_iterator<string>{iss},istream_iterator<string> ());
        complex<double> coeff = Read_coeff(tokens[0]);
        for(int i=1; i < tokens.size()-1; i++){
            // Determining whether dagger or not:
            double creation = 1.0 , Anfactor = 1.0 , Crfactor = 1.0;
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
                pair<bool , int> Anfound = Find_number(number, CrString);
                // Finding a diagonal term and accounting for the permutation factor for bringing together a n_i = c^_i c_i .
                if(Anfound.first){
                    Anfactor = Anfactor*pow(-1.0 , CrString.size() + AnString.size() - Anfound.second - 1);
                    diagonal.push_back(number);
                    CrString.erase(CrString.begin() + Anfound.second);
                }
                else{
                    AnString.push_back(number);
                }
            }
            else{
                Crfactor = pow(-1.0 , AnString.size());
                CrString.push_back(number);
            }
            perm[number] += creation;
            coeff = coeff*Crfactor*Anfactor;
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
                pair<bool , int> Pconjfound = Perm_found(Perm_conjugate(perm) , Data.Permutations);
                if(Pconjfound.first){
                    pair<bool, int> Dfound = Diag_found(diagonal , Data.Diagonals[Pconjfound.second]);
                    // If the diagonal of a conjugate operator is not found, we add it to the list
                    // Otherwise, we ignore and only ensure that the coefficients are added with their hermitian counterparts to ensure hermiticity!
                    pair<int, int> cran_lengths = Perm_lengths(perm);
                    double sgn_factor = pow(-1.0 ,cran_lengths.first + cran_lengths.second);
                    complex<double> coeffconj = conj(coeff)*sgn_factor;
                    if(!Dfound.first){
                        Data.Coeffs[Pconjfound.second].push_back(coeffconj);
                        Data.Diagonals[Pconjfound.second].push_back(diagonal);
                    }
                    else if(coeffconj != Data.Coeffs[Pconjfound.second][Dfound.second]){
                        // Adding the conjugate to the list if there's a permutation conjugate that is not a full conjugate of an 
                        //      existing term!
                        Data.Coeffs[Pconjfound.second][Dfound.second] += coeffconj;
                    }
                }
                else{
                    Data.Coeffs.push_back({coeff});
                    Data.Diagonals.push_back({diagonal});
                    Data.Permutations.push_back(perm);
                }
            }
        }
    }
    inputFile.close();
    return Data;
}

int main(int argc , char* argv[]) {
    vector<vector<double>> PermMatrix;
    NumberOps D0;
    vector<NumberOps> Diags;
    vector<vector<complex<double>>> Cs;
    DPdataElec ElecData;

    string filename(argv[1]);
    ElecData = Data_to_Perms(filename);

    PermMatrix = Transpose(ElecData.Permutations);
    Diags = ElecData.Diagonals;
    Cs = ElecData.Coeffs;
    D0 = ElecData.D0;
    vector<complex<double>> C0 = ElecData.Coeffs0;
    vector<vector<int>> Cycles = Nullspace(PermMatrix);


    cout << "Running the conjugate coeffs" << endl;

    vector<vector<complex<double>>> Cs_conj = Conjugate_coeffs(Cs , Transpose(PermMatrix));

    cout << endl;
    cout << "The permutations are: "  << endl;
    Print_matrix(PermMatrix);
    cout << endl;

    cout << "The diagonals are: "  << endl;
    for(int i=0; i < Diags.size(); ++i){
        cout << "Diag_" << i << " is " << endl;
        Print_matrix(Diags[i]);
    }
    cout << endl;

    cout << "The coefficients are: " << endl;
    Print_matrix(Cs);
    cout << endl;

    cout << "The conjugate coefficients are: " << endl;
    Print_matrix(Cs_conj);
    cout << endl;

    cout << "The purely diagonal terms are: " << endl;
    Print_matrix(D0);
    cout << endl;

    cout << "The coefficient for the purely diagonal are: " << endl;
    for(int i = 0; i < C0.size(); ++i){
        cout << C0[i] << endl;
    }
    cout << endl;

    // Minimizing cycles lengths and printing them!
    Cycles = Transpose(Cycles);
    while(Cycle_minimize(Cycles));
    cout << "The cycles after minimization are: " << endl;
    Cycles = Transpose(Cycles);
    Print_Matrix(Cycles);

    return 0;
}