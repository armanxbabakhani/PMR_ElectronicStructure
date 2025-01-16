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
typedef vector<vector<int>> Matrix;

struct DPdataElec{
    vector<vector<complex<double>>> Coeffs , CoeffsConj;
    vector<NumberOps> Diagonals;
    CrAnOps Permutations;
    vector<complex<double>> Coeffs0;
    vector<vector<int>> D0;  
};

template<typename T>
void print_matrix(const vector<vector<T>>& matrix , string name) {
    int m = matrix.size();
    for (int i = 0; i < m; i++) {
        cout << name << "[" << i << "] = ";
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
        vector<double> Perms_i(Perms[i].begin() , Perms[i].end());
        pair<int, int> cran_lengths = Perm_lengths( Perms_i );
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
        bool StateKilled = false;
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
                pair<bool , int> Crfound = Find_number(number, AnString);
                // If two excitations of the same orbital occurs, state is killed (fermionic statistics)!
                if(!Crfound.first){
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
                    StateKilled = true;
                }
            }
            else{
                pair<bool , int> Anfound = Find_number(number, CrString);
                if(!Anfound.first){
                    Crfactor = pow(-1.0 , AnString.size());
                    CrString.push_back(number);
                }
                // If two annihilations of the same orbital occurs, state is killed (fermionic statistics)!
                else{
                    StateKilled = true;
                }
            }

            // Checking if any of the annihilations or excitations have occurred repeatedly to kill the entire state:
            if(!StateKilled){
                perm[number] += creation;
                coeff = coeff*Crfactor*Anfactor;
            }
            else{
                coeff = 0.0;
                perm = vector<double> (NumOfParticles+1 , 0.0);
                break;
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
                pair<bool , int> Pconjfound = Perm_found(Perm_conjugate(perm) , Data.Permutations);
                if(Pconjfound.first){
                    pair<bool, int> Dfound = Diag_found(diagonal , Data.Diagonals[Pconjfound.second]);
                    // If the diagonal of a conjugate operator is not found, we add it to the list
                    // Otherwise, we ignore and only ensure that the coefficients are added with their hermitian counterparts to ensure hermiticity!
                    pair<int, int> cran_lengths = Perm_lengths(perm);
                    double sgn_factor = pow(-1.0 , cran_lengths.first + cran_lengths.second);
                    complex<double> coeffconj = conj(coeff)*sgn_factor;
                    //cout << "The conjugate coefficient is " << coeffconj << endl;
                    //cout << endl;
                    if(!Dfound.first){
                        Data.Coeffs[Pconjfound.second].push_back(coeffconj);
                        Data.Diagonals[Pconjfound.second].push_back(diagonal);
                    }
                    /*
                    else if(coeffconj != Data.Coeffs[Pconjfound.second][Dfound.second]){
                        // Adding the conjugate to the list if there's a permutation conjugate that is not a full conjugate of an 
                        //      existing term!
                        //Data.Coeffs[Pconjfound.second][Dfound.second] += coeffconj;
                    }*/
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

int determinant(const Matrix& matrix) {
    int n = matrix.size();

    // Base case for 1x1 matrix
    if (n == 1) {
        return matrix[0][0];
    }

    // Base case for 2x2 matrix
    if (n == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }

    int det = 0;

    // Recursive case: expand along the first row
    for (int col = 0; col < n; ++col) {
        // Create a submatrix by excluding the first row and the current column
        vector<vector<int>> submatrix(n - 1, vector<int>(n - 1));
        for (int i = 1; i < n; ++i) {
            int subCol = 0;
            for (int j = 0; j < n; ++j) {
                if (j == col) continue; // Skip the current column
                submatrix[i - 1][subCol++] = matrix[i][j];
            }
        }

        // Compute the cofactor and add to determinant
        int sign = (col % 2 == 0) ? 1 : -1;
        det += sign * matrix[0][col] * determinant(submatrix);
    }

    return det;
}

void remove_wrong_cycles(vector<vector<int>> &NullspaceInt, const vector<vector<int>> &PermutationMatrices) {
    // Ensure dimensions match between NullspaceInt and PermutationMatrices
    cout << "The dimension of NullspaceInt: " << NullspaceInt.size() << " and " << NullspaceInt[0].size() << endl;

    cout << "The dimension of Permutation: " << PermutationMatrices.size() << " and " << PermutationMatrices[0].size() << endl;
    if (NullspaceInt.empty() || PermutationMatrices.empty() || NullspaceInt[0].size() != PermutationMatrices.size()) {
        cerr << "Dimension mismatch between NullspaceInt and PermutationMatrices." << endl;
        return;
    }

    // Resultant vector for checking each condition
    vector<vector<int>> validVectors;

    for (const auto &nullVec : NullspaceInt) {
        // Initialize a vector<int> to store the sum result for this nullVec
        vector<int> sumResult(PermutationMatrices[0].size(), 0);

        for (size_t i = 0; i < nullVec.size(); ++i) {
            // Check if the coefficient in nullVec is 0, -1, or 1
            if (nullVec[i] != 0 && nullVec[i] != 1 && nullVec[i] != -1) {
                cerr << "Invalid coefficient found in NullspaceInt." << endl;
                continue;
            }

            // Add the scaled PermutationMatrices[i] to sumResult
            for (size_t j = 0; j < PermutationMatrices[i].size(); ++j) {
                sumResult[j] += nullVec[i] * PermutationMatrices[i][j];
            }
        }

        // Check if sumResult is a vector with all entries zero
        bool isValid = all_of(sumResult.begin(), sumResult.end(), [](int val) { return val == 0; });

        if (isValid) {
            validVectors.push_back(nullVec); // Keep this vector as it satisfies the condition
        }
    }

    // Update NullspaceInt to only include valid vectors
    NullspaceInt = validVectors;
}


Matrix find_all_electronic_cycles(Matrix PermutationMatrices){
    // The permutation matrices contain 1 , 0 , or -1 , so we will convert this into mod 3 linear algebra
    //  by replacing the -1 with 2

    // Replacing -1 with 2
    Matrix PermMatrixMod3 = PermutationMatrices;
    for(int i=0; i < PermutationMatrices.size(); ++i){
        for(int j=0; j < PermutationMatrices[i].size(); ++j){
            if(PermutationMatrices[i][j] == -1){
                PermMatrixMod3[i][j] = 2;
            }
        }
    }
    cout << "Right before Modp nullspace computation! " << endl;

    Matrix NullspaceInt , Nullspace = Modp_nullspace(PermMatrixMod3 , 3);

    cout << "Mod nullspace is computed! " << endl;

    // Converting the 2 back to -1:
    for(int i=0; i < Nullspace.size(); ++i){
        for(int j=0; j < Nullspace[i].size(); ++j){
            if(Nullspace[i][j] == 2){
                Nullspace[i][j] = -1;
            }
        }
    }

    for(int i = 0; i < Nullspace.size() ; i++){
        if(CheckOnesInt(Nullspace[i])){
            vector<int> nullint = Nullspace[i];
            //vector<int> nullint(nullspace[i].begin() , nullspace[i].end());
            NullspaceInt.push_back(nullint);
        }
    }

    // Remove the wrong cycles:
    NullspaceInt = Transpose(NullspaceInt);
    remove_wrong_cycles(NullspaceInt , Transpose(PermutationMatrices));

    return NullspaceInt;
}

int main(int argc , char* argv[]) {

    NumberOps D0;
    vector<NumberOps> Diags;
    vector<vector<complex<double>>> Cs;
    DPdataElec ElecData;

    string filename(argv[1]);
    ElecData = Data_to_Perms(filename);
    vector<vector<double>> PermMatrixDouble = Transpose(ElecData.Permutations);
    vector<vector<int>> PermMatrix(PermMatrixDouble.size() , vector<int>(PermMatrixDouble[0].size() , 0));
    
    for (size_t i = 0; i < PermMatrixDouble.size(); ++i) {
        for(size_t j =0; j < PermMatrixDouble[i].size(); j++){
            PermMatrix[i][j] = static_cast<int>(PermMatrixDouble[i][j]); // Converts double to int (truncates the decimal part)
        }
    }

    Diags = ElecData.Diagonals;
    Cs = ElecData.Coeffs;
    D0 = ElecData.D0;
    vector<complex<double>> C0 = ElecData.Coeffs0;
    cout << "Starting to compute the nullspace! " << endl;
    // vector<vector<int>> Cycles = Nullspace(PermMatrixDouble);
    vector<vector<int>> Cycles = find_all_electronic_cycles(PermMatrix);
    cout << "Nullspace computation is done! " << endl;
    Print_Matrix(Cycles);

    //vector<vector<complex<double>>> Cs_conj = Conjugate_coeffs(Cs , Transpose(PermMatrixDouble));

    cout << endl;
    cout << "The permutations are: "  << endl;
    print_matrix(Transpose(PermMatrix) , "Permutations");
    cout << endl;

    // cout << "The diagonals are: "  << endl;
    // for(int i=0; i < Diags.size(); ++i){
    //     cout << "Diag_" << i << " is " << endl;
    //     Print_matrix(Diags[i]);
    // }
    // cout << endl;

    // cout << "The coefficients are: " << endl;
    // Print_matrix(Cs);
    // cout << endl;

    // cout << "The conjugate coefficients are: " << endl;
    // Print_matrix(Cs_conj);
    // cout << endl;

    // cout << "The purely diagonal terms are: " << endl;
    // Print_matrix(D0);
    // cout << endl;

    // cout << "The coefficient for the purely diagonal are: " << endl;
    // for(int i = 0; i < C0.size(); ++i){
    //     cout << C0[i] << endl;
    // }
    // cout << endl;

    // Minimizing cycles lengths and printing them!
    for(int i=0; i < Cycles.size(); ++i){
        for(int j=0; j < Cycles[i].size(); ++j){
            if(Cycles[i][j] == 2){
                Cycles[i][j] = -1;
            }
        }
    }

    Cycles = Transpose(Cycles);
    //while(Cycle_minimize(Cycles));
    //cout << "The cycles after minimization are: " << endl;
    //Cycles = Transpose(Cycles);
    print_matrix(Transpose(Cycles) , "cycles");

    return 0;
}