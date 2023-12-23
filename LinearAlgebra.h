#include <iostream>
#include <vector>

using namespace std;
// Function to perform Gaussian elimination with back-substitution
void Gauss_Elim(vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int i = 0; i < cols; ++i) {
        // Find pivot row
        int pivotRow = -1, lead;
        for (int j = i; j < rows; ++j) {
            if (matrix[j][i] != 0) {
                pivotRow = j;
                lead = i;
                break;
            }
        }

        // If no pivot found, continue to the next column
        if (pivotRow == -1) {
            continue;
        }

        /*
        for(int k = 0; k < cols; k++){
            matrix[i][k] = matrix[i][k]/matrix[i][i];
        }*/

        // Swap rows to make the pivot element non-zero
        swap(matrix[i], matrix[pivotRow]);

        // Eliminate other rows
        for (int j = 0; j < rows; ++j) {
            if (j != i && matrix[j][i] != 0) {
                double factor = matrix[j][i] / matrix[i][i];
                for (int k = i; k < cols; ++k) {
                    matrix[j][k] -= factor * matrix[i][k];
                }
            }
        }
        double factor = matrix[i][i];
        if(factor != 0){
            for(int k = i; k < cols; k++){
                matrix[i][k] = matrix[i][k]/factor;
            }
        }
    }
}

// Function to find the nullspace basis of a matrix
vector<vector<double>> Nullspace_old(vector<vector<double>> matrix) {
    // Create an augmented matrix [A | I]
    int rows = matrix.size();
    int cols = matrix[0].size();
    vector<vector<double>> AugmentedMatrix = matrix;

    // Apply Gaussian elimination with back-substitution to the augmented matrix
    Gauss_Elim(AugmentedMatrix);

    // Extract the nullspace basis from the right side of the augmented matrix
    vector<vector<double>> NullspaceBasis(rows, vector<double>(rows, 0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < rows; ++j) {
            NullspaceBasis[j][i] = -AugmentedMatrix[j][cols + i];
        }
    }

    return NullspaceBasis;
}

template<typename T>
vector<vector<T>> Transpose(vector<vector<T>> A){
    vector<vector<T>> A_transpose;
    int rowsA = A.size();
    if(rowsA == 0){
        return A_transpose;
    }
    else{
        int colsA = A[0].size();
        for(int j=0; j < colsA ; j++){
            vector<T> A_transpose_j(rowsA , 0);
            for(int i = 0; i < rowsA ; i++){
                A_transpose_j[i] = A[i][j];
            }
            A_transpose.push_back(A_transpose_j);
        }
        return A_transpose;
    }
}

bool Num_found_inVec(int num , vector<int> Vec){
    bool found = false;
    for (int i=0; i < Vec.size(); i++){
        if(num == Vec[i]){
            found = true;
            break;
        }
    }
    return found;
}

bool CheckOnes(vector<double> vec){
    for(int i = 0; i < vec.size(); i++){
        if(abs(abs(vec[i])-1) < 1E-6 || abs(vec[i]) < 1E-6){
        }
        else{
            return false;
        }
    }
    return true;
}

vector<vector<int>> Nullspace(vector<vector<double>> A) {
    vector<vector<double>> nullspace;
    vector<vector<int>> nullspaceInt;
    int colsA = A[0].size() , rowsA = A.size() , rank=0 , nullspaceDim;
    vector<int> PivotCols , PivotRows;

    Gauss_Elim(A);
    //Determine pivot columns:
    for(int i = 0; i < rowsA; i++){
        for(int j = 0; j < colsA; j++){
            if(abs(A[i][j]-1) < 1E-6){
                PivotCols.push_back(j);
                PivotRows.push_back(i);
                rank++;
                break;
            }
        }
    }

    nullspaceDim = colsA - rank;

    for(int i=0; i < colsA; i++){
        if(!Num_found_inVec(i, PivotCols)){ 
            vector<double> BasisVec(A[0].size() , 0);
            BasisVec[i] = 1;
            // Mark the rows that have a non-zero entry in this column:
            for(int k = 0; k < rank; k++){
                BasisVec[PivotCols[k]] = -1.0*A[PivotRows[k]][i];
            }
            nullspace.push_back(BasisVec);
        }
    }
    //nullspace = Transpose(nullspace);

    for(int i = 0; i < nullspace.size() ; i++){
        if(CheckOnes(nullspace[i])){
            cout << "Inside the checkones!! " << endl;
            vector<int> nullint(nullspace[i].begin() , nullspace[i].end());
            nullspaceInt.push_back(nullint);
        }
    }

    return Transpose(nullspaceInt);
}


