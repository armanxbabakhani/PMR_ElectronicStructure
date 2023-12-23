#include <iostream>
#include <vector>
#include "LinearAlgebra.h"

using namespace std;

int main() {
    // Example usage
    vector<vector<double>> matrix = {
        {1, 1, 0},
        {-1, 0, -1},
        {0, -1, 1}
    };

    vector<vector<double>> matrixRed = matrix;
    Gauss_Elim(matrixRed);
    vector<vector<int>> nullspaceBasis = Nullspace(matrix);

    // Print the result
    cout << "The reduced matrix is:" << endl;
    for (const auto& row : matrixRed) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    cout << endl;

    // Print the result
    cout << "Nullspace of the matrix:" << endl;
    for (const auto& row : nullspaceBasis) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    return 0;
}