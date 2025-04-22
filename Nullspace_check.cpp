#include <iostream>
#include <vector>
#include <set>
#include <numeric>
#include <cmath>
#include "LinearAlgebra.h"

using namespace std;

// Function to compute LCM of two numbers
int lcm(int a, int b) {
    return abs(a * b) / gcd(a, b);
}

// Function to compute LCM of denominators in a vector of doubles
int lcm_of_denominators(const vector<double>& vec) {
    int denom_lcm = 1;
    for (double val : vec) {
        double frac_part = val - static_cast<int>(val);
        if (fabs(frac_part) > 1e-8) {
            int denom = 1;
            while (fabs(round(val * denom) - val * denom) > 1e-6 && denom < 10000) {
                denom++;
            }
            denom_lcm = lcm(denom_lcm, denom);
        }
    }
    return denom_lcm;
}

// Convert vector<double> to scaled integer vector
vector<int> convert_to_integer_vector(const vector<double>& vec) {
    int scale = lcm_of_denominators(vec);
    vector<int> result;
    for (double val : vec) {
        result.push_back(static_cast<int>(round(val * scale)));
    }
    return result;
}

// Remove duplicates from a list of integer vectors
vector<vector<int>> remove_duplicates(const vector<vector<int>>& vecs) {
    set<vector<int>> unique_set(vecs.begin(), vecs.end());
    return vector<vector<int>>(unique_set.begin(), unique_set.end());
}

int main() {
    // Example matrix with values -1, 0, 1
    vector<vector<double>> A = {
        {1 , -1 , 0, 0, 0},
        {-1 , 1 , -1, 1, 0},
        {-1 , 1 , 1, -1, 1},
        {0 , 0 , 1, -1, -1}
    };

    A = Transpose(A);
    // Compute nullspace using the original Nullspace function
    vector<vector<int>> nullspace = Nullspace(A);

    cout << "Initial Nullspace (±1 only filter):" << endl;
    Print_Matrix(nullspace);

    // Optional: Recompute using full Gaussian elimination and convert to integers
    Gauss_Elim(A);
    vector<vector<double>> raw_nullspace;
    int cols = A[0].size(), rows = A.size();
    vector<int> pivots;

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            if (fabs(A[i][j] - 1.0) < 1e-6) {
                pivots.push_back(j);
                break;
            }

    for (int i = 0; i < cols; ++i) {
        if (!Num_found_inVec(i, pivots)) {
            vector<double> basis(cols, 0.0);
            basis[i] = 1.0;
            for (int k = 0; k < pivots.size(); ++k) {
                basis[pivots[k]] = -A[k][i];
            }
            raw_nullspace.push_back(basis);
        }
    }

    // Convert and deduplicate
    vector<vector<int>> integer_nullspace;
    for (auto vec : raw_nullspace)
        integer_nullspace.push_back(convert_to_integer_vector(vec));

    integer_nullspace = remove_duplicates(integer_nullspace);

    cout << "\nInteger-scaled, deduplicated nullspace vectors:" << endl;
    Print_Matrix(integer_nullspace);

    return 0;
}
