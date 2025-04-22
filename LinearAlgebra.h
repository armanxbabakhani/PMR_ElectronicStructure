#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
typedef vector<vector<int>> Matrix;

template<typename T>
void Print_Matrix(const vector<vector<T>>& matrix) {
    int m = matrix.size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            cout << matrix[i][j] << "  ";
        }
        cout << endl;
        cout << endl;
    }
}

void Gauss_Elim(vector<vector<double>>& A) {
    int rowsA = A.size() , colsA = A[0].size() , lead = 0;

    for (int r = 0; r < rowsA; r++) {
        if (lead >= colsA) {
            break;
        }
        int i = r;
        while (A[i][lead] == 0) {
            i++;
            if (i == rowsA) {
                i = r;
                lead++;
                if (lead == colsA) {
                    return;
                }
            }
        }
        
        swap(A[i], A[r]);

        int lv = A[r][lead]; // Multiplier is the number that multiplies the leading term to 1.
        for (int j = 0; j < colsA; j++) {
            A[r][j] = A[r][j]/lv;
        }
        for (int i = 0; i < rowsA; i++) {
            if (i != r) {
                int lv2 = A[i][lead];
                for (int j = 0; j < colsA; j++) {
                    A[i][j] -= (A[r][j] * lv2);
                }
            }
        }
        lead++;
    }
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

bool All_zeros(Matrix A){
    for(int i = 0; i < A.size(); i++){
        for(int j = 0; j < A[0].size(); j++){
            if(A[i][j] != 0){
                return false;
            }
        }
    }
    return true;
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

bool CheckOnesInt(vector<int> vec){
    for(int i = 0; i < vec.size(); i++){
        if(abs(abs(vec[i])-1) < 1E-6 || abs(vec[i]) < 1E-6){
        }
        else{
            return false;
        }
    }
    return true;
}

vector<int> Round_toint(vector<double> vec){
    vector<int> output;
    for(int i = 0; i < vec.size() ; ++i){
        int num = static_cast<int>(round(vec[i]));
        output.push_back(num);
    }
    return output;
}

int modInverse(int a, int mod) {
    // Assumes mod is a prime number
    int m = mod - 2;  // Fermat's little theorem
    int result = 1;
    while (m) {
        if (m & 1)
            result = (result * a) % mod;
        a = (a * a) % mod;
        m >>= 1;
    }
    return result;
}

Matrix Modp_nullspace(const Matrix& A, int p) {
    int n = A.size();        // number of rows
    int m = A[0].size();     // number of columns
    Matrix B = A;            // Make a copy of A to transform into RREF
    Matrix nullspace;
    
    for(int i = 0; i < n ; ++i){
        for(int j = 0; j < m; ++j){
            B[i][j] = B[i][j]%p;
        }
    }
    
    if(All_zeros(B)){
        return {{0}};
    }


    // Step 1: Transform B into RREF form mod p
    int row = 0;
    vector<int> lead_col(m, -1);  // Track leading columns
    for (int col = 0; col < m && row < n; col++) {
        int sel = row;
        // Find a non-zero entry in the current column
        for (int i = row; i < n; i++) {
            if (B[i][col] % p != 0) {
                sel = i;
                break;
            }
        }

        if (B[sel][col] % p == 0)
            continue;  // Skip this column if no non-zero entries

        // Swap rows if necessary
        if (sel != row) {
            swap(B[sel], B[row]);
        }

        // Scale row to make the leading coefficient 1
        int inv = modInverse(B[row][col], p);
        for (int j = 0; j < m; j++) {
            B[row][j] = (B[row][j] * inv) % p;
            if (B[row][j] < 0) B[row][j] += p;
        }

        // Zero out all other entries in this column
        for (int i = 0; i < n; i++) {
            if (i != row) {
                int c = B[i][col];
                for (int j = 0; j < m; j++) {
                    B[i][j] = (B[i][j] - c * B[row][j] % p + p) % p;
                }
            }
        }

        lead_col[row] = col;  // Mark the leading column
        row++;
    }

    // Step 2: Identify free variables and form nullspace basis vectors
    nullspace.clear();
    vector<int> free_var_index;
    for (int j = 0; j < m; j++) {
        if (find(lead_col.begin(), lead_col.end(), j) == lead_col.end()) {
            free_var_index.push_back(j);  // j is a free variable
        }
    }

    for (int idx : free_var_index) {
        vector<int> basis_vec(m, 0);
        basis_vec[idx] = 1;  // Set the free variable
        for (int i = 0; i < row; i++) {
            if (lead_col[i] != -1) {
                basis_vec[lead_col[i]] = -B[i][idx];
                if (basis_vec[lead_col[i]] < 0)
                    basis_vec[lead_col[i]] += p;
            }
        }
        nullspace.push_back(basis_vec);
    }
    return Transpose(nullspace);
}

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


vector<vector<int>> Nullspace(vector<vector<double>> A) {
    vector<vector<int>> nullspaceInt;
    vector<vector<double>> nullspace;
    int colsA = A[0].size(), rowsA = A.size();
    vector<int> PivotCols, PivotRows;

    // Step 1: Gaussian Elimination
    Gauss_Elim(A);

    // Step 2: Identify pivot columns
    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsA; j++) {
            if (abs(A[i][j] - 1.0) < 1E-6) {
                PivotCols.push_back(j);
                PivotRows.push_back(i);
                break;
            }
        }
    }

    // Step 3: Identify free variables and build rational nullspace
    for (int j = 0; j < colsA; j++) {
        if (!Num_found_inVec(j, PivotCols)) {
            vector<double> basis(colsA, 0.0);
            basis[j] = 1.0;
            for (int k = 0; k < PivotCols.size(); k++) {
                basis[PivotCols[k]] = -A[PivotRows[k]][j];
            }
            nullspace.push_back(basis);
        }
    }

    // Step 4: Convert each basis vector to scaled integers
    for (const auto& vec : nullspace) {
        int lcm_denom = 1;
        for (double val : vec) {
            double frac_part = val - round(val);
            if (abs(frac_part) > 1e-8) {
                int denom = 1;
                while (abs(round(val * denom) - val * denom) > 1e-6 && denom < 10000)
                    denom++;
                lcm_denom = lcm(lcm_denom, denom);
            }
        }

        vector<int> int_vec;
        for (double val : vec) {
            int_vec.push_back(static_cast<int>(round(val * lcm_denom)));
        }
        nullspaceInt.push_back(int_vec);
    }

    // Step 5: Remove duplicates
    set<vector<int>> unique_basis(nullspaceInt.begin(), nullspaceInt.end());
    return Transpose(vector<vector<int>>(unique_basis.begin(), unique_basis.end()));
}


vector<int> Vec_add(vector<int> v1 , vector<int> v2){
    if(v1.size() != v2.size()){
        cerr << "Vec_add Error: The sizes of the two input vectors are different!" << endl;
        return {0};
    }
    vector<int> sumv12;
    for(int i = 0; i < v1.size() ; ++i){
        sumv12.push_back(v1[i] + v2[i]);
    }
    return sumv12;
}

vector<int> Vec_sub(vector<int> v1 , vector<int> v2){
    if(v1.size() != v2.size()){
        cerr << "Vec_add Error: The sizes of the two input vectors are different!" << endl;
        return {0};
    }
    vector<int> sumv12;
    for(int i = 0; i < v1.size() ; ++i){
        sumv12.push_back(v1[i] - v2[i]);
    }
    return sumv12;
}

int Abs_sum(vector<int> vec){
    int sum = 0;
    for(int i = 0; i < vec.size(); i++){
        sum += abs(vec[i]);
    }
    return sum;
}

bool compare_null(const vector<int> & a, const vector<int> & b){
	return Abs_sum(a) < Abs_sum(b);
}

int Cycle_minimize(vector<vector<int>>& null_eigs){ 
	int nullsize, k, m, null_k, changes_made = 0; vector<int> curr;
	nullsize = null_eigs.size();
	sort(null_eigs.begin(), null_eigs.end(), compare_null);
	for(k = nullsize-1; k > 0 ; k--){
		null_k = Abs_sum(null_eigs[k]);
		for(m = 0 ; m < k; m++){
			curr = Vec_sub( null_eigs[k] , null_eigs[m]);
			if(Abs_sum(curr) < null_k){
				null_eigs[k] = curr; changes_made = 1;
				break;
			}
		}
	}
	return changes_made;
}
