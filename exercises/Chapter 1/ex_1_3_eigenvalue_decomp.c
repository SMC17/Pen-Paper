/*
 * Exercise 1.3: Eigenvalue Decomposition
 * 
 * Proves properties of eigenvalues and eigenvectors.
 * Demonstrates A = UΛV^T decomposition and matrix inversion via eigendecomposition.
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#define TOL 1e-10
#define N 3  // Matrix size for tests

// Matrix-vector multiply: result = A * v
void matvec(double *result, double A[][N], const double *v, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
}

// Matrix-matrix multiply: C = A * B
void matmul(double C[][N], double A[][N], double B[][N], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Check if Au = λu (eigenvector equation)
bool is_eigenpair(double A[][N], const double *u, double lambda, int n) {
    double Au[N], lambda_u[N];
    matvec(Au, A, u, n);
    for (int i = 0; i < n; i++) {
        lambda_u[i] = lambda * u[i];
    }
    
    for (int i = 0; i < n; i++) {
        if (fabs(Au[i] - lambda_u[i]) > TOL) return false;
    }
    return true;
}

// Vector addition: result = alpha*u1 + beta*u2
void vec_linear_combo(double *result, double alpha, const double *u1, 
                      double beta, const double *u2, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = alpha * u1[i] + beta * u2[i];
    }
}

// Create diagonal matrix from eigenvalues
void make_diagonal(double Lambda[][N], const double *eigenvals, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Lambda[i][j] = (i == j) ? eigenvals[i] : 0.0;
        }
    }
}

// Matrix from column vectors: A = [v1 v2 v3]
void matrix_from_columns(double A[][N], double v[][N], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = v[j][i];  // v[j] is j-th column
        }
    }
}

// Simple 3x3 matrix inverse (for testing)
bool invert_3x3(double inv[][N], double A[][N]) {
    double det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
               - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
               + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
    
    if (fabs(det) < TOL) return false;
    
    inv[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1]) / det;
    inv[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2]) / det;
    inv[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1]) / det;
    inv[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2]) / det;
    inv[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0]) / det;
    inv[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2]) / det;
    inv[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0]) / det;
    inv[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1]) / det;
    inv[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0]) / det;
    
    return true;
}

// Compare matrices
bool matrices_equal(double A[][N], double B[][N], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(A[i][j] - B[i][j]) > TOL) return false;
        }
    }
    return true;
}

void print_matrix(const char *name, double A[][N], int n) {
    printf("%s =\n", name);
    for (int i = 0; i < n; i++) {
        printf("  [");
        for (int j = 0; j < n; j++) {
            printf("%7.3f", A[i][j]);
        }
        printf(" ]\n");
    }
}

// Test (a): Linear combination of eigenvectors
bool test_1_3_a(void) {
    printf("\n=== Test 1.3(a): Linear combo of eigenvectors ===\n");
    
    // Simple matrix with known eigenvectors
    double A[N][N] = {
        {2, 0, 0},
        {0, 2, 0},
        {0, 0, 3}
    };
    
    // Two eigenvectors with eigenvalue 2
    double u1[N] = {1, 0, 0};
    double u2[N] = {0, 1, 0};
    double lambda = 2.0;
    
    // Linear combination: u = 3*u1 + 2*u2
    double u[N];
    vec_linear_combo(u, 3.0, u1, 2.0, u2, N);
    
    printf("Matrix: diagonal with λ=2,2,3\n");
    printf("u1=[%.1f,%.1f,%.1f] with λ=%.1f: %s\n", 
           u1[0], u1[1], u1[2], lambda, 
           is_eigenpair(A, u1, lambda, N) ? "✓" : "✗");
    printf("u2=[%.1f,%.1f,%.1f] with λ=%.1f: %s\n", 
           u2[0], u2[1], u2[2], lambda,
           is_eigenpair(A, u2, lambda, N) ? "✓" : "✗");
    printf("u=3*u1+2*u2=[%.1f,%.1f,%.1f] with λ=%.1f: %s\n", 
           u[0], u[1], u[2], lambda,
           is_eigenpair(A, u, lambda, N) ? "✓" : "✗");
    
    return is_eigenpair(A, u, lambda, N);
}

// Test (b): Matrix form AU = UΛ
bool test_1_3_b(void) {
    printf("\n=== Test 1.3(b): Matrix form AU = UΛ ===\n");
    
    // Diagonal matrix (eigenvalues on diagonal)
    double A[N][N] = {
        {4, 0, 0},
        {0, 2, 0},
        {0, 0, 1}
    };
    
    // Eigenvectors (identity for diagonal matrix)
    double U_cols[N][N] = {
        {1, 0, 0},  // u1
        {0, 1, 0},  // u2
        {0, 0, 1}   // u3
    };
    double U[N][N];
    matrix_from_columns(U, U_cols, N);
    
    double eigenvals[N] = {4, 2, 1};
    double Lambda[N][N];
    make_diagonal(Lambda, eigenvals, N);
    
    // Compute AU and UΛ
    double AU[N][N], U_Lambda[N][N];
    matmul(AU, A, U, N);
    matmul(U_Lambda, U, Lambda, N);
    
    print_matrix("AU", AU, N);
    print_matrix("UΛ", U_Lambda, N);
    
    bool equal = matrices_equal(AU, U_Lambda, N);
    printf("AU = UΛ: %s\n", equal ? "✓" : "✗");
    
    return equal;
}

// Test (c): Eigenvalue decomposition A = UΛV^T and A^(-1)
bool test_1_3_c(void) {
    printf("\n=== Test 1.3(c): A = UΛV^T decomposition ===\n");
    
    // Use a diagonal matrix (makes V=U for simplicity)
    double A_orig[N][N] = {
        {3, 0, 0},
        {0, 2, 0},
        {0, 0, 1}
    };
    
    double U_cols[N][N] = {{1,0,0}, {0,1,0}, {0,0,1}};
    double U[N][N], V[N][N];  // For diagonal, V=U
    matrix_from_columns(U, U_cols, N);
    matrix_from_columns(V, U_cols, N);
    
    double eigenvals[N] = {3, 2, 1};
    double Lambda[N][N], Lambda_inv[N][N];
    make_diagonal(Lambda, eigenvals, N);
    
    // Λ^(-1)
    double inv_eigenvals[N] = {1.0/3, 1.0/2, 1.0/1};
    make_diagonal(Lambda_inv, inv_eigenvals, N);
    
    // Reconstruct A = UΛV^T
    double temp[N][N], A_recon[N][N];
    matmul(temp, U, Lambda, N);
    matmul(A_recon, temp, V, N);  // V is already transposed for diagonal
    
    printf("Original A:\n");
    print_matrix("A", A_orig, N);
    printf("\nReconstructed A = UΛV^T:\n");
    print_matrix("A_recon", A_recon, N);
    
    bool decomp_ok = matrices_equal(A_orig, A_recon, N);
    printf("\nA = UΛV^T: %s\n", decomp_ok ? "✓" : "✗");
    
    // Test A^(-1) = UΛ^(-1)V^T
    double A_inv_direct[N][N], A_inv_eigen[N][N];
    invert_3x3(A_inv_direct, A_orig);
    
    matmul(temp, U, Lambda_inv, N);
    matmul(A_inv_eigen, temp, V, N);
    
    printf("\nDirect inverse:\n");
    print_matrix("A^(-1)", A_inv_direct, N);
    printf("\nEigen inverse UΛ^(-1)V^T:\n");
    print_matrix("A^(-1)_eigen", A_inv_eigen, N);
    
    bool inv_ok = matrices_equal(A_inv_direct, A_inv_eigen, N);
    printf("\nA^(-1) = UΛ^(-1)V^T: %s\n", inv_ok ? "✓" : "✗");
    
    return decomp_ok && inv_ok;
}

// Test (c) sum form: A = Σᵢ λᵢuᵢvᵢᵀ (PDF equations 1.8, 1.9)
bool test_1_3_c_sum_form(void) {
    printf("\n=== Test 1.3(c): Sum form A = Σᵢ λᵢuᵢvᵢᵀ ===\n");
    
    // Use diagonal matrix for simplicity (V=U)
    double A_orig[N][N] = {
        {3, 0, 0},
        {0, 2, 0},
        {0, 0, 1}
    };
    
    // Eigenvectors and eigenvalues
    double u_cols[N][N] = {{1,0,0}, {0,1,0}, {0,0,1}};
    double eigenvals[N] = {3, 2, 1};
    
    // Reconstruct A = Σᵢ λᵢuᵢvᵢᵀ
    // For each i: add λᵢ * uᵢ * vᵢᵀ to result
    double A_sum[N][N] = {{0}};
    for (int i = 0; i < N; i++) {
        // λᵢ * uᵢ * vᵢᵀ (outer product)
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                A_sum[row][col] += eigenvals[i] * u_cols[i][row] * u_cols[i][col];
            }
        }
    }
    
    printf("Original A:\n");
    print_matrix("A", A_orig, N);
    printf("\nReconstructed A = Σᵢ λᵢuᵢvᵢᵀ:\n");
    print_matrix("A_sum", A_sum, N);
    
    bool sum_ok = matrices_equal(A_orig, A_sum, N);
    printf("\nA = Σᵢ λᵢuᵢvᵢᵀ: %s\n", sum_ok ? "✓" : "✗");
    
    // Now test A^(-1) = Σᵢ (1/λᵢ)uᵢvᵢᵀ
    double A_inv_sum[N][N] = {{0}};
    for (int i = 0; i < N; i++) {
        // (1/λᵢ) * uᵢ * vᵢᵀ (outer product)
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                A_inv_sum[row][col] += (1.0/eigenvals[i]) * u_cols[i][row] * u_cols[i][col];
            }
        }
    }
    
    double A_inv_direct[N][N];
    invert_3x3(A_inv_direct, A_orig);
    
    printf("\nDirect inverse:\n");
    print_matrix("A^(-1)", A_inv_direct, N);
    printf("\nSum form A^(-1) = Σᵢ (1/λᵢ)uᵢvᵢᵀ:\n");
    print_matrix("A^(-1)_sum", A_inv_sum, N);
    
    bool inv_sum_ok = matrices_equal(A_inv_direct, A_inv_sum, N);
    printf("\nA^(-1) = Σᵢ (1/λᵢ)uᵢvᵢᵀ: %s\n", inv_sum_ok ? "✓" : "✗");
    
    return sum_ok && inv_sum_ok;
}

// Edge case: zero eigenvalue (singular matrix)
bool test_edge_singular(void) {
    printf("\n=== Edge case: Singular matrix (zero eigenvalue) ===\n");
    
    double A[N][N] = {
        {0, 0, 0},
        {0, 2, 0},
        {0, 0, 3}
    };
    
    double u[N] = {1, 0, 0};
    double lambda = 0.0;
    
    printf("Eigenvector u=[%.1f,%.1f,%.1f] with λ=%.1f\n", 
           u[0], u[1], u[2], lambda);
    
    bool is_eigen = is_eigenpair(A, u, lambda, N);
    printf("Au = λu: %s\n", is_eigen ? "✓" : "✗");
    printf("Note: Matrix is singular (not invertible)\n");
    
    return is_eigen;
}

int main(void) {
    int passed = 0, total = 0;
    
    total++; if (test_1_3_a()) passed++;
    total++; if (test_1_3_b()) passed++;
    total++; if (test_1_3_c()) passed++;
    total++; if (test_1_3_c_sum_form()) passed++;
    total++; if (test_edge_singular()) passed++;
    
    printf("\n%s: %d/%d tests passed\n", 
           passed == total ? "✓ SUCCESS" : "✗ FAILED", passed, total);
    return passed == total ? 0 : 1;
}
