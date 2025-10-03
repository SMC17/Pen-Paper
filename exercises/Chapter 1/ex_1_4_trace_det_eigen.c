/*
 * Exercise 1.4: Trace, Determinants and Eigenvalues
 * 
 * Proves fundamental relationships:
 * - tr(A) = Σλᵢ (trace equals sum of eigenvalues)
 * - det(A) = ∏λᵢ (determinant equals product of eigenvalues)
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define TOL 1e-10
#define N 3

// Trace: sum of diagonal elements
double trace(double A[][N], int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += A[i][i];
    }
    return sum;
}

// Determinant (3x3 only for simplicity)
double det3(double A[][N]) {
    return A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1])
         - A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0])
         + A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);
}

// Sum of eigenvalues
double sum_eigenvalues(const double *lambda, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += lambda[i];
    }
    return sum;
}

// Product of eigenvalues
double product_eigenvalues(const double *lambda, int n) {
    double prod = 1.0;
    for (int i = 0; i < n; i++) {
        prod *= lambda[i];
    }
    return prod;
}

void print_matrix(const char *name, double A[][N], int n) {
    printf("%s =\n", name);
    for (int i = 0; i < n; i++) {
        printf("  [");
        for (int j = 0; j < n; j++) {
            printf("%6.2f", A[i][j]);
        }
        printf(" ]\n");
    }
}

// Test (a): tr(A) = Σλᵢ
bool test_1_4_a(void) {
    printf("\n=== Test 1.4(a): tr(A) = Σλᵢ ===\n");
    
    // Diagonal matrix (eigenvalues on diagonal)
    double A[N][N] = {
        {4, 0, 0},
        {0, 2, 0},
        {0, 0, 1}
    };
    double eigenvals[N] = {4, 2, 1};
    
    print_matrix("A", A, N);
    
    double tr_A = trace(A, N);
    double sum_lambda = sum_eigenvalues(eigenvals, N);
    
    printf("\ntr(A) = %.6f\n", tr_A);
    printf("Σλᵢ   = %.6f\n", sum_lambda);
    printf("Difference: %.2e\n", fabs(tr_A - sum_lambda));
    
    bool pass = fabs(tr_A - sum_lambda) < TOL;
    printf("tr(A) = Σλᵢ: %s\n", pass ? "✓" : "✗");
    
    return pass;
}

// Test (a) with non-diagonal matrix
bool test_1_4_a_general(void) {
    printf("\n=== Test 1.4(a): Non-diagonal case ===\n");
    
    // Matrix with known eigenvalues (constructed example)
    // For a 2x2 rotation-like matrix scaled
    double A[N][N] = {
        {3, 0, 0},
        {0, 1, 2},
        {0, 2, 4}
    };
    
    // Eigenvalues: 3, 0, 5 (can verify: λ₁=3, block eigenvals 0,5)
    double eigenvals[N] = {3, 0, 5};
    
    print_matrix("A", A, N);
    
    double tr_A = trace(A, N);
    double sum_lambda = sum_eigenvalues(eigenvals, N);
    
    printf("\ntr(A) = %.6f\n", tr_A);
    printf("Σλᵢ   = %.6f\n", sum_lambda);
    printf("Difference: %.2e\n", fabs(tr_A - sum_lambda));
    
    bool pass = fabs(tr_A - sum_lambda) < TOL;
    printf("tr(A) = Σλᵢ: %s\n", pass ? "✓" : "✗");
    
    return pass;
}

// Test (b): det(A) = ∏λᵢ
bool test_1_4_b(void) {
    printf("\n=== Test 1.4(b): det(A) = ∏λᵢ ===\n");
    
    // Diagonal matrix
    double A[N][N] = {
        {3, 0, 0},
        {0, 2, 0},
        {0, 0, 5}
    };
    double eigenvals[N] = {3, 2, 5};
    
    print_matrix("A", A, N);
    
    double det_A = det3(A);
    double prod_lambda = product_eigenvalues(eigenvals, N);
    
    printf("\ndet(A) = %.6f\n", det_A);
    printf("∏λᵢ    = %.6f\n", prod_lambda);
    printf("Difference: %.2e\n", fabs(det_A - prod_lambda));
    
    bool pass = fabs(det_A - prod_lambda) < TOL;
    printf("det(A) = ∏λᵢ: %s\n", pass ? "✓" : "✗");
    
    return pass;
}

// Test both properties together
bool test_combined(void) {
    printf("\n=== Test: Combined properties ===\n");
    
    double A[N][N] = {
        {2, 0, 0},
        {0, 3, 0},
        {0, 0, 4}
    };
    double eigenvals[N] = {2, 3, 4};
    
    print_matrix("A", A, N);
    printf("Eigenvalues: λ₁=%.0f, λ₂=%.0f, λ₃=%.0f\n", 
           eigenvals[0], eigenvals[1], eigenvals[2]);
    
    double tr_A = trace(A, N);
    double sum_lambda = sum_eigenvalues(eigenvals, N);
    double det_A = det3(A);
    double prod_lambda = product_eigenvalues(eigenvals, N);
    
    printf("\nTrace:      tr(A)=%.0f, Σλᵢ=%.0f", tr_A, sum_lambda);
    bool trace_ok = fabs(tr_A - sum_lambda) < TOL;
    printf(" %s\n", trace_ok ? "✓" : "✗");
    
    printf("Determinant: det(A)=%.0f, ∏λᵢ=%.0f", det_A, prod_lambda);
    bool det_ok = fabs(det_A - prod_lambda) < TOL;
    printf(" %s\n", det_ok ? "✓" : "✗");
    
    return trace_ok && det_ok;
}

// Edge case: singular matrix (zero eigenvalue)
bool test_edge_singular(void) {
    printf("\n=== Edge case: Singular matrix ===\n");
    
    double A[N][N] = {
        {0, 0, 0},
        {0, 2, 0},
        {0, 0, 3}
    };
    double eigenvals[N] = {0, 2, 3};
    
    print_matrix("A", A, N);
    printf("Eigenvalues: λ₁=%.0f, λ₂=%.0f, λ₃=%.0f\n", 
           eigenvals[0], eigenvals[1], eigenvals[2]);
    
    double det_A = det3(A);
    double prod_lambda = product_eigenvalues(eigenvals, N);
    
    printf("\ndet(A) = %.6f\n", det_A);
    printf("∏λᵢ    = %.6f\n", prod_lambda);
    printf("Note: Zero eigenvalue → zero determinant → singular\n");
    
    bool pass = fabs(det_A) < TOL && fabs(prod_lambda) < TOL;
    printf("det(A) = ∏λᵢ = 0: %s\n", pass ? "✓" : "✗");
    
    return pass;
}

// Identity matrix test
bool test_identity(void) {
    printf("\n=== Test: Identity matrix ===\n");
    
    double I[N][N] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    double eigenvals[N] = {1, 1, 1};
    
    printf("Identity matrix I (all eigenvalues = 1)\n");
    
    double tr_I = trace(I, N);
    double sum_lambda = sum_eigenvalues(eigenvals, N);
    double det_I = det3(I);
    double prod_lambda = product_eigenvalues(eigenvals, N);
    
    printf("tr(I)=%.0f, Σλᵢ=%.0f ", tr_I, sum_lambda);
    bool trace_ok = fabs(tr_I - sum_lambda) < TOL;
    printf("%s\n", trace_ok ? "✓" : "✗");
    
    printf("det(I)=%.0f, ∏λᵢ=%.0f ", det_I, prod_lambda);
    bool det_ok = fabs(det_I - prod_lambda) < TOL && fabs(det_I - 1.0) < TOL;
    printf("%s\n", det_ok ? "✓" : "✗");
    
    return trace_ok && det_ok;
}

int main(void) {
    int passed = 0, total = 0;
    
    total++; if (test_1_4_a()) passed++;
    total++; if (test_1_4_a_general()) passed++;
    total++; if (test_1_4_b()) passed++;
    total++; if (test_combined()) passed++;
    total++; if (test_edge_singular()) passed++;
    total++; if (test_identity()) passed++;
    
    printf("\n%s: %d/%d tests passed\n", 
           passed == total ? "✓ SUCCESS" : "✗ FAILED", passed, total);
    return passed == total ? 0 : 1;
}
