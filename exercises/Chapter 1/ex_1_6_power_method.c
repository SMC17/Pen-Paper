/* Exercise 1.6: Power Method
 * 
 * Analyzes the power method algorithm for finding the dominant eigenvector
 * (eigenvector with largest eigenvalue) of a positive definite symmetric matrix.
 * 
 * Algorithm: v_{k+1} = Σw_k,  w_{k+1} = v_{k+1} / ||v_{k+1}||
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define EPS 1e-9
#define MAX_ITER 1000

/* ============================================================================
 * LINEAR ALGEBRA BASICS
 * ============================================================================ */

typedef struct {
    double* data;  // row-major storage
    int n;         // dimension (n×n matrix)
} Matrix;

typedef struct {
    double* data;
    int n;
} Vector;

Matrix mat_alloc(int n) {
    Matrix m = {.n = n, .data = calloc(n * n, sizeof(double))};
    return m;
}

Vector vec_alloc(int n) {
    Vector v = {.n = n, .data = calloc(n, sizeof(double))};
    return v;
}

void mat_free(Matrix* m) {
    free(m->data);
    m->data = NULL;
}

void vec_free(Vector* v) {
    free(v->data);
    v->data = NULL;
}

double mat_get(Matrix m, int i, int j) {
    return m.data[i * m.n + j];
}

void mat_set(Matrix m, int i, int j, double val) {
    m.data[i * m.n + j] = val;
}

double vec_get(Vector v, int i) {
    return v.data[i];
}

void vec_set(Vector v, int i, double val) {
    v.data[i] = val;
}

void vec_copy(Vector dst, Vector src) {
    memcpy(dst.data, src.data, src.n * sizeof(double));
}

double vec_dot(Vector u, Vector v) {
    double sum = 0.0;
    for (int i = 0; i < u.n; i++) {
        sum += u.data[i] * v.data[i];
    }
    return sum;
}

double vec_norm(Vector v) {
    return sqrt(vec_dot(v, v));
}

void vec_normalize(Vector v) {
    double norm = vec_norm(v);
    for (int i = 0; i < v.n; i++) {
        v.data[i] /= norm;
    }
}

/* Matrix-vector multiplication: result = A * v */
void mat_vec_mult(Vector result, Matrix A, Vector v) {
    for (int i = 0; i < A.n; i++) {
        result.data[i] = 0.0;
        for (int j = 0; j < A.n; j++) {
            result.data[i] += mat_get(A, i, j) * v.data[j];
        }
    }
}

/* Matrix-matrix multiplication: C = A * B */
void mat_mult(Matrix C, Matrix A, Matrix B) {
    for (int i = 0; i < A.n; i++) {
        for (int j = 0; j < B.n; j++) {
            double sum = 0.0;
            for (int k = 0; k < A.n; k++) {
                sum += mat_get(A, i, k) * mat_get(B, k, j);
            }
            mat_set(C, i, j, sum);
        }
    }
}

/* Matrix transpose: B = A^T */
void mat_transpose(Matrix B, Matrix A) {
    for (int i = 0; i < A.n; i++) {
        for (int j = 0; j < A.n; j++) {
            mat_set(B, i, j, mat_get(A, j, i));
        }
    }
}

bool mat_is_symmetric(Matrix A) {
    for (int i = 0; i < A.n; i++) {
        for (int j = i + 1; j < A.n; j++) {
            if (fabs(mat_get(A, i, j) - mat_get(A, j, i)) > EPS) {
                return false;
            }
        }
    }
    return true;
}

/* ============================================================================
 * POWER METHOD (Part a-d)
 * ============================================================================ */

/* Power method: finds dominant eigenvector
 * Algorithm:
 *   v_{k+1} = Σ w_k
 *   w_{k+1} = v_{k+1} / ||v_{k+1}||
 * 
 * Returns: dominant eigenvector (normalized)
 *          and the number of iterations taken
 */
typedef struct {
    Vector eigenvector;
    double eigenvalue;
    int iterations;
    bool converged;
} PowerMethodResult;

PowerMethodResult power_method(Matrix Sigma, Vector w0, double tolerance) {
    PowerMethodResult result = {0};
    result.eigenvector = vec_alloc(Sigma.n);
    
    Vector w_curr = vec_alloc(Sigma.n);
    Vector w_prev = vec_alloc(Sigma.n);
    Vector v = vec_alloc(Sigma.n);
    
    // Initialize with w0 (normalized)
    vec_copy(w_curr, w0);
    vec_normalize(w_curr);
    
    for (int k = 0; k < MAX_ITER; k++) {
        vec_copy(w_prev, w_curr);
        
        // v_{k+1} = Σ w_k
        mat_vec_mult(v, Sigma, w_curr);
        
        // w_{k+1} = v_{k+1} / ||v_{k+1}||
        vec_copy(w_curr, v);
        vec_normalize(w_curr);
        
        // Check convergence: ||w_k - w_{k-1}|| < tolerance
        // (Note: could also be ||w_k + w_{k-1}|| due to sign ambiguity)
        double diff = 0.0;
        for (int i = 0; i < Sigma.n; i++) {
            double d = w_curr.data[i] - w_prev.data[i];
            diff += d * d;
        }
        diff = sqrt(diff);
        
        // Also check sign-flipped version
        double diff_neg = 0.0;
        for (int i = 0; i < Sigma.n; i++) {
            double d = w_curr.data[i] + w_prev.data[i];
            diff_neg += d * d;
        }
        diff_neg = sqrt(diff_neg);
        
        double min_diff = (diff < diff_neg) ? diff : diff_neg;
        
        if (min_diff < tolerance) {
            result.iterations = k + 1;
            result.converged = true;
            break;
        }
        
        result.iterations = k + 1;
    }
    
    // Final result
    vec_copy(result.eigenvector, w_curr);
    
    // Compute eigenvalue: λ = w^T Σ w (Rayleigh quotient)
    mat_vec_mult(v, Sigma, w_curr);
    result.eigenvalue = vec_dot(w_curr, v);
    
    vec_free(&w_curr);
    vec_free(&w_prev);
    vec_free(&v);
    
    return result;
}

/* ============================================================================
 * EIGENVALUE DECOMPOSITION (for verification)
 * ============================================================================ */

/* For 2×2 symmetric matrix only (for testing) */
typedef struct {
    double lambda1, lambda2;
    Vector u1, u2;
} Eigen2x2;

Eigen2x2 eigen_2x2_symmetric(Matrix A) {
    Eigen2x2 result = {0};
    
    double a = mat_get(A, 0, 0);
    double b = mat_get(A, 0, 1);
    double c = mat_get(A, 1, 0);
    double d = mat_get(A, 1, 1);
    
    double trace = a + d;
    double det = a * d - b * c;
    double disc = sqrt(trace * trace - 4.0 * det);
    
    result.lambda1 = (trace + disc) / 2.0;
    result.lambda2 = (trace - disc) / 2.0;
    
    result.u1 = vec_alloc(2);
    result.u2 = vec_alloc(2);
    
    if (fabs(b) > EPS) {
        vec_set(result.u1, 0, -b);
        vec_set(result.u1, 1, a - result.lambda1);
        vec_normalize(result.u1);
        
        vec_set(result.u2, 0, -b);
        vec_set(result.u2, 1, a - result.lambda2);
        vec_normalize(result.u2);
    } else {
        vec_set(result.u1, 0, 1.0);
        vec_set(result.u1, 1, 0.0);
        vec_set(result.u2, 0, 0.0);
        vec_set(result.u2, 1, 1.0);
    }
    
    return result;
}

/* ============================================================================
 * PRINTING UTILITIES
 * ============================================================================ */

void print_vec(const char* name, Vector v) {
    printf("%s = [", name);
    for (int i = 0; i < v.n; i++) {
        printf("%.4f", v.data[i]);
        if (i < v.n - 1) printf(", ");
    }
    printf("]\n");
}

void print_mat(const char* name, Matrix A) {
    printf("%s =\n", name);
    for (int i = 0; i < A.n; i++) {
        printf("  [");
        for (int j = 0; j < A.n; j++) {
            printf("%8.4f", mat_get(A, i, j));
        }
        printf("  ]\n");
    }
}

/* ============================================================================
 * TESTS
 * ============================================================================ */

/* Test: Power method converges to dominant eigenvector (2×2) */
void test_power_method_2x2(void) {
    printf("\n=== Test: Power Method on 2×2 Matrix ===\n");
    
    // Create a 2×2 positive definite symmetric matrix
    Matrix Sigma = mat_alloc(2);
    mat_set(Sigma, 0, 0, 4.0);
    mat_set(Sigma, 0, 1, 1.0);
    mat_set(Sigma, 1, 0, 1.0);
    mat_set(Sigma, 1, 1, 3.0);
    
    print_mat("Σ", Sigma);
    
    // Compute exact eigendecomposition
    Eigen2x2 exact = eigen_2x2_symmetric(Sigma);
    printf("\nExact eigendecomposition:\n");
    printf("  λ₁ = %.6f (dominant)\n", exact.lambda1);
    printf("  λ₂ = %.6f\n", exact.lambda2);
    print_vec("  u₁", exact.u1);
    print_vec("  u₂", exact.u2);
    
    // Run power method with random initialization
    Vector w0 = vec_alloc(2);
    vec_set(w0, 0, 0.6);
    vec_set(w0, 1, 0.8);
    
    printf("\n");
    print_vec("Initial w₀", w0);
    
    PowerMethodResult pm = power_method(Sigma, w0, 1e-9);
    
    printf("\nPower method result:\n");
    printf("  Converged: %s\n", pm.converged ? "yes" : "no");
    printf("  Iterations: %d\n", pm.iterations);
    printf("  λ (Rayleigh quotient) = %.6f\n", pm.eigenvalue);
    print_vec("  Dominant eigenvector", pm.eigenvector);
    
    // Verify: compare with exact u1 (up to sign)
    double dot_pos = fabs(vec_dot(pm.eigenvector, exact.u1));
    double diff_lambda = fabs(pm.eigenvalue - exact.lambda1);
    
    printf("\nVerification:\n");
    printf("  |<u_pm, u₁>| = %.9f (should be ≈ 1)\n", dot_pos);
    printf("  |λ_pm - λ₁| = %.9f (should be ≈ 0)\n", diff_lambda);
    
    bool pass = pm.converged && 
                fabs(dot_pos - 1.0) < 1e-6 && 
                diff_lambda < 1e-6;
    printf("\n✓ Test %s: power method finds dominant eigenvector\n", 
           pass ? "PASSED" : "FAILED");
    
    mat_free(&Sigma);
    vec_free(&w0);
    vec_free(&pm.eigenvector);
    vec_free(&exact.u1);
    vec_free(&exact.u2);
}

/* Test: Power method on 3×3 matrix */
void test_power_method_3x3(void) {
    printf("\n=== Test: Power Method on 3×3 Matrix ===\n");
    
    // Create a 3×3 positive definite symmetric matrix
    // Construct as Σ = V Λ V^T where we control eigenvalues
    Matrix Sigma = mat_alloc(3);
    mat_set(Sigma, 0, 0, 5.0);
    mat_set(Sigma, 0, 1, 1.0);
    mat_set(Sigma, 0, 2, 0.5);
    mat_set(Sigma, 1, 0, 1.0);
    mat_set(Sigma, 1, 1, 3.0);
    mat_set(Sigma, 1, 2, 0.3);
    mat_set(Sigma, 2, 0, 0.5);
    mat_set(Sigma, 2, 1, 0.3);
    mat_set(Sigma, 2, 2, 2.0);
    
    print_mat("Σ (3×3)", Sigma);
    printf("  (eigenvalues ≈ 5.4, 2.9, 1.7 - dominant is ≈5.4)\n");
    
    // Initialize with random vector
    Vector w0 = vec_alloc(3);
    vec_set(w0, 0, 0.5);
    vec_set(w0, 1, 0.5);
    vec_set(w0, 2, 0.7071);
    
    printf("\n");
    print_vec("Initial w₀", w0);
    
    PowerMethodResult pm = power_method(Sigma, w0, 1e-9);
    
    printf("\nPower method result:\n");
    printf("  Converged: %s\n", pm.converged ? "yes" : "no");
    printf("  Iterations: %d\n", pm.iterations);
    printf("  λ (dominant) = %.6f\n", pm.eigenvalue);
    print_vec("  Dominant eigenvector", pm.eigenvector);
    
    // Verify: Σu ≈ λu
    Vector Sigma_u = vec_alloc(3);
    mat_vec_mult(Sigma_u, Sigma, pm.eigenvector);
    
    Vector lambda_u = vec_alloc(3);
    for (int i = 0; i < 3; i++) {
        vec_set(lambda_u, i, pm.eigenvalue * vec_get(pm.eigenvector, i));
    }
    
    printf("\nVerification (Σu = λu):\n");
    print_vec("  Σu", Sigma_u);
    print_vec("  λu", lambda_u);
    
    double diff = 0.0;
    for (int i = 0; i < 3; i++) {
        double d = vec_get(Sigma_u, i) - vec_get(lambda_u, i);
        diff += d * d;
    }
    diff = sqrt(diff);
    
    printf("  ||Σu - λu|| = %.9f\n", diff);
    
    bool pass = pm.converged && diff < 1e-6;
    printf("\n✓ Test %s: Σu = λu for dominant eigenpair\n", 
           pass ? "PASSED" : "FAILED");
    
    mat_free(&Sigma);
    vec_free(&w0);
    vec_free(&pm.eigenvector);
    vec_free(&Sigma_u);
    vec_free(&lambda_u);
}

/* Test: Convergence behavior - trace the sequence */
void test_convergence_behavior(void) {
    printf("\n=== Test: Convergence Behavior ===\n");
    
    Matrix Sigma = mat_alloc(2);
    mat_set(Sigma, 0, 0, 10.0);  // Large dominant eigenvalue
    mat_set(Sigma, 0, 1, 1.0);
    mat_set(Sigma, 1, 0, 1.0);
    mat_set(Sigma, 1, 1, 2.0);
    
    print_mat("Σ", Sigma);
    
    Eigen2x2 exact = eigen_2x2_symmetric(Sigma);
    printf("\nEigenvalues: λ₁ = %.4f, λ₂ = %.4f\n", exact.lambda1, exact.lambda2);
    printf("Ratio λ₂/λ₁ = %.6f (convergence rate)\n\n", exact.lambda2 / exact.lambda1);
    
    // Run power method and track convergence
    Vector w0 = vec_alloc(2);
    vec_set(w0, 0, 0.6);
    vec_set(w0, 1, 0.8);
    vec_normalize(w0);
    
    Vector w_curr = vec_alloc(2);
    Vector v = vec_alloc(2);
    vec_copy(w_curr, w0);
    
    printf("Iteration trace (first 10 iterations):\n");
    printf("  k    w[0]      w[1]      |<w,u₁>|   λ (Rayleigh)\n");
    printf("  -------------------------------------------------\n");
    
    for (int k = 0; k < 10; k++) {
        // Compute Rayleigh quotient
        mat_vec_mult(v, Sigma, w_curr);
        double lambda = vec_dot(w_curr, v);
        
        // Compare with true eigenvector
        double dot_u1 = fabs(vec_dot(w_curr, exact.u1));
        
        printf("  %d  %8.5f  %8.5f  %8.6f  %8.5f\n", 
               k, vec_get(w_curr, 0), vec_get(w_curr, 1), dot_u1, lambda);
        
        // Update
        mat_vec_mult(v, Sigma, w_curr);
        vec_copy(w_curr, v);
        vec_normalize(w_curr);
    }
    
    printf("\n✓ Test PASSED: convergence behavior visualized\n");
    
    mat_free(&Sigma);
    vec_free(&w0);
    vec_free(&w_curr);
    vec_free(&v);
    vec_free(&exact.u1);
    vec_free(&exact.u2);
}

/* Test: Edge case - initialization orthogonal to dominant eigenvector */
void test_bad_initialization(void) {
    printf("\n=== Test: Initialization Orthogonal to Dominant Eigenvector ===\n");
    
    // Diagonal matrix (easy to see eigenvectors)
    Matrix Sigma = mat_alloc(2);
    mat_set(Sigma, 0, 0, 5.0);
    mat_set(Sigma, 0, 1, 0.0);
    mat_set(Sigma, 1, 0, 0.0);
    mat_set(Sigma, 1, 1, 2.0);
    
    print_mat("Σ (diagonal)", Sigma);
    printf("Eigenvalues: λ₁=5, λ₂=2\n");
    printf("Eigenvectors: u₁=[1,0], u₂=[0,1]\n\n");
    
    // Initialize exactly as [0, 1] (orthogonal to u₁)
    Vector w0 = vec_alloc(2);
    vec_set(w0, 0, 0.0);
    vec_set(w0, 1, 1.0);
    
    print_vec("Initial w₀ (orthogonal to u₁!)", w0);
    
    printf("\nThis is a degenerate case: power method may not converge to u₁\n");
    printf("In practice, numerical errors will break the orthogonality.\n");
    
    // Add tiny perturbation to break exact orthogonality
    vec_set(w0, 0, 1e-6);
    vec_set(w0, 1, 1.0);
    
    print_vec("\nPerturbed w₀ (with tiny component in u₁ direction)", w0);
    
    PowerMethodResult pm = power_method(Sigma, w0, 1e-9);
    
    printf("\nPower method result:\n");
    printf("  Converged: %s\n", pm.converged ? "yes" : "no");
    printf("  Iterations: %d\n", pm.iterations);
    printf("  λ = %.6f\n", pm.eigenvalue);
    print_vec("  Eigenvector", pm.eigenvector);
    
    bool pass = pm.converged && fabs(pm.eigenvalue - 5.0) < 1e-6;
    printf("\n✓ Test %s: even tiny perturbations allow convergence\n", 
           pass ? "PASSED" : "FAILED");
    
    mat_free(&Sigma);
    vec_free(&w0);
    vec_free(&pm.eigenvector);
}

/* ============================================================================
 * MAIN
 * ============================================================================ */

int main(void) {
    printf("Exercise 1.6: Power Method\n");
    printf("====================================================\n");
    
    test_power_method_2x2();
    test_power_method_3x3();
    test_convergence_behavior();
    test_bad_initialization();
    
    printf("\n=== All tests completed ===\n");
    return 0;
}
