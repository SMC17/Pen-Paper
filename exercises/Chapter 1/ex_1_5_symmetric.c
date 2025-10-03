/* Exercise 1.5: Eigenvalue decomposition for symmetric matrices
 * 
 * Part (a): Show that eigenvectors of a symmetric matrix with distinct
 *           eigenvalues are orthogonal.
 * Part (b): Show that positive definiteness <=> all eigenvalues > 0,
 *           and that positive definite matrices are invertible.
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define EPS 1e-9

/* ============================================================================
 * BASIC LINEAR ALGEBRA
 * ============================================================================ */

typedef struct {
    double x, y;
} Vec2;

typedef struct {
    double a, b, c, d;  // [a b; c d]
} Mat2;

double dot(Vec2 u, Vec2 v) {
    return u.x * v.x + u.y * v.y;
}

double norm(Vec2 v) {
    return sqrt(dot(v, v));
}

Vec2 mat_vec_mult(Mat2 A, Vec2 v) {
    return (Vec2){
        .x = A.a * v.x + A.b * v.y,
        .y = A.c * v.x + A.d * v.y
    };
}

double det(Mat2 A) {
    return A.a * A.d - A.b * A.c;
}

double trace(Mat2 A) {
    return A.a + A.d;
}

bool is_symmetric(Mat2 A) {
    return fabs(A.b - A.c) < EPS;
}

/* ============================================================================
 * EIGENVALUE DECOMPOSITION FOR 2x2 SYMMETRIC MATRICES
 * ============================================================================ */

typedef struct {
    double lambda1, lambda2;  // eigenvalues
    Vec2 u1, u2;              // eigenvectors (normalized)
    bool distinct;             // whether eigenvalues are distinct
} EigenDecomp;

/* For 2x2 matrix, eigenvalues are roots of characteristic polynomial:
 * det(A - λI) = λ² - trace(A)λ + det(A) = 0
 * λ = (trace ± sqrt(trace² - 4*det)) / 2
 */
EigenDecomp eigen_decomp_2x2_symmetric(Mat2 A) {
    EigenDecomp result = {0};
    
    if (!is_symmetric(A)) {
        fprintf(stderr, "Warning: matrix is not symmetric\n");
    }
    
    double tr = trace(A);
    double dt = det(A);
    double discriminant = tr * tr - 4.0 * dt;
    
    if (discriminant < 0) {
        fprintf(stderr, "Error: complex eigenvalues (not expected for symmetric matrix)\n");
        return result;
    }
    
    double sqrt_disc = sqrt(discriminant);
    result.lambda1 = (tr + sqrt_disc) / 2.0;
    result.lambda2 = (tr - sqrt_disc) / 2.0;
    result.distinct = fabs(result.lambda1 - result.lambda2) > EPS;
    
    /* Find eigenvector for lambda1: (A - lambda1*I)u1 = 0
     * [a-λ₁   b  ] [x]   [0]
     * [c    d-λ₁] [y] = [0]
     * 
     * If b ≠ 0: from first row, x = -b/(a-λ₁) * y, choose y=1
     * If b = 0: eigenvector is [1, 0] or [0, 1]
     */
    if (fabs(A.b) > EPS) {
        double x1 = -A.b;
        double y1 = A.a - result.lambda1;
        double n1 = sqrt(x1*x1 + y1*y1);
        result.u1 = (Vec2){x1/n1, y1/n1};
        
        double x2 = -A.b;
        double y2 = A.a - result.lambda2;
        double n2 = sqrt(x2*x2 + y2*y2);
        result.u2 = (Vec2){x2/n2, y2/n2};
    } else {
        // Diagonal matrix
        result.u1 = (Vec2){1.0, 0.0};
        result.u2 = (Vec2){0.0, 1.0};
    }
    
    return result;
}

/* ============================================================================
 * POSITIVE DEFINITENESS
 * ============================================================================ */

/* Check if v^T A v > 0 for a specific vector v */
double quadratic_form(Mat2 A, Vec2 v) {
    Vec2 Av = mat_vec_mult(A, v);
    return dot(v, Av);
}

/* A matrix is positive definite if v^T A v > 0 for all v ≠ 0
 * For 2x2 symmetric matrices, we can check:
 * - All eigenvalues > 0, OR
 * - Sample test: check on several non-zero vectors
 */
bool is_positive_definite_by_eigenvalues(EigenDecomp eig) {
    return eig.lambda1 > EPS && eig.lambda2 > EPS;
}

bool is_positive_definite_by_sampling(Mat2 A) {
    /* Test on a set of sample vectors
     * This is not exhaustive but good for testing
     */
    Vec2 test_vectors[] = {
        {1.0, 0.0},
        {0.0, 1.0},
        {1.0, 1.0},
        {1.0, -1.0},
        {0.5, 0.866},  // 30 degrees
        {0.707, 0.707}, // 45 degrees
    };
    
    for (int i = 0; i < 6; i++) {
        if (quadratic_form(A, test_vectors[i]) <= EPS) {
            return false;
        }
    }
    return true;
}

/* ============================================================================
 * MATRIX INVERSION
 * ============================================================================ */

/* Inverse of 2x2 matrix A = [a b; c d]
 * A^(-1) = (1/det(A)) * [d -b; -c a]
 */
Mat2 inverse(Mat2 A) {
    double d = det(A);
    if (fabs(d) < EPS) {
        fprintf(stderr, "Error: matrix is singular (not invertible)\n");
        return (Mat2){0};
    }
    return (Mat2){
        .a = A.d / d,
        .b = -A.b / d,
        .c = -A.c / d,
        .d = A.a / d
    };
}

Mat2 mat_mult(Mat2 A, Mat2 B) {
    return (Mat2){
        .a = A.a * B.a + A.b * B.c,
        .b = A.a * B.b + A.b * B.d,
        .c = A.c * B.a + A.d * B.c,
        .d = A.c * B.b + A.d * B.d
    };
}

bool is_identity(Mat2 A) {
    return fabs(A.a - 1.0) < EPS && 
           fabs(A.b) < EPS && 
           fabs(A.c) < EPS && 
           fabs(A.d - 1.0) < EPS;
}

/* ============================================================================
 * TESTS
 * ============================================================================ */

void print_vec(const char* name, Vec2 v) {
    printf("%s = [%.4f, %.4f]\n", name, v.x, v.y);
}

void print_mat(const char* name, Mat2 A) {
    printf("%s = [%.4f  %.4f]\n", name, A.a, A.b);
    printf("     [%.4f  %.4f]\n", A.c, A.d);
}

void print_eigen(EigenDecomp eig) {
    printf("λ₁ = %.6f, λ₂ = %.6f\n", eig.lambda1, eig.lambda2);
    print_vec("u₁", eig.u1);
    print_vec("u₂", eig.u2);
    printf("Distinct eigenvalues: %s\n", eig.distinct ? "yes" : "no");
}

/* Test Part (a): Eigenvectors with distinct eigenvalues are orthogonal */
void test_orthogonality_of_eigenvectors(void) {
    printf("\n=== Test: Orthogonality of Eigenvectors (Part a) ===\n");
    
    Mat2 A = {.a = 3.0, .b = 1.0, .c = 1.0, .d = 3.0};  // symmetric
    print_mat("A", A);
    printf("\n");
    
    EigenDecomp eig = eigen_decomp_2x2_symmetric(A);
    print_eigen(eig);
    printf("\n");
    
    // Verify Au₁ = λ₁u₁
    Vec2 Au1 = mat_vec_mult(A, eig.u1);
    Vec2 lambda1_u1 = {eig.lambda1 * eig.u1.x, eig.lambda1 * eig.u1.y};
    double diff1 = sqrt(pow(Au1.x - lambda1_u1.x, 2) + pow(Au1.y - lambda1_u1.y, 2));
    printf("Verify Au₁ = λ₁u₁:\n");
    print_vec("  Au₁", Au1);
    print_vec("  λ₁u₁", lambda1_u1);
    printf("  Difference: %.9f\n", diff1);
    
    // Verify Au₂ = λ₂u₂
    Vec2 Au2 = mat_vec_mult(A, eig.u2);
    Vec2 lambda2_u2 = {eig.lambda2 * eig.u2.x, eig.lambda2 * eig.u2.y};
    double diff2 = sqrt(pow(Au2.x - lambda2_u2.x, 2) + pow(Au2.y - lambda2_u2.y, 2));
    printf("\nVerify Au₂ = λ₂u₂:\n");
    print_vec("  Au₂", Au2);
    print_vec("  λ₂u₂", lambda2_u2);
    printf("  Difference: %.9f\n", diff2);
    
    // Check orthogonality
    double dot_product = dot(eig.u1, eig.u2);
    printf("\nOrthogonality check:\n");
    printf("  u₁ · u₂ = %.9f\n", dot_product);
    printf("  |u₁ · u₂| = %.9f\n", fabs(dot_product));
    
    bool pass = eig.distinct && 
                diff1 < EPS && 
                diff2 < EPS && 
                fabs(dot_product) < EPS;
    printf("\n✓ Test %s: eigenvectors are orthogonal\n", pass ? "PASSED" : "FAILED");
}

/* Test Part (b): Positive definiteness <=> all eigenvalues > 0 */
void test_positive_definiteness(void) {
    printf("\n=== Test: Positive Definiteness (Part b) ===\n");
    
    // Positive definite matrix
    Mat2 A_pos = {.a = 4.0, .b = 1.0, .c = 1.0, .d = 3.0};
    print_mat("A (positive definite)", A_pos);
    
    EigenDecomp eig_pos = eigen_decomp_2x2_symmetric(A_pos);
    printf("  λ₁ = %.6f, λ₂ = %.6f\n", eig_pos.lambda1, eig_pos.lambda2);
    
    bool pd_by_eig = is_positive_definite_by_eigenvalues(eig_pos);
    bool pd_by_sampling = is_positive_definite_by_sampling(A_pos);
    
    printf("  Positive definite by eigenvalues: %s\n", pd_by_eig ? "yes" : "no");
    printf("  Positive definite by sampling: %s\n", pd_by_sampling ? "yes" : "no");
    printf("  Match: %s\n", (pd_by_eig == pd_by_sampling) ? "yes" : "no");
    
    // NOT positive definite (one negative eigenvalue)
    printf("\n");
    Mat2 A_neg = {.a = 1.0, .b = 2.0, .c = 2.0, .d = 1.0};
    print_mat("B (not positive definite)", A_neg);
    
    EigenDecomp eig_neg = eigen_decomp_2x2_symmetric(A_neg);
    printf("  λ₁ = %.6f, λ₂ = %.6f\n", eig_neg.lambda1, eig_neg.lambda2);
    
    bool nd_by_eig = is_positive_definite_by_eigenvalues(eig_neg);
    bool nd_by_sampling = is_positive_definite_by_sampling(A_neg);
    
    printf("  Positive definite by eigenvalues: %s\n", nd_by_eig ? "yes" : "no");
    printf("  Positive definite by sampling: %s\n", nd_by_sampling ? "yes" : "no");
    printf("  Match: %s\n", (nd_by_eig == nd_by_sampling) ? "yes" : "no");
    
    // Test specific vector that makes it negative
    Vec2 v_test = {1.0, -1.0};  // Should give negative quadratic form for A_neg
    double qf = quadratic_form(A_neg, v_test);
    printf("  v^T B v for v=[1,-1]: %.6f (should be negative)\n", qf);
    
    bool pass = (pd_by_eig == pd_by_sampling) && 
                (nd_by_eig == nd_by_sampling) &&
                pd_by_eig && !nd_by_eig &&
                qf < 0;
    printf("\n✓ Test %s: positive definiteness <=> all eigenvalues > 0\n", 
           pass ? "PASSED" : "FAILED");
}

/* Test Part (b): Positive definite => invertible */
void test_invertibility(void) {
    printf("\n=== Test: Positive Definite => Invertible (Part b) ===\n");
    
    Mat2 A = {.a = 4.0, .b = 1.0, .c = 1.0, .d = 3.0};
    print_mat("A (positive definite)", A);
    
    EigenDecomp eig = eigen_decomp_2x2_symmetric(A);
    printf("  λ₁ = %.6f, λ₂ = %.6f (both > 0)\n", eig.lambda1, eig.lambda2);
    
    Mat2 A_inv = inverse(A);
    print_mat("\nA^(-1)", A_inv);
    
    // Verify A * A^(-1) = I
    Mat2 product = mat_mult(A, A_inv);
    printf("\n");
    print_mat("A * A^(-1)", product);
    
    bool is_I = is_identity(product);
    printf("Is identity: %s\n", is_I ? "yes" : "no");
    
    bool pass = is_positive_definite_by_eigenvalues(eig) && is_I;
    printf("\n✓ Test %s: positive definite matrix is invertible\n", 
           pass ? "PASSED" : "FAILED");
}

/* Edge case: degenerate (not positive definite) */
void test_non_invertible(void) {
    printf("\n=== Test: Non-positive-definite => May Not Be Invertible ===\n");
    
    // Singular matrix (determinant = 0)
    Mat2 A = {.a = 1.0, .b = 2.0, .c = 2.0, .d = 4.0};
    print_mat("A (singular)", A);
    printf("  det(A) = %.6f\n", det(A));
    
    EigenDecomp eig = eigen_decomp_2x2_symmetric(A);
    printf("  λ₁ = %.6f, λ₂ = %.6f\n", eig.lambda1, eig.lambda2);
    printf("  Has zero eigenvalue: %s\n", 
           (fabs(eig.lambda2) < EPS) ? "yes" : "no");
    
    bool pd = is_positive_definite_by_eigenvalues(eig);
    printf("  Positive definite: %s\n", pd ? "yes" : "no");
    printf("  Invertible: %s\n", (fabs(det(A)) > EPS) ? "yes" : "no");
    
    bool pass = !pd && fabs(det(A)) < EPS;
    printf("\n✓ Test %s: zero eigenvalue => not invertible\n", 
           pass ? "PASSED" : "FAILED");
}

/* ============================================================================
 * MAIN
 * ============================================================================ */

int main(void) {
    printf("Exercise 1.5: Eigenvalue Decomposition for Symmetric Matrices\n");
    printf("================================================================\n");
    
    test_orthogonality_of_eigenvectors();
    test_positive_definiteness();
    test_invertibility();
    test_non_invertible();
    
    printf("\n=== All tests completed ===\n");
    return 0;
}
