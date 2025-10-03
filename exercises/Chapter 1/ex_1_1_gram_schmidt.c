/*
 * Exercise 1.1: Gram-Schmidt Orthogonalization
 * 
 * Proves that the Gram-Schmidt process produces orthogonal vectors.
 * Minimal code - only what's needed to demonstrate the mathematics.
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define TOL 1e-10

// Dot product: u^T * v
double dot(const double *u, const double *v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += u[i] * v[i];
    return sum;
}

// Gram-Schmidt: u_i = a_i - sum_{j<i} (u_j^T a_i / u_j^T u_j) u_j
void gram_schmidt(double u[][3], double a[][3], int k) {
    for (int i = 0; i < k; i++) {
        // u[i] = a[i]
        for (int d = 0; d < 3; d++) u[i][d] = a[i][d];
        
        // Subtract projections
        for (int j = 0; j < i; j++) {
            double proj = dot(u[j], a[i], 3) / dot(u[j], u[j], 3);
            for (int d = 0; d < 3; d++) u[i][d] -= proj * u[j][d];
        }
    }
}

void print_vec(const char *name, double v[3]) {
    printf("%s = [%.4f, %.4f, %.4f]\n", name, v[0], v[1], v[2]);
}

bool test_orthogonal(double u[][3], int k, const char *test_name) {
    printf("\n%s:\n", test_name);
    for (int i = 0; i < k; i++) print_vec(i == 0 ? "u1" : i == 1 ? "u2" : "u3", u[i]);
    
    for (int i = 0; i < k; i++) {
        for (int j = i + 1; j < k; j++) {
            double d = dot(u[i], u[j], 3);
            printf("  u%d^T u%d = %.2e\n", i+1, j+1, d);
            if (fabs(d) > TOL) return false;
        }
    }
    return true;
}

int main(void) {
    int passed = 0, total = 0;
    
    // Test (a): Basic case
    {
        total++;
        double a[][3] = {{1,0,0}, {1,1,0}};
        double u[2][3];
        gram_schmidt(u, a, 2);
        if (test_orthogonal(u, 2, "Test 1.1(a): Two vectors")) passed++;
    }
    
    // Test (a) edge: a2 = α*a1 produces zero
    {
        total++;
        double a[][3] = {{1,2,3}, {2,4,6}};
        double u[2][3];
        gram_schmidt(u, a, 2);
        double norm = sqrt(dot(u[1], u[1], 3));
        printf("\nTest 1.1(a) edge: a2=2*a1\n");
        printf("  |u2| = %.2e (should be ~0)\n", norm);
        if (fabs(norm) < TOL) passed++;
    }
    
    // Test (b): Linear combination representation
    // Show that v = αa₁ + βa₂ can be written as γ₁u₁ + γ₂u₂
    {
        total++;
        double a[][3] = {{1,0,0}, {1,1,0}};
        double u[2][3];
        gram_schmidt(u, a, 2);  // u1={1,0,0}, u2={0,1,0}
        
        // Create linear combination v = 2*a1 + 3*a2 = {5, 3, 0}
        double alpha = 2.0, beta = 3.0;
        double v[3] = {
            alpha * a[0][0] + beta * a[1][0],
            alpha * a[0][1] + beta * a[1][1],
            alpha * a[0][2] + beta * a[1][2]
        };
        
        // According to PDF (S.1.12): γ₁ = α + β(u₁ᵀa₂)/(u₁ᵀu₁), γ₂ = β
        double u1_dot_a2 = dot(u[0], a[1], 3);
        double u1_dot_u1 = dot(u[0], u[0], 3);
        double gamma1 = alpha + beta * (u1_dot_a2 / u1_dot_u1);
        double gamma2 = beta;
        
        // Compute v_reconstructed = γ₁u₁ + γ₂u₂
        double v_recon[3] = {
            gamma1 * u[0][0] + gamma2 * u[1][0],
            gamma1 * u[0][1] + gamma2 * u[1][1],
            gamma1 * u[0][2] + gamma2 * u[1][2]
        };
        
        // Check if v = v_reconstructed
        double diff = 0.0;
        for (int i = 0; i < 3; i++) {
            double d = v[i] - v_recon[i];
            diff += d * d;
        }
        diff = sqrt(diff);
        
        printf("\nTest 1.1(b): Linear combination representation\n");
        printf("  v = %.0f*a1 + %.0f*a2 = ", alpha, beta);
        print_vec("v", v);
        printf("  Represented as %.2f*u1 + %.2f*u2\n", gamma1, gamma2);
        print_vec("  v_recon", v_recon);
        printf("  ||v - v_recon|| = %.2e\n", diff);
        
        if (fabs(diff) < TOL) passed++;
    }
    
    // Test (c): Three vectors
    {
        total++;
        double a[][3] = {{1,0,0}, {1,1,0}, {1,1,1}};
        double u[3][3];
        gram_schmidt(u, a, 3);
        if (test_orthogonal(u, 3, "Test 1.1(c): Three vectors")) passed++;
    }
    
    // Test (d): Linear combinations for k vectors
    // Show that v = α₁a₁ + α₂a₂ + α₃a₃ can be written as β₁u₁ + β₂u₂ + β₃u₃
    {
        total++;
        double a[][3] = {{1,0,0}, {1,1,0}, {1,1,1}};
        double u[3][3];
        gram_schmidt(u, a, 3);  // u1={1,0,0}, u2={0,1,0}, u3={0,0,1}
        
        // Create linear combination v = 2*a1 + 3*a2 + 4*a3 = {9, 7, 4}
        double alpha[3] = {2.0, 3.0, 4.0};
        double v[3] = {
            alpha[0] * a[0][0] + alpha[1] * a[1][0] + alpha[2] * a[2][0],
            alpha[0] * a[0][1] + alpha[1] * a[1][1] + alpha[2] * a[2][1],
            alpha[0] * a[0][2] + alpha[1] * a[1][2] + alpha[2] * a[2][2]
        };
        
        // Find coefficients βᵢ by solving v = Σ βᵢuᵢ
        // Since u is orthogonal basis: βᵢ = (uᵢᵀv) / (uᵢᵀuᵢ)
        double beta[3];
        for (int i = 0; i < 3; i++) {
            beta[i] = dot(u[i], v, 3) / dot(u[i], u[i], 3);
        }
        
        // Compute v_reconstructed = Σ βᵢuᵢ
        double v_recon[3] = {0, 0, 0};
        for (int i = 0; i < 3; i++) {
            for (int d = 0; d < 3; d++) {
                v_recon[d] += beta[i] * u[i][d];
            }
        }
        
        // Check if v = v_reconstructed
        double diff = 0.0;
        for (int i = 0; i < 3; i++) {
            double d = v[i] - v_recon[i];
            diff += d * d;
        }
        diff = sqrt(diff);
        
        printf("\nTest 1.1(d): Linear combo for k=3 vectors\n");
        printf("  v = %.0f*a1 + %.0f*a2 + %.0f*a3 = ", alpha[0], alpha[1], alpha[2]);
        print_vec("v", v);
        printf("  Represented as %.2f*u1 + %.2f*u2 + %.2f*u3\n", beta[0], beta[1], beta[2]);
        print_vec("  v_recon", v_recon);
        printf("  ||v - v_recon|| = %.2e\n", diff);
        
        if (fabs(diff) < TOL) passed++;
    }
    
    // Test (e): Linearly dependent
    {
        total++;
        double a[][3] = {{1,0,0}, {0,1,0}, {2,3,0}};  // a3 = 2*a1 + 3*a2
        double u[3][3];
        gram_schmidt(u, a, 3);
        double norm = sqrt(dot(u[2], u[2], 3));
        printf("\nTest 1.1(e): a3 = 2*a1 + 3*a2\n");
        printf("  |u3| = %.2e (should be ~0)\n", norm);
        if (fabs(norm) < TOL) passed++;
    }
    
    printf("\n%s: %d/%d tests passed\n", 
           passed == total ? "✓ SUCCESS" : "✗ FAILED", passed, total);
    return passed == total ? 0 : 1;
}
