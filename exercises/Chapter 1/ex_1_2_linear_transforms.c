/*
 * Exercise 1.2: Linear Transforms
 * 
 * Proves relationships between parallelogram area, determinant, and linear transforms.
 * Demonstrates that det(A) represents area scaling under transformation.
 */

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define TOL 1e-10

// Dot product
double dot(const double *a, const double *b) {
    return a[0]*b[0] + a[1]*b[1];
}

// 2x2 determinant: det(A) = a11*a22 - a12*a21
double det2(double a11, double a12, double a21, double a22) {
    return a11*a22 - a12*a21;
}

// Area formula from (a): S² = (a₂ᵀa₂)(a₁ᵀa₁) - (a₂ᵀa₁)²
double area_gram_schmidt(const double *a1, const double *a2) {
    double a1_dot_a1 = dot(a1, a1);
    double a2_dot_a2 = dot(a2, a2);
    double a1_dot_a2 = dot(a1, a2);
    return sqrt(a1_dot_a1 * a2_dot_a2 - a1_dot_a2 * a1_dot_a2);
}

// Area via determinant from (b): S = |det A|
double area_determinant(const double *a1, const double *a2) {
    return fabs(det2(a1[0], a2[0], a1[1], a2[1]));
}

// Transform rectangle by matrix A
void transform_rect(double *y, const double *x, 
                    double a11, double a12, double a21, double a22) {
    y[0] = a11*x[0] + a12*x[1];
    y[1] = a21*x[0] + a22*x[1];
}

// Test (a) and (b): Area formulas are equivalent
bool test_1_2_ab(void) {
    printf("\n=== Test 1.2(a,b): Area of parallelogram ===\n");
    
    // Test case 1: Simple vectors
    double a1[] = {3.0, 0.0};
    double a2[] = {1.0, 2.0};
    
    double s_gram = area_gram_schmidt(a1, a2);
    double s_det = area_determinant(a1, a2);
    
    printf("Vectors: a1=[%.1f,%.1f], a2=[%.1f,%.1f]\n", 
           a1[0], a1[1], a2[0], a2[1]);
    printf("  Area (Gram-Schmidt): %.6f\n", s_gram);
    printf("  Area (determinant):  %.6f\n", s_det);
    printf("  Difference: %.2e\n", fabs(s_gram - s_det));
    
    if (fabs(s_gram - s_det) > TOL) return false;
    
    // Test case 2: Another example
    double b1[] = {2.0, 1.0};
    double b2[] = {1.0, 3.0};
    
    s_gram = area_gram_schmidt(b1, b2);
    s_det = area_determinant(b1, b2);
    
    printf("\nVectors: b1=[%.1f,%.1f], b2=[%.1f,%.1f]\n", 
           b1[0], b1[1], b2[0], b2[1]);
    printf("  Area (Gram-Schmidt): %.6f\n", s_gram);
    printf("  Area (determinant):  %.6f\n", s_det);
    printf("  Difference: %.2e\n", fabs(s_gram - s_det));
    
    return fabs(s_gram - s_det) < TOL;
}

// Test (c): Linear transform scales area by |det A|
bool test_1_2_c(void) {
    printf("\n=== Test 1.2(c): Area scaling under transform ===\n");
    
    // Rectangle Ux: [0,Δ1] × [0,Δ2]
    double delta1 = 2.0, delta2 = 3.0;
    double area_x = delta1 * delta2;
    
    // Matrix A
    double a11 = 2.0, a12 = 1.0;
    double a21 = 1.0, a22 = 3.0;
    double det_a = det2(a11, a12, a21, a22);
    
    // Transform corners of rectangle
    double corners[4][2] = {{0,0}, {delta1,0}, {0,delta2}, {delta1,delta2}};
    double transformed[4][2];
    
    for (int i = 0; i < 4; i++) {
        transform_rect(transformed[i], corners[i], a11, a12, a21, a22);
    }
    
    // Transformed rectangle spans parallelogram with vectors Δ1*a1 and Δ2*a2
    double scaled_a1[] = {delta1*a11, delta1*a21};
    double scaled_a2[] = {delta2*a12, delta2*a22};
    double area_y = area_determinant(scaled_a1, scaled_a2);
    
    double expected_area = area_x * fabs(det_a);
    
    printf("Rectangle Ux: [0,%.1f]×[0,%.1f], area=%.1f\n", delta1, delta2, area_x);
    printf("Matrix A: [[%.1f,%.1f],[%.1f,%.1f]], det(A)=%.1f\n", 
           a11, a12, a21, a22, det_a);
    printf("  Parallelogram Uy area: %.6f\n", area_y);
    printf("  Expected (area_x × |det A|): %.6f\n", expected_area);
    printf("  Difference: %.2e\n", fabs(area_y - expected_area));
    
    return fabs(area_y - expected_area) < TOL;
}

// Test (d): Change of variables intuition
bool test_1_2_d(void) {
    printf("\n=== Test 1.2(d): Change of variables ===\n");
    printf("Intuition: ∫_Uy f(y)dy = ∫_Ux f(Ax)|det A|dx\n\n");
    
    // Simple test: integrate constant function f(x)=1 over rectangle
    // Should give area of transformed region
    
    double delta1 = 2.0, delta2 = 1.5;
    double area_x = delta1 * delta2;
    
    double a11 = 3.0, a12 = 0.5;
    double a21 = 0.5, a22 = 2.0;
    double det_a = det2(a11, a12, a21, a22);
    
    // For f(x)=1: ∫_Ux f(Ax)|det A|dx = |det A| × area(Ux)
    double integral_value = fabs(det_a) * area_x;
    
    // This should equal area(Uy)
    double scaled_a1[] = {delta1*a11, delta1*a21};
    double scaled_a2[] = {delta2*a12, delta2*a22};
    double area_y = area_determinant(scaled_a1, scaled_a2);
    
    printf("For constant f(x)=1:\n");
    printf("  ∫_Ux |det A| dx = |det A| × area(Ux) = %.1f × %.1f = %.6f\n", 
           fabs(det_a), area_x, integral_value);
    printf("  Area(Uy) = %.6f\n", area_y);
    printf("  Difference: %.2e\n", fabs(integral_value - area_y));
    printf("\nThe |det A| term compensates for volume change under transformation.\n");
    
    return fabs(integral_value - area_y) < TOL;
}

// Edge case: degenerate matrix (det=0)
bool test_edge_degenerate(void) {
    printf("\n=== Edge case: Degenerate matrix ===\n");
    
    // Linearly dependent vectors: a2 = 2*a1
    double a1[] = {1.0, 2.0};
    double a2[] = {2.0, 4.0};
    
    double area = area_determinant(a1, a2);
    double det = det2(a1[0], a2[0], a1[1], a2[1]);
    
    printf("Vectors: a1=[%.1f,%.1f], a2=[%.1f,%.1f]\n", 
           a1[0], a1[1], a2[0], a2[1]);
    printf("  det(A) = %.6f\n", det);
    printf("  Area = %.6f\n", area);
    printf("  Expected: 0 (degenerate)\n");
    
    return fabs(area) < TOL && fabs(det) < TOL;
}

int main(void) {
    int passed = 0, total = 0;
    
    total++; if (test_1_2_ab()) passed++;
    total++; if (test_1_2_c()) passed++;
    total++; if (test_1_2_d()) passed++;
    total++; if (test_edge_degenerate()) passed++;
    
    printf("\n%s: %d/%d tests passed\n", 
           passed == total ? "✓ SUCCESS" : "✗ FAILED", passed, total);
    return passed == total ? 0 : 1;
}
