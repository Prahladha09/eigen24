#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define limit 1e-9
#define maxiter 100

void eigenvalues(double *A, double *V, int n) {
    int iterations = 0;
    double max;

    do {
        max = 0.0;
        int p, q;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (fabs(A[i * n + j]) > max) {
                    max = fabs(A[i * n + j]);
                    p = i;
                    q = j;
                }
            }
        }

        if (max < limit) {
            break;
        }
        double theta = 0.5 * atan2(2 * A[p * n + q], A[q * n + q] - A[p * n + p]);
        double c = cos(theta);
        double s = sin(theta);
        double App = A[p * n + p];
        double Aqq = A[q * n + q];
        double Apq = A[p * n + q];

        A[p * n + p] = c * c * App + s * s * Aqq - 2 * s * c * Apq;
        A[q * n + q] = s * s * App + c * c * Aqq + 2 * s * c * Apq;
        A[p * n + q] = A[q * n + p] = 0.0;
        for (int i = 0; i < n; i++) {
            if (i != p && i != q) {
                double Aip = A[i * n + p];
                double Aiq = A[i * n + q];
                A[i * n + p] = A[p * n + i] = c * Aip - s * Aiq;
                A[i * n + q] = A[q * n + i] = s * Aip + c * Aiq;
            }
            double Vip = V[i * n + p];
            double Viq = V[i * n + q];
            V[i * n + p] = c * Vip - s * Viq;
            V[i * n + q] = s * Vip + c * Viq;
        }
        iterations++;
    } while (max > limit && iterations < maxiter);

    if (iterations >= maxiter) {
        printf("Maximum iterations reached without full convergence.\n");
    }
}

void printMatrix(double *matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", matrix[i * n + j]);
        }
        printf("\n");
    }
}

int main() {
    int n;
    printf("Enter the dimension of the matrix: ");
    scanf("%d", &n);
    double *A = (double *)malloc(n * n * sizeof(double));
    double *V = (double *)malloc(n * n * sizeof(double));
    printf("Enter the elements of the %dx%d symmetric matrix:\n", n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &A[i * n + j]);
            if (i == j) {
                V[i * n + j] = 1.0;
            } else {
                V[i * n + j] = 0.0;
            }
        }
    }

    eigenvalues(A, V, n);

    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%lf ", A[i * n + i]);
    }
    printf("\n\nEigenvectors:\n");
    printMatrix(V, n);
    free(A);
    free(V);
    return 0;
}

