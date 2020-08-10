/* Program Gauss-Seidel
   Solution of a system of linear equations by Gauss-Seidel's
   iteration method. Assume that the coefficient matrix satisfies
   the condition of convergence.*/

#include<stdio.h>

#include<math.h>

#include<stdlib.h>

int main() {
  float a[10][10], b[10], x[10], xn[10], epp = 0.00001, sum;
  int i, j, n, flag;

  printf("Enter number of variables: ");
  scanf("%d", & n);
  printf("Enter the coefficients row-wise: ");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      scanf("%f", & a[i][j]);
  printf("Enter right hand vectors: ");
  for (i = 0; i < n; i++)
    scanf("%f", & b[i]);
  for (i = 0; i < n; i++)
    x[i] = 0; //initialize

  /* testing of diagonal dominance may be included here from 
     the program of Gauss-Jacobi's method */
  do {
    for (i = 0; i < n; i++) {
      sum = b[i];
      for (j = 0; j < n; j++) {
        if (j < i)
          sum -= a[i][j] * xn[j];
        else if (j > i)
          sum -= a[i][j] * x[j];
        xn[i] = sum / a[i][j];
      }
    }
    flag = 0; // indicates |x[i]-xn[i]|<epp for all i
    for (i = 0; i < n; i++)
      if (fabs(x[i] - xn[i]) > epp)
        flag = 1;
    if (flag == 1)
      for (i = 0; i < n; i++)
        x[i] = xn[i]; // reset x[i]	
  } while (flag == 1);
  printf("Solution is \n");
  for (i = 0; i < n; i++)
    printf("%8.5f ", xn[i]);
}