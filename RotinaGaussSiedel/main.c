#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int main() {
  float a[10][10], b[10], x[10], xn[10], epp, sum;
  int i, j, n, flag, maxIt, itAtual=0;

  printf("precisao desejada(e): ");
  scanf("%f", &epp);

  printf("numero maximo de iteracoes: ");
  scanf("%d", &maxIt);

  printf("ordem do sistema: ");
  scanf("%d", &n);
  printf("matriz dos coeficientes: ");
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
	    scanf("%f", & a[i][j]);
	}
  printf("vetor dos termos independentes: ");
  for (i = 0; i < n; i++)
    scanf("%f", & b[i]);
  for (i = 0; i < n; i++)
    x[i] = 0; 

  do {
  	itAtual = itAtual + 1;
    for (i = 0; i < n; i++) {
      sum = b[i];
      for (j = 0; j < n; j++) {
        if (j < i) {
          sum -= a[i][j] * xn[j];
        }
        else if (j > i) {
          sum -= a[i][j] * x[j];
        }
        xn[i] = sum / a[i][j];
      }
    }
    flag = 0;
    for (i = 0; i < n; i++) {
      if (fabs(x[i] - xn[i]) > epp) {
        flag = 1;
      }
    }
    if (flag == 1){ 
      for (i = 0; i < n; i++) {
        x[i] = xn[i]; 
      }
    }
    if (itAtual == maxIt) {
      printf("num interacoes maximo atingido ");
  	  printf("%d ", itAtual);
      printf("/n A solucao atingida ate esta iteracao eh: \n");	 
	  for (i = 0; i < n; i++)
   		printf("%8.5f ", xn[i]);

      return 0;
    }
  } while (flag == 1);

  printf("A solucao e \n");
  for (i = 0; i < n; i++)
    printf("%8.5f ", xn[i]);
}
