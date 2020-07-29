/* Programa para achar matriz superior e inferior */ 
#include<conio.h>

#include<stdio.h>

#include <stdbool.h>

#define clrscr() printf("\e[1;1H\e[2J")

int main() {
  bool lower;
  int rows, cols, r, c, matrix[10][10];
  clrscr(); /*Limpa a Tela*/
  printf("Se voce quiser a matriz triangular superior digite 0, se quiser a inferior digite qualquer numero diferente de 0: ");
  scanf("%d", & lower);
  printf("Digite o numero de linhas da matriz: ");
  scanf("%d", & rows);
  printf("\n");
  printf("Digite o numero de colunas da matriz: ");
  scanf("%d", & cols);
  printf("\n");
  printf("Digite os elementos da matriz: \n");
  for (r = 0; r < rows; r++) {
    for (c = 0; c < cols; c++) {
      printf("m[%d][%d] ", r + 1, c + 1);
      scanf("%d", & matrix[r][c]);
    }
  }
  if (lower) {
    printf("\n A matriz inferior e: ");
    for (r = 0; r < rows; r++) {
      printf("\n");
      for (c = 0; c < cols; c++) {
        if (r >= c) {
          printf("%d\t ", matrix[r][c]);
        } else {
          printf("0");
          printf("\t");
        }
      }
    }
  } else {
    printf("\n\n A matriz superior e: ");
    for (r = 0; r < rows; r++) {
      printf("\n");
      for (c = 0; c < cols; c++) {
        if (r > c) {
          printf("0");
          printf("\t");
        } else {
          printf("%d\t ", matrix[r][c]);

        }
      }
    }
  }

  getch();
  return 0;
}