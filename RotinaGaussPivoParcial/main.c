#include <stdio.h>
#include <stdlib.h>

void mostraMatriz(int n, float mat[n][n]){
    int i,j;
    for(i=0;i<n;i++){
        printf("|");
        for(j=0;j<n;j++){
            printf(" %.2f", mat[i][j]);
        }
        printf(" |\n");
    }
}

void gaussPivo(int n, float A[n][n], float b[n]){

    int i, j, k, l;
    float max, aux, auxb, pivo, somat, X[n], m[n][n];

    for (k = 0; k < n-1; k++){
        max = abs(A[k][k]);
        l = k;
        // Procura maior pivo
        for (i = k+1; i < n; i++){
          if (abs(A[i][k]) > max){
        max = abs(A[i][k]);
        l = i;
          }
        }
        // Verifica se o maior é o pivo
        if (l != k){
          // Troca linha atual pela linha do maior pivo
          for (j = k; j <= n; j++){
            aux = A[k][j];
            auxb = b[k];
            A[k][j] = A[l][j];
            A[l][j] = aux;
            b[k] = b[l];
            b[l] = auxb;
          }
        }
        // Metodo de Gauss Após pivotagem
        for (i = k+1; i < n; i++){
          m[i][k]= - (A[i][k]/A[k][k]);
          A[i][k] = 0;
          for (j = k+1; j <= n; j++){
            A[i][j] = A[i][j] + m[i][k] * A[k][j];
          }
          b[i] = b[i] + (m[i][k] * b[k]);
        }
    }

    printf("\n\n----Matriz Escalonada----\n");
    mostraMatriz(n, A);

    printf("\n--b--\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", b[i]);
    }

    // SUBST REGRESSIVA
    printf("\n--x--\n");
    X[n-1] = b[n-1]/A[n-1][n-1];
    for(i=n-2;i>=0;i--){
        somat = 0;
        for(j=i+1;j<n;j++){
            somat += A[i][j]*X[j];
        }
        X[i] = (b[i] - somat )/A[i][i];
    }
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", X[i]);
    }

}


int main (void) {


    int n, i, j, k, metodo=99;

    printf("Digite a ordem N da matriz => ");
    scanf("%d", &n);

    float A[n][n], b[n], X[n];

    //Lendo Matriz A
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("Digite o valor da posicao A[%d][%d]=",i+1,j+1);
            scanf("%f",&A[i][j]);
        }
    }
    printf("----Matriz A----\n");
    mostraMatriz(n, A);

    printf("\n Entre com os termos independentes (b):\n");
    for(i=0; i<n;i++){
        printf("b[%d]:", i+1);
        scanf("%f", &b[i]);
    }
    printf("\n\n--b--\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", b[i]);
    }

	gaussPivo(n, A, b);

	return 0;
}
