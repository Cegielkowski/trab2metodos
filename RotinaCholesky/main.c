#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

void mostraMatrizT(int n, float mat[n][n]){
    int i,j;
    for(i=0;i<n;i++){
        printf("|");
        for(j=0;j<n;j++){
            printf(" %.2f", mat[j][i]);
        }
        printf(" |\n");
    }
}

float cholesky(int n, float A[n][n], float b[n]){
    int i, j, k;
    float G[n][n], X[n], Y[n], somat;

    for (k = 0; k < n; k++){
        for (i = 0; i <= k; i++){
            somat = 0;
            for (j = 0; j < i; j++){
                somat += G[i][j] * G[k][j];
            }

            if (i == k){
                if(A[i][i] - somat <= 0){
                    printf("Erro matriz nao definida positiva");
                    exit(0);
                }
                G[i][i] = sqrt(A[i][i] - somat);
            }
            else
            G[k][i] = 1.0 / G[i][i] * (A[k][i] - somat);
        }
    }
    printf("\n----Matriz G----\n");
    mostraMatriz(n, G);
    printf("\n----Matriz Gt----\n");
    mostraMatrizT(n, G);

    //Resolvendo G.y = b
    printf("\n--y--\n");
    for(i = 0; i < n; ++i) {
        somat = 0;

        for(j = 0; j < i; ++j){
            somat += G[i][j] * Y[j];
        }
        Y[i] = (b[i] - somat) / G[i][i];
    }
    for(i=0;i<n;i++){
        printf("|%.2f|\n", Y[i]);
    }

    //Resolvendo Gt.x = y
    printf("\n--x--\n");
    for (i=n-1;i>=0;i--) {
        somat=0;
        for(j=i+1;j<n;j++){
            somat += G[j][i] * X[j];
        }
        X[i] = (Y[i] - somat) / G[i][i];
    }
    for(i=0;i<n;i++){
        printf("|%.2f|\n", X[i]);
    }

}

int main(){
    int n, i, j;

    printf("Digite a ordem N da matriz => ");
    scanf("%d", &n);

    float A[n][n], b[n];
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
    cholesky(n, A, b);

    return 0;
}
