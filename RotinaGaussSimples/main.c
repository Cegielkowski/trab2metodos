#include <stdio.h>
#include <stdlib.h>

void mostraMatriz(int n, float mat[n][n]){
    int i,j;
    for(i=0;i<n;i++){
        printf("| ");
        for(j=0;j<n;j++){
            printf("%5.2f", mat[i][j]);
        }
        printf(" |\n");
    }
}

void gaussSimples(int n, float A[n][n], float b[n]){

    int i, j, k;
    float somat, m[n][n], X[n];
	for (k = 0; k < n - 1; k++) {
		for (i = k + 1; i < n; i++) {
			m[i][k]= - (A[i][k]/A[k][k]);
			for (j = k; j < n; j++) {
				A[i][j] = A[i][j] + (m[i][k] * A[k][j]);
			}
			b[i] = b[i] + (m[i][k] * b[k]);
		}
	}

    printf("\n\n----Matriz Escalonada----\n");
    mostraMatriz(n, A);

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
        printf("|%5.2f|\n", X[i]);
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

	gaussSimples(n, A, b);

	return 0;
}
