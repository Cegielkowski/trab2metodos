#include <stdio.h>
#include <stdlib.h>

void mostraMatriz(int n, float mat[n][n]){
    int i,j;
    for(i=0;i<n;i++){
        printf("|");
        for(j=0;j<n;j++){
            printf(" %.2f",mat[i][j]);
        }
        printf(" |\n");
    }
}

void decomposicaoLU(int n, float mat[n][n], float b[n]){
    int i,j;
    float M[n][n],aux[n][n], Y[n], X[n], somat;

    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            M[i][j]=0;
            aux[i][j]=mat[i][j];
        }
    }
    for(j=0;j<n;j++){
        for(i=j+1;i<n;i++){
            if(aux[i][j]!=0){
                M[i][j]=aux[i][j]/aux[j][j];
                int c=0;
                for(c=j;c<n;c++){
                    aux[i][c]=aux[i][c]+aux[j][c]*(-1*(M[i][j]));
                }
            }
        }

    }
    printf("\n\nRESULTADO\n");
    for(i=0;i<n;i++) M[i][i]=1;
    printf("\n----Matriz L----\n");
    mostraMatriz(n, M);

    printf("\n----Matriz U----\n");
    mostraMatriz(n, aux);


    // SUBST PROGRESSIVA
    printf("\n--y--\n");
    Y[0]=b[0]/M[0][0];
    for(i=1;i<n;i++){
        somat = 0;
        for(j=0;j<i;j++){
            somat += M[i][j]*Y[j];
        }
        Y[i] = (b[i] - somat )/M[i][i];
    }
    for(i=0;i<n;i++){
        printf("|%.2f|\n", Y[i]);
    }



    // SUBST REGRESSIVA
    printf("\n--x--\n");
    X[n-1] = Y[n-1]/aux[n-1][n-1];
    for(i=n-2;i>=0;i--){
        somat = 0;
        for(j=i+1;j<n;j++){
            somat += aux[i][j]*X[j];
        }
        X[i] = (Y[i] - somat )/aux[i][i];
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

    decomposicaoLU(n, A, b);

    return 0;
}
