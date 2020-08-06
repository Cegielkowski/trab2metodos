/*
TRABALHO 2 DISCIPLINA DE M�TODOS NUM�RICOS COMPUTACIONAIS
UNESP Bauru - 1 SEM 2020
AMANDA MEIRA
ARTHUR CIPOLARI 151022071
LUCAS CEGIELKOWSKI
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
                    printf("\nErro matriz nao definida positiva\n\n");
                    return 0;
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

void gaussSimples(int n, float A[n][n], float b[n]){

    int i, j, k;
    float m[n][n];
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

    printf("\n--x--\n");
    for(i=0;i<n;i++){
        printf("|%.2f|\n", b[i]);
    }

}

int main(void){

    int n, i, j, metodo=99;

	printf("Programa desenvolvido para o Trabalho 2 da\ndisciplina de Metodos Numericos Computacionais\nda UNESP Bauru - 1 Sem 2020\n\n\n");

    printf("Digite a ordem N da matriz => ");
    scanf("%d", &n);

    float A[n][n], b[n];

    while (metodo!=0){
        printf("\n\n=============MENU===============");
        printf("\n================================");
        printf("\n===  1-                      ===");
        printf("\n===  2-                      ===");
        printf("\n===  3-                      ===");
        printf("\n===  4-Decomposicao LU       ===");
        printf("\n===  5-Cholesky              ===");
        printf("\n===  6-Gauss Simples         ===");
        printf("\n===  0-Para Sair             ===");
        printf("\n================================");
        printf("\n================================");
        printf("\nEscolha: ");
        scanf ("%d", &metodo);
        if(metodo==0){
            printf ("\nPressione qualquer tecla para fechar...\n");
            exit(0);
        }
        switch (metodo){
        case 1:
            printf("Nao implementado");
            break;

        case 2:
            printf("Nao implementado");
            break;

        case 3:
            printf("Nao implementado");
            break;

        case 4:
            printf("\nMETODO DECOMPOSICAO LU\n\n");
            //Lendo Matriz A
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    printf("Digite o valor da posicao A[%d][%d]=",i+1,j+1);
                    scanf("%f",&A[i][j]);
                }
            }

            printf("\nEntre com os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

            //Print matriz e termo independente
            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);
            printf("\n\n--b--\n");
            for(i=0;i<n;i++){
                printf("|%.2f|\n", b[i]);
            }

            decomposicaoLU(n, A, b);
            break;

        case 5:
            printf("\nMETODO CHOLESKY ESCOLHIDO\n\n");

            //Lendo Matriz A
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    printf("Digite o valor da posicao A[%d][%d]=",i+1,j+1);
                    scanf("%f",&A[i][j]);
                }
            }

            printf("\nEntre com os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

            //Print matriz e termo independente
            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);
            printf("\n\n--b--\n");
            for(i=0;i<n;i++){
                printf("|%.2f|\n", b[i]);
            }

            cholesky(n, A, b);

            break;
        case 6:
            printf("\nMETODO GAUSS SIMPLES ESCOLHIDO\n\n");

            //Lendo Matriz A
            for(i=0;i<n;i++){
                for(j=0;j<n;j++){
                    printf("Digite o valor da posicao A[%d][%d]=",i+1,j+1);
                    scanf("%f",&A[i][j]);
                }
            }

            printf("\nEntre com os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

            //Print matriz e termo independente
            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);
            printf("\n\n--b--\n");
            for(i=0;i<n;i++){
                printf("|%.2f|\n", b[i]);
            }

            gaussSimples(n, A, b);

            break;

        default:
            printf("\nOpcao Invalida!\nEscolha novamente...\n\n");
        }
    }
    return 0;
}


