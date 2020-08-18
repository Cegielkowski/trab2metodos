/*
TRABALHO 2 DISCIPLINA DE M?TODOS NUM?RICOS COMPUTACIONAIS
UNESP Bauru - 1 SEM 2020
AMANDA MEIRA        151020191
ARTHUR CIPOLARI     151022071
LUCAS CEGIELKOWSKI  161025978
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <conio.h>
#define NMAX 100

int lerMatriz(float A[NMAX][NMAX]){
    int i, j, n;
    printf("\nDigite a ordem N da matriz\n");
	scanf("%d", &n);
    //Lendo Matriz A
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("Digite o valor da posicao A[%d][%d]=",i+1,j+1);
            scanf("%f",&A[i][j]);
        }
    }
    return n;
}

void mostraMatriz(int n, float mat[NMAX][NMAX]){
    int i,j;
    for(i=0;i<n;i++){
        printf("|");
        for(j=0;j<n;j++){
            printf(" %5.2f",mat[i][j]);
        }
        printf(" |\n");
    }
}

void calculaMatrizT(int n, float G[NMAX][NMAX], float Gt[NMAX][NMAX]){
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            Gt[i][j] = G[j][i];
        }
    }
}

void minor(float b[NMAX][NMAX],float a[NMAX][NMAX],int i,int n){
	int j,l,h=0,k=0;
	for(l=1;l<n;l++)
		for( j=0;j<n;j++){
			if(j == i)
				continue;
			b[h][k] = a[l][j];
			k++;
			if(k == (n-1)){
				h++;
				k=0;
			}
		}
}

float determinante(int n, float a[NMAX][NMAX]){
	int i;
	float b[NMAX][NMAX], sum=0;
	if (n == 1)
return a[0][0];
	else if(n == 2)
return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
	else
		for(i=0;i<n;i++){
			minor(b,a,i,n);
			sum = (float) (sum+a[0][i]*pow(-1,i)*determinante((n-1), b));
		}
return sum;
}

float triangularSuperior(int n, float A[NMAX][NMAX], float b[NMAX], float X[NMAX]){
    int i, j;
    float somat;

    X[n-1] = b[n-1]/A[n-1][n-1];
    for(i=n-2;i>=0;i--){
        somat = 0;
        for(j=i+1;j<n;j++){
            somat += A[i][j]*X[j];
        }
        X[i] = (b[i] - somat )/A[i][i];
    }
}

float triangularInferior(int n, float A[NMAX][NMAX], float b[NMAX], float X[NMAX]){
    int i, j;
    float somat;

    for(i = 0; i < n; ++i) {
        somat = 0;
        for(j = 0; j < i; ++j){
            somat += A[i][j] * X[j];
        }
        X[i] = (b[i] - somat) / A[i][i];
    }
}

void decomposicaoLU(int n, float mat[NMAX][NMAX], float b[NMAX]){
    int i,j;
    float M[NMAX][NMAX],aux[NMAX][NMAX], Y[NMAX], X[NMAX], somat;

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

    triangularInferior(n, M, b, Y);
    printf("\n----y----\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", Y[i]);
    }

    triangularSuperior(n, aux, Y, X);
    printf("\n----x----\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", X[i]);
    }

}

float cholesky(int n, float A[NMAX][NMAX], float b[NMAX]){
    int i, j, k;
    float G[NMAX][NMAX], Gt[NMAX][NMAX], X[NMAX], Y[NMAX], somat;
    float *Xt;

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
    calculaMatrizT(n, G, Gt);
    printf("\n----Matriz Gt----\n");
    mostraMatriz(n, Gt);

    //Resolvendo G.y = b
    //calcula substituicao progressiva
    triangularInferior(n, G, b, Y);
    printf("\n----y----\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", Y[i]);
    }

    //Resolvendo Gt.x = y
    //calcula substituicao regressiva
    triangularSuperior(n, Gt, Y, X);
    printf("\n----x----\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", X[i]);
    }

}

void gaussSimples(int n, float A[NMAX][NMAX], float b[NMAX]){

    int i, j, k;
    float somat, m[NMAX][NMAX], X[NMAX];
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

    triangularSuperior(n, A, b, X);
    printf("\n----x----\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", X[i]);
    }

}

void gaussPivoteamentoParcial(int n, float A[NMAX][NMAX], float b[NMAX]){

    int i, j, k, l;
    float max, aux, auxb, pivo, somat, X[NMAX], m[NMAX][NMAX];

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

    printf("\n----b----\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", b[i]);
    }

    triangularSuperior(n, A, b, X);
    printf("\n----x----\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", X[i]);
    }

}

void matrizInversa(int n, float a[NMAX][NMAX]) {
    float b[NMAX][NMAX], c[NMAX][NMAX], d[NMAX][NMAX], deter;
	int l,h,m,k,i,j;
	deter = (float) determinante(n, a);

	printf("\n\n Determinante da matriz: %.4f ",deter);

	if(deter == 0)
		printf("\n Matriz inversa da matriz digitada nao e possivel\n");
	else if(n == 1)
		d[0][0] = 1;
	else{
        for (h=0;h<n;h++){
            for (l=0;l<n;l++){
                m=0;
                k=0;
                for (i=0;i<n;i++)
                    for (j=0;j<n;j++)
                        if (i != h && j != l){
                            b[m][k]=a[i][j];
                            if (k<(n-2))
                                k++;
                            else{
                                k=0;
                                m++;
                            }
                        }
                c[h][l] = pow(-1,(h+l))*determinante((n-1), b);	// c = cofator da Matriz
            }
        }
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
                b[i][j] = c[j][i];
            }
        }
        for (i=0;i<n;i++){
            for (j=0;j<n;j++){
                d[i][j] = b[i][j]/deter;	// d =  Matriz Inversa
            }
        }
	}
    printf("\n----Matriz Inversa----\n");
	mostraMatriz(n, d);
}

void JacobiRichardson(int n, float A[NMAX][NMAX], float b[NMAX], float aproximacao[NMAX], int iteracoes, float erro){
	float  dinversa[NMAX][NMAX], matrizR[NMAX][NMAX], temporaria[NMAX], matrizResult[NMAX];
	int linha, coluna;

	void multiplicacao(float matrizA[NMAX][NMAX], float matrizB[NMAX])
	{
		int i, j;
		for (i = 0; i < n; i++)
		{
			matrizResult[i]=0;
			for (j = 0; j < n; j++)
				matrizResult[i]=matrizResult[i]+matrizA[i][j]*matrizB[j];
		}
	}

	for (linha = 0; linha < n; linha++)
		for (coluna = 0; coluna < n; coluna++)
		{
			if (linha == coluna)
				dinversa[linha][coluna] = 1/A[linha][coluna];
			else
				dinversa[linha][coluna] = 0;
		}
		for (linha = 0; linha < n; linha ++)
			for (coluna = 0; coluna < n; coluna++)
			{
				if (linha == coluna)
					matrizR[linha][coluna] = 0;
				else
				if (linha != coluna)
					matrizR[linha][coluna] = A[linha][coluna];
			}

	int i = 1;
	int j;
	while (i <= iteracoes)
	{
		multiplicacao(matrizR,aproximacao);
		for (linha = 0; linha < n; linha ++)
			temporaria[linha] = b[linha] - matrizResult[linha];
		multiplicacao(dinversa, temporaria);
		for (j = 0; j < n; j++)
			aproximacao[j] = matrizResult[j];
		printf("Os valores de x apos a interacao %d sao\n",i);
		for (linha = 0; linha < n; linha++)
			printf("%.3f\n", aproximacao[linha]);
		i++;
	}
	printf("\nNumero maximo de iteracoes atingido!");
}

void gaussSeidel(int n, float a[NMAX][NMAX], float b[NMAX], float x[NMAX], float epp, int maxIt) {
    float xn[NMAX], sum, t, e;
    int i, j, itAtual;

        printf("-----------------Tabela--Iteracoes----------------\n");
        for(itAtual=1;itAtual<=maxIt;itAtual++){
            for(i=0;i<n;i++){
                sum=0;
                for(j=0;j<n;j++)
                if(j!=i){
                  sum+=a[i][j]*xn[j];
                }
                t=(b[i]-sum)/a[i][i];
                e=fabs(xn[i]-t);
                xn[i]=t;
            }
            if(e<epp){
                printf("\nConverge em %d iteracoes\n", itAtual);
                printf("\nErro: %f < %f\n", e, epp);
                for(i=0;i<n;i++)
                printf("a[%3d]=%7.4f\n", i+1,xn[i]);
                return;
            }

            printf(" %5d\t",itAtual);
            for(i=0;i<n;i++)
            printf(" %9.4f\t",xn[i]);
            printf("\n");
            if(itAtual>=maxIt){
                printf("\nMaximo de %d iteracoes atingido\n", itAtual);
                printf("\nErro alcancado: %f \n", e);
                for(i=0;i<n;i++)
                printf("x[%3d]=%7.4f\n", i+1,xn[i]);
                return;
            }
        }
}

int main(void){

    int n, i, j, maxIt, metodo=99;

	printf("Programa desenvolvido para o Trabalho 2 da\ndisciplina de Metodos Numericos Computacionais\nda UNESP Bauru - 1 Sem 2020\n\n\n");

	// Definindo um tamanho fixo para a matriz no inicio do programa
	// Apenas para dispensar uso de alocacao dinamica e/ou um 'nMAX' pre-definido
	// Consumindo recursos sem necessidade
    //printf("Digite a ordem N da matriz => ");
    //scanf("%d", &n);

    float A[NMAX][NMAX], b[NMAX], X[NMAX],aproximacao[NMAX], erro;

    while (metodo!=0){
        printf("\n\n=============MENU===============");
        printf("\n================================");
        printf("\n===  1-Determinante          ===");
        printf("\n===  2-Triangular Inferior   ===");
        printf("\n===  3-Triangular Superior   ===");
        printf("\n===  4-Decomposicao LU       ===");
        printf("\n===  5-Cholesky              ===");
        printf("\n===  6-Gauss Simples         ===");
        printf("\n===  7-Gauss Pivo Parcial    ===");
        printf("\n===  8-Matriz Inversa        ===");
        printf("\n===  9-Jacobi Richardson     ===");
        printf("\n===  10-Gauss Seidel         ===");
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
            printf("\nMETODO DETERMINANTE ESCOLHIDO\n\n");

            n = lerMatriz(A);

            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);

            printf("\nA determinante da matriz e': %5.2f\n", determinante(n, A));

            break;

        case 2:
			printf("\nMETODO SISTEMA TRIANGULAR INFERIOR ESCOLHIDO\n\n");

            n = lerMatriz(A);

           if(A[1][n-1] != 0 ){
                printf("\n Matriz nao e' triangular inferior\n");
                break;
            }

            printf("\nDigite os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

            //Print matriz e termo independente
            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);
            printf("\n\n----b----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", b[i]);
            }

            triangularInferior(n, A, b, X);
            printf("\n----x----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", X[i]);
            }
            break;

        case 3:
			printf("\nMETODO SISTEMA TRIANGULAR SUPERIOR ESCOLHIDO\n\n");

            n = lerMatriz(A);

            if(A[n-1][1] != 0 ){
                printf("\n Matriz nao e' triangular superior\n");
                break;
            }


            printf("\nDigite os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

            //Print matriz e termo independente
            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);
            printf("\n\n----b----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", b[i]);
            }

            triangularSuperior(n, A, b, X);
            printf("\n----x----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", X[i]);
            }
            break;

        case 4:
            printf("\nMETODO DECOMPOSICAO LU\n\n");

            n = lerMatriz(A);

            printf("\nDigite os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

            //Print matriz e termo independente
            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);
            printf("\n\n----b----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", b[i]);
            }

            decomposicaoLU(n, A, b);
            break;

        case 5:
            printf("\nMETODO CHOLESKY ESCOLHIDO\n\n");

            n = lerMatriz(A);

            printf("\nDigite os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

            //Print matriz e termo independente
            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);
            printf("\n\n----b----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", b[i]);
            }

            cholesky(n, A, b);

            break;
        case 6:
            printf("\nMETODO GAUSS SIMPLES ESCOLHIDO\n\n");

            n = lerMatriz(A);

            printf("\nDigite os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

            //Print matriz e termo independente
            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);
            printf("\n\n----b----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", b[i]);
            }

            gaussSimples(n, A, b);

            break;

        case 7:
            printf("\nMETODO GAUSS PIVOTEAMENTO PARCIAL ESCOLHIDO\n\n");

            n = lerMatriz(A);

            printf("\nDigite os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

            //Print matriz e termo independente
            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);
            printf("\n\n----b----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", b[i]);
            }

            gaussPivoteamentoParcial(n, A, b);
            break;

        case 8:
            printf("\nMETODO MATRIZ INVERSA ESCOLHIDO\n\n");

            n = lerMatriz(A);

            printf("\n----Matriz A----\n");
            mostraMatriz(n, A);

            matrizInversa(n, A);
            break;

        case 9:
        	printf("\nMETODO JACOBI RICHARDSON ESCOLHIDO\n\n");

			//Lendo Matriz A
            n = lerMatriz(A);

            printf("\nDigite os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

			printf("\nDigite a aproximacao inicial (x0):\n");
            for(i=0; i<n;i++){
                printf("x0[%d]:", i+1);
                scanf("%f", &aproximacao[i]);
            }

            printf("\nDigite o numero maximo de iteracoes: \n");
			scanf("%d", &maxIt);

			printf("\nDigite a precisao (e): \n");
			scanf("%f", &erro);

			JacobiRichardson(n, A, b, aproximacao, maxIt, erro);
			break;

        case 10:
            printf("\nMETODO GAUSS SEIDEL ESCOLHIDO\n\n");

			//Lendo Matriz A
            n = lerMatriz(A);

            printf("\nDigite os termos independentes (b):\n");
            for(i=0; i<n;i++){
                printf("b[%d]:", i+1);
                scanf("%f", &b[i]);
            }

			printf("\nDigite a aproximacao inicial (x0):\n");
            for(i=0; i<n;i++){
                printf("x0[%d]:", i+1);
                scanf("%f", &aproximacao[i]);
            }

            printf("\nDigite o numero maximo de iteracoes: \n");
			scanf("%d", &maxIt);

			printf("\nDigite a precisao (e): \n");
			scanf("%f", &erro);

            gaussSeidel(n, A, b, aproximacao, erro, maxIt);
            break;

        default:
            printf("\nOpcao Invalida!\nEscolha novamente...\n\n");
        }
    }

    return 0;
}
