/*
TRABALHO 2 DISCIPLINA DE M?TODOS NUM?RICOS COMPUTACIONAIS
UNESP Bauru - 1 SEM 2020
AMANDA MEIRA
ARTHUR CIPOLARI 151022071
LUCAS CEGIELKOWSKI 161025978
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <conio.h>

int cin(float a[100][100]){
	int i,j,n;
	printf("\n Digite o tamanho da matriz N*N : ");
	scanf("%d",&n);
	printf("\n--------------------------\n");
	for(i=0;i<n;i++)
		for(j=0;j<n;j++){
			printf(" Matriz[%d][%d] : ",i+1,j+1);
			scanf("%f",&a[i][j]);
		}
	printf("\n----------------------------------------------------\n");
return n;
}

//-----------------------------------------------------
// Mostra a matriz: 
void cout(float a[100][100],int n,int show){
	int i,j;
	if(show == 1)
		for(i=0;i < n;i++){
			for(j=0;j < n;j++)
				printf(" %.2f \t",a[i][j]);
			printf("\n");
		}
	else if(show == 2){
		printf("\n\n A Matriz inversa eh : \n\n");
		for (i=0;i<n;i++){
			for (j=0;j<n;j++)
				printf(" %.4f \t",a[i][j]);
			printf("\n");
		}
	}
}

void minor(float b[100][100],float a[100][100],int i,int n){
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

//---------------------------------------------------
//	calculate determinte of matrix
float det(float a[100][100],int n){
	int i;
	float b[100][100],sum=0;
	if (n == 1)
return a[0][0];
	else if(n == 2)
return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
	else
		for(i=0;i<n;i++){
			minor(b,a,i,n);	
			sum = (float) (sum+a[0][i]*pow(-1,i)*det(b,(n-1)));	
		}
return sum;
}

//---------------------------------------------------
//	calcula a matriz transposta
void transpose(float c[100][100],float d[100][100],float n,float det){
	int i,j;
	float b[100][100];
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			b[i][j] = c[j][i];
	for (i=0;i<n;i++)
		for (j=0;j<n;j++)
			d[i][j] = b[i][j]/det;	// array d[][] =  matriz inversa
}

//---------------------------------------------------
//	calcula cofator da matriz
void cofactor(float a[100][100],float d[100][100],float n,float determinte){
	float b[100][100],c[100][100];
	int l,h,m,k,i,j;
	for (h=0;h<n;h++)
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
			c[h][l] = pow(-1,(h+l))*det(b,(n-1));	// c = cofactor Matriz
		}
	transpose(c,d,n,determinte);	
}

//---------------------------------------------------
//	calcula o inveso da matriz
void inverse(float a[100][100],float d[100][100],int n,float det){
	if(det == 0)
		printf("\n Matriz inversa da matriz digitada nao e possivel\n");
	else if(n == 1)
		d[0][0] = 1;
	else
		cofactor(a,d,n,det);	
}

//---------------------------------------------------
void matrizInversa() {
	int i,j,n;
	float a[100][100],d[100][100],deter;
	printf("\n Programa para achar o inveso da matrix\n\n"); 
	n = cin(a);	
	int print_matrix = 1;
	cout(a,n,print_matrix);	
	deter = (float) det(a,n);
		printf("----------------------------------------------------\n");
		printf("\n\n Determinante da matriz: %.4f ",deter);
		printf("\n\n-----------------------\n");
	inverse(a,d,n,deter);	
	int print_inverse = 2;
	cout(d,n,print_inverse);
		printf("\n\n==============================* Fim *==============================\n");
	getchar();    
    return 0;
}
 
void gaussSiedel() {
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

void superiorInferior() {
  bool lower;
  int rows, cols, r, c, matrix[10][10];

  printf(" Se voce quiser a matriz triangular superior digite 0, se quiser a inferior digite qualquer numero diferente de 0: ");
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
}

void determinante()	{
	int m = 0;
    double **a = 0;    
    int i = 0, j = 0, k = 0;	
    double factor = 0;	
    double temp = 0;	
    int count = 0;	

    printf("dimensao => ");
    scanf("%d", &m);

    a = (double **) calloc(m, sizeof(double *));
    for(i = 0; i < m; i++)
    {
        a[i] = (double *) calloc(m, sizeof(double));
    }

    printf("\n\nEntre com o conteudo da matriz\n\n");
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < m; j++)
        {
            printf("A[%d ; %d] => ", i+1, j+1);
            scanf("%lf", &a[i][j]);
        }
    }

    // mostra a matriz
    printf("\nMatriz digitada:\n");
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < m; j++)
        {
            printf("%8.3f ", a[i][j]);
        }
        printf("\n");
    }

    // faz a transformação em um triangulo...
    for(i = 0; i < m - 1; i++)
    {
        if(a[i][i] == 0)
        {
            for(k = i; k < m; k++)
            {
                if(a[k][i] != 0)
                {
                    for(j = 0; j < m; j++)
                    {
                        temp = a[i][j];
                        a[i][j] = a[k][j];
                        a[k][j] = temp;
                    }
                    k = m;
                }
            }
            count++;
        }

        if(a[i][i] != 0)
        {
            for(k = i + 1; k < m; k++)
            {
                factor = -1.0 * a[k][i] /  a[i][i];
                for(j = i; j < m; j++)
                {
                    a[k][j] = a[k][j] + (factor * a[i][j]);
                }
            }
        }
    }

    temp = 1.0;
    // Calcula o determinante
    for(i = 0; i < m; i++)
        temp *= a[i][i];

    printf("\nDeterminante:\n");
    if(count % 2 == 0)
        printf("%8.3f \n", temp);
    else
        printf("%8.3f \n", -1.0 * temp);
}

void mostraMatriz(int n, float mat[n][n]){
    int i,j;
    for(i=0;i<n;i++){
        printf("|");
        for(j=0;j<n;j++){
            printf(" %5.2f",mat[i][j]);
        }
        printf(" |\n");
    }
}

void mostraMatrizT(int n, float mat[n][n]){
    int i,j;
    for(i=0;i<n;i++){
        printf("|");
        for(j=0;j<n;j++){
            printf(" %5.2f", mat[j][i]);
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
    printf("\n----y----\n");
    Y[0]=b[0]/M[0][0];
    for(i=1;i<n;i++){
        somat = 0;
        for(j=0;j<i;j++){
            somat += M[i][j]*Y[j];
        }
        Y[i] = (b[i] - somat )/M[i][i];
    }
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", Y[i]);
    }



    // SUBST REGRESSIVA
    printf("\n----x----\n");
    X[n-1] = Y[n-1]/aux[n-1][n-1];
    for(i=n-2;i>=0;i--){
        somat = 0;
        for(j=i+1;j<n;j++){
            somat += aux[i][j]*X[j];
        }
        X[i] = (Y[i] - somat )/aux[i][i];
    }
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", X[i]);
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
    printf("\n----y----\n");
    for(i = 0; i < n; ++i) {
        somat = 0;

        for(j = 0; j < i; ++j){
            somat += G[i][j] * Y[j];
        }
        Y[i] = (b[i] - somat) / G[i][i];
    }
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", Y[i]);
    }

    //Resolvendo Gt.x = y
    printf("\n----x----\n");
    for (i=n-1;i>=0;i--) {
        somat=0;
        for(j=i+1;j<n;j++){
            somat += G[j][i] * X[j];
        }
        X[i] = (Y[i] - somat) / G[i][i];
    }
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", X[i]);
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
    printf("\n---x---\n");
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

void gaussPivoteamentoParcial(int n, float A[n][n], float b[n]){

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

    printf("\n----b----\n");
    for(i=0;i<n;i++){
        printf("| %5.2f |\n", b[i]);
    }

    // SUBST REGRESSIVA
    printf("\n----x----\n");
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

int main(void){

    int n, i, j, metodo=99;

	printf("Programa desenvolvido para o Trabalho 2 da\ndisciplina de Metodos Numericos Computacionais\nda UNESP Bauru - 1 Sem 2020\n\n\n");

	// Definindo um tamanho fixo para a matriz no inicio do programa
	// Apenas para dispensar uso de alocacao dinamica e/ou um 'nMAX' pre-definido
	// Consumindo recursos sem necessidade
    printf("Digite a ordem N da matriz => ");
    scanf("%d", &n);

    float A[n][n], b[n];

    while (metodo!=0){
        printf("\n\n=============MENU===============");
        printf("\n========Decompos========================");
        printf("\n===  1-Determinante         ===");
        printf("\n===  2-Sistematriangulares(inf/sup) ===");
        printf("\n===  3-Gauss Siedel       ===");
        printf("\n===  4-Decomposicao LU       ===");
        printf("\n===  5-Cholesky              ===");
        printf("\n===  6-Gauss Simples         ===");
        printf("\n===  7-Gauss Pivo Parcial    ===");
        printf("\n===  8-Matriz inversa        ===");
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
			determinante();
            break;

        case 2:
			superiorInferior();
            break;

        case 3:
			gaussSiedel();
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
            printf("\n\n----b----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", b[i]);
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
            printf("\n\n----b----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", b[i]);
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
            printf("\n\n----b----\n");
            for(i=0;i<n;i++){
                printf("| %5.2f |\n", b[i]);
            }

            gaussSimples(n, A, b);

            break;
        
        case 8:
            matrizInversa();
            break;

        default:
            printf("\nOpcao Invalida!\nEscolha novamente...\n\n");
        }
    }
    
    return 0;
}
