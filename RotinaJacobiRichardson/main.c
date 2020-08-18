//Método Jacobi Richardson
#include <stdio.h>

float coeficientes[10][10], dinversa[10][10], aproximacao[10][1], independentes[10][1], erro, matrizR[10][10], temporaria[10][1], matrizResult[10][1];
int ordem, linha, coluna, iteracoes;

void multiplicacao(float matrizA[][10], float matrizB[][1])
{
	int i, j;
	for (i = 0; i < ordem; i++)
	{
		matrizResult[i][0]=0;
		for (j = 0; j < ordem; j++)
		matrizResult[i][0]=matrizResult[i][0]+matrizA[i][j]*matrizB[j][0];
	}
}

int main()
{
	printf("Digite a ordem do sistema: \n");
	scanf("%d",&ordem);
	
	printf("Digite a matriz dos coeficientes: \n");
	for (linha = 0; linha < ordem; linha++)
		for (coluna = 0; coluna < ordem; coluna++)
			scanf("%f",&coeficientes[linha][coluna]);
			
	printf("Digite o vetor de termos independentes: \n");
	for (linha = 0; linha < ordem; linha++)
		scanf("%f", &independentes[linha][0]);
			
	printf("Digite a aproximacao inicial: \n");
	for (linha = 0; linha < ordem; linha++)
		scanf("%f", &aproximacao[linha][0]);

	printf("Digite o número maximo de iteracoes: \n");
	scanf("%d", &iteracoes);
	
	printf("Digite a precisao (e): \n");
	scanf("%f", &erro);
	
	for (linha = 0; linha < ordem; linha++)
		for (coluna = 0; coluna < ordem; coluna++)
		{
			if (linha == coluna)
				dinversa[linha][coluna] = 1/coeficientes[linha][coluna];
			else
				dinversa[linha][coluna] = 0;
		}
		for (linha = 0; linha < ordem; linha ++)
			for (coluna = 0; coluna < ordem; coluna++)
			{
				if (linha == coluna)
					matrizR[linha][coluna] = 0;
				else
				if (linha != coluna)
					matrizR[linha][coluna] = coeficientes[linha][coluna];
			}
	
	int i = 1;
	int j;
	while (i <= iteracoes)
	{
		multiplicacao(matrizR,aproximacao);
		for (linha = 0; linha < ordem; linha ++)
			temporaria[linha][0] = independentes[linha][0] - matrizResult[linha][0];
		multiplicacao(dinversa, temporaria);
		for (j = 0; j < ordem; j++)
			aproximacao[j][0] = matrizResult[j][0];
		printf("Os valores de x apos a interacao %d sao\n",i);
		for (linha = 0; linha < ordem; linha++)
			printf("%.3f\n", aproximacao[linha][0]);
		i++;	
	}
	printf("\nNumero maximo de iteracoes atingido!");
	getch();
}
