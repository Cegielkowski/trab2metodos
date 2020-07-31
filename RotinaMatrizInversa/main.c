#include <stdio.h>
#include <math.h>
//----------------------------------------------------------
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
int main(void){
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