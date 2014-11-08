
#include<stdio.h>
#include<stdlib.h>

#include"createarr2.h"

/*
********************************START FUNCTION REALIZE********************************
*/




void** create_array2(int m,int n,int size)
{
int j,flag=1;
void** pp;

pp=(void**)malloc(m*sizeof(void*));

if(pp==NULL)
{
flag=0;
perror("\nmemory allocate error.\n");
}

for(j=0;j<m;j++)
{
pp[j]=(void*)malloc(n*size);
if(pp[j]==NULL){flag=0;}
}

return pp;

}



void free_array2(void** pp,int m)
{
	int j;
 for(j=0;j<m;j++)
 {
  if(pp[j]!=NULL)
  {
  free(pp[j]);
  }
 }

 free(pp);


}




 void ****create_array4(int m,int n,int p,int q,int size)
 {
 int i,j,k,flag=1;
 void**** p4;
  if ((p4 = (void ****)malloc(m * sizeof(void***))) == NULL ) {
    fprintf(stderr, "mb3dat: can't allocate memory");
  }

  for (i = 0; i < n; i++) 
  {
    if ((p4[i] = (void***)malloc(n*sizeof(void**))) == NULL)
	{
      fprintf(stderr, "mb3dat: can't allocate memory");
    }

    for (j = 0; j < p; j++)
	{
      if ((p4[i][j] = (void**)malloc(p*sizeof(void*))) == NULL) 
	  {
        fprintf(stderr, "mb3dat: can't allocate memory");
      }

     for(k=0;k<q;k++)
	 {
      if ((p4[i][j][k] = (void *)malloc(q*size)) == NULL) 
	  {
        fprintf(stderr, "mb3dat: can't allocate memory");
      }
	 }
	}
  }
return p4;
 }

void free_array4(double ****p4,int m,int n,int p,int q)
{
	int i,j,k;

	for(i = 0; i<m; i++)
		for(j=0;j<n;j++)
			for(k=0;k<p;k++)
			{
			  if(p4[j]!=NULL)
				{
					free(p4[i][j][k]);
				}
	
			}
}



double ***Array3D(int columns, int rows, int floors)
{
  double ***x;
  int i, j;
                                                                                                      
  if ((x = (double ***)malloc(columns * sizeof(double**))) == NULL ) {
    fprintf(stderr, "mb3dat: can't allocate memory");
  }
  for (i = 0; i < columns; i++) {
    if ((x[i] = (double **)malloc(rows * sizeof(double*))) == NULL) {
      fprintf(stderr, "mb3dat: can't allocate memory");
    }
    for (j = 0; j < rows; j++)
      if ((x[i][j] = (double *)malloc(floors * sizeof(double))) == NULL) {
        fprintf(stderr, "mb3dat: can't allocate memory");
      }
  }
  return x;
}



void FreeArray3D(double ***uu,int columns, int rows, int floors)
{
   int i, j;
                                                                                                      
   for (i = 0; i < columns; i++)
     {
        for (j = 0; j < rows; j++){
           free(uu[i][j]);
        }
        free(uu[i]);
     }
   free(uu);
   return;
} // end FreeArray3D()



/*
main()
{
int i,j,k,l,a=1;
int **pp,****p4;



pp=(int**)create_array2(3,4,sizeof(int));
 p4=(int****)create_array4(3,3,3,3,sizeof(int));

for(i=0;i<3;i++)
for(j=0;j<3;j++)
{
for(k=0;k<3;k++)
for(l=0;l<3;l++)
{
		p4[i][j][k][l]=a++;
		printf("%4d,",p4[i][j][k][l]);

}
putchar('\n');
}
printf("\n");

for(i=0;i<3;i++)
{
	for(j=0;j<4;j++)
	{
		pp[i][j]=i;
		printf("%d,",pp[i][j]);
	
	}
}
free_array2(pp,3);

}

*/


/*
************************************************************************
*/

