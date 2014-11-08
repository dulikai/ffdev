
#include<math.h>
#include<stdio.h>
#include <stdlib.h>
#include "matrixop.h"






/***********************************************************************************
+++R+++++++++++++++++++++++++++++++++| for ERI |+++++++++++++++++++++++++++++++++++++++
---------------------------------------------------------------------------------*/



  int num_sht(int i,int j)
  {
   int m,it,jt;

   it=i;jt=j;
   if(i>j)
   {
     it=j;jt=i;
   }
    m=jt*(jt+1)/2+it;

   return m;
  }

  int num4_sht(int i,int j,int m,int n)
  {
   int t,p,q,ix;

   if(i>j)
   {
    t=i;i=j;j=t;
   }
   if(m>n)
   {
    t=m;m=n;n=t;
   }

   p=j*(j+1)/2+i;
   q=n*(n+1)/2+m;

   if(p>q)
   {
    t=p;p=q;q=t;
   }

   ix=q*(q+1)/2+p;

   return ix;
  }












/*
==============================================================================================
********************************** OPERATIONS ON MATRIX **************************************
///////////////////////////////////////////////////////////////////////////////////////////////
----------------------------------------------------------------------------------------------
*/
/*--------------------------------------------------------------------------------------------
++++R+++++++++++++++++++++++++++++ADDITION+++SUBTRATION+++++++++++++++++++++++++++++++++++++++++
----------------------------------------------------------------------------------------------*/

 void matrix_zero(double **matrix1,int n)
 {
  int i,j;

  for(i=0;i<n;i++)
	  for(j=0;j<n;j++)
	  {
	    matrix1[i][j]=0.0;
	  }
 }


 void matrix_ab(double **matrix_a,double **matrix_b,double **matrix_c,int n,int sign)
 {
  int i,j;
  double t;

  for(i=0;i<n;i++)
  {
   for(j=0;j<n;j++)
   {
	t=matrix_b[i][j];

	if(sign<0)
	{
	 t=-t;
	}

    matrix_c[i][j]=matrix_a[i][j]+t;
   }
  }
 }

 void matrix_self_ab(double **matrix_a,double **matrix_b,int n,int sign)
 {
  int i,j;
  double t;

  for(i=0;i<n;i++)
  {
   for(j=0;j<n;j++)
   {
	t=matrix_b[i][j];

	if(sign<0)
	{
	 t=-t;
	}

    matrix_a[i][j] = matrix_a[i][j]+t;
   }
  }
   
 }





/*--------------------------------------------------------------------------------------------
++++R++++++++++++++++++++++++++++++multiplicationI+++++++++++++++++++++++++++++++++++++++++
----------------------------------------------------------------------------------------------*/
/*
 square matrix only.
 */
 void matrix_multi(double **matrix_a,double **matrix_b,double **matrix_c,int n)
 {
  int i,j,k;

  for(i=0;i<n;i++)
  
	  for(j=0;j<n;j++)
	  {
	   matrix_c[i][j]=0.0;

	   for(k=0;k<n;k++)
	   {
		   matrix_c[i][j]+=matrix_a[i][k]*matrix_b[k][j];
		   /* Cij = Aik * Bkj */
	   }
	  }
 }

 void matrix_multi_const(double **a,double **b,int n,double lambda)
 {
  int i,j;

  for(i=0;i<n;i++)
	  for(j=0;j<n;j++)
	  {
	   a[i][j]=lambda*b[i][j];
	  }
 }


/*--------------------------------------------------------------------------------------------
+++++R++++++++++++++++++++++++++multiplicationI+++++++++++++++++++++++++++++++++++++++++
----------------------------------------------------------------------------------------------*/


 void matrix_turn(double **matrix_a,int n)
 {
  int i,j;
  double t;

  for(i=0;i<n;i++)
	  for(j=i;j<n;j++)
	  {
		  t=matrix_a[i][j];
		  matrix_a[i][j]=matrix_a[j][i];
		  matrix_a[j][i]=t;
	  }
 }


 void matrix_turn_copy(double **matrix_a,double **matrix_b,int n)
 {
  int i,j;

  for(i=0;i<n;i++)
	  for(j=0;j<n;j++)
	  {
		  matrix_b[j][i]=matrix_a[i][j];
/*		  matrix_a[i][j]=matrix_a[j][i];
		  matrix_a[j][i]=t;
*/
	  }
 }


 double matrix_elem_multi_add(double **matrix_a,double **matrix_b,int n)
 {
	int i,j;
	double sum=0.0;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			sum+=matrix_a[i][j] * matrix_b[i][j];
		}
	return sum;
 }


 double matrix_trace(double **matrix_a,int n)
 {
	int i;
	double sum=0.0;

	for(i=0;i<n;i++)
	{
		sum+=matrix_a[i][i];
	}
	return sum;
 }

 /*
  in:double **mat_a, double **mat_b, int n
  out: return tr(mat_a * mat_b);
 */

 double matrix_m2_trace(double **mat_a, double **mat_b, int n)
 {
  int i,j;
  double emt;

  emt = 0.0;

  for(i = 0; i < n; i++)
	  for(j = 0;j < n; j++)
	  {
        emt += mat_a[i][j] * mat_b[j][i];
	  }

	  return emt;
 }

 void matrix_copy(double **matrix_a,double **matrix_b,int n)
 {
	int i,j;

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			matrix_b[i][j]=matrix_a[i][j];
		}
		return ;
 }


/*
  in:double *vec_ina, double *vec_inb, int n
  out: return dot value;
  used for dot multi of vector,double *vec_ina, double *vec_inb  can be different of the same;
 */
 double vector_dot(double *vec_ina, double *vec_inb, int n)
 {
   int i;
   double sum;

   sum = 0.0;

   for(i = 0; i < n; i++)
	   sum += vec_ina[i] * vec_inb[i];

   return sum;
 }



/*
 in:double *vec_in, double **mat_in,int n,int vec_left_flag
 out:double *vec_out,
if vec_left_flag == 0,mean mat_in * vec_in
else mean vec_in * mat_in

  warning : mat must be square matrix.
 */
 void vec_mat_multi(double *vec_in, double **mat_in,  double *vec_out, int n,int vec_left_flag)
 {
  int i,j;

/*
  get vec * mat
  */
  if(vec_left_flag)
  {
	for(i = 0; i < n; i++)
	{
		vec_out[i] = 0.0;
		for(j = 0; j< n; j++)
		{
			 vec_out[i] += vec_in[j] * mat_in[j][i];
		}
	}
  }
  else
  {
	for(i = 0; i < n; i++)
	{
		vec_out[i] = 0.0;
		for(j = 0; j< n; j++)
		{
			vec_out[i] += mat_in[i][j] * vec_in[j] ;
		}
	}

  }
 }

/*
  in: double *vec_ina, double *vec_inb, double *vec_out, int n, int flag
  out:double *vec_out,

cite:
  if(flag==0)the vec_out = vec_ina - vec_inb;
  else vec_out = vec_ina + vec_inb;
 */
 void vec_add_sub(double *vec_ina, double *vec_inb, double *vec_out, int n, int flag)
 {
  int i;

  if(flag)
  {
     for(i = 0;i < n; i++)
	     vec_out[i] = vec_ina[i] + vec_inb[i];
  }

  else
  {
    for(i = 0; i < n; i++)
		vec_out[i] = vec_ina[i] - vec_inb[i];
  }

  return;
 }





/*
 this is copyed from a lib,which calculate matrix int the form of 1D array,
 i will changed it int 2D array,int the following function.

  Gauss-Jordan method 

  in: double a[],int n
  out: double a[];
  purpose:get the inverse of a matrix.
*/

  int brinv(double a[],int n)
  { int *is,*js,i,j,k,l,u,v;
    double d,p;
    is=(int *)malloc(n*sizeof(int));
    js=(int *)malloc(n*sizeof(int));
    for (k=0; k<=n-1; k++)
      { d=0.0;
        for (i=k; i<=n-1; i++)
        for (j=k; j<=n-1; j++)
          { l=i*n+j; p=fabs(a[l]);
            if (p>d) { d=p; is[k]=i; js[k]=j;}
          }
        if (d+1.0==1.0)
          { free(is); free(js); printf("err**not inv\n");
            return(0);
          }
        if (is[k]!=k)
          for (j=0; j<=n-1; j++)
            { u=k*n+j; v=is[k]*n+j;
              p=a[u]; a[u]=a[v]; a[v]=p;
            }
        if (js[k]!=k)
          for (i=0; i<=n-1; i++)
            { u=i*n+k; v=i*n+js[k];
              p=a[u]; a[u]=a[v]; a[v]=p;
            }
        l=k*n+k;
        a[l]=1.0/a[l];
        for (j=0; j<=n-1; j++)
          if (j!=k)
            { u=k*n+j; a[u]=a[u]*a[l];}
        for (i=0; i<=n-1; i++)
          if (i!=k)
            for (j=0; j<=n-1; j++)
              if (j!=k)
                { u=i*n+j;
                  a[u]=a[u]-a[i*n+k]*a[k*n+j];
                }
        for (i=0; i<=n-1; i++)
          if (i!=k)
            { u=i*n+k; a[u]=-a[u]*a[l];}
      }
    for (k=n-1; k>=0; k--)
      { if (js[k]!=k)
          for (j=0; j<=n-1; j++)
            { u=k*n+j; v=js[k]*n+j;
              p=a[u]; a[u]=a[v]; a[v]=p;
            }
        if (is[k]!=k)
          for (i=0; i<=n-1; i++)
            { u=i*n+k; v=i*n+is[k];
              p=a[u]; a[u]=a[v]; a[v]=p;
            }
      }
    free(is); free(js);
    return(1);
  }

/*
  converter of the functon 'brinv'
  */

int matrix_inv(double **mat,int n)
{
	int i,j,k;
	double *a;

	a=(double *)malloc(n*n* sizeof(double));

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			k=i*n+j;
			a[k] = mat[i][j];
		}

	brinv(a,n);

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			k=i*n+j;
			mat[i][j] = a[k] ;
		}

	free(a);

	return 1;
}




/*
  converter of the functon 'brinv'
  */

int matrix_copy_inv(double **mat,double **invmat,int n)
{
	int i,j,k;
	double *a;

	a=(double *)malloc(n*n* sizeof(double));

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			k=i*n+j;
			a[k] = mat[i][j];
		}

	brinv(a,n);

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			k=i*n+j;
			invmat[i][j] = a[k] ;
		}

	free(a);

	return 1;
}






/*
Gauss-Jordan complete pivoting Elimination 

  AX=B

in:  double a[], matrix A;
	 double b[],matrix B,right vector;
	 n: rank of the linear function;
	 m: number of const right vector

out:double a[], matrix A: destroyed.
    double b[],matrix B : root of the linear function;

*/

  int agjdn(double a[],double b[],int n,int m)
  { int *js,l,k,i,j,is,p,q;
    double d,t;
    js=(int *)malloc(n*sizeof(int));
    l=1;
    for (k=0;k<=n-1;k++)
      { d=0.0;
        for (i=k;i<=n-1;i++)
          for (j=k;j<=n-1;j++)
            { t=fabs(a[i*n+j]);
              if (t>d) { d=t; js[k]=j; is=i;}
            }
        if (d+1.0==1.0) l=0;
        else
          { if (js[k]!=k)
              for (i=0;i<=n-1;i++)
                { p=i*n+k; q=i*n+js[k];
                  t=a[p]; a[p]=a[q]; a[q]=t;
                }
            if (is!=k)
              { for (j=k;j<=n-1;j++)
                  { p=k*n+j; q=is*n+j;
                    t=a[p]; a[p]=a[q]; a[q]=t;
                  }
                for (j=0;j<=m-1;j++)
                  { p=k*m+j; q=is*m+j;
                    t=b[p]; b[p]=b[q]; b[q]=t;
                  }
              }
          }
        if (l==0)
          { free(js); printf("fail\n");
            return(0);
          }
        d=a[k*n+k];
        for (j=k+1;j<=n-1;j++)
          { p=k*n+j; a[p]=a[p]/d;}
        for (j=0;j<=m-1;j++)
          { p=k*m+j; b[p]=b[p]/d;}
        for (j=k+1;j<=n-1;j++)
          for (i=0;i<=n-1;i++)
            { p=i*n+j;
              if (i!=k)
                a[p]=a[p]-a[i*n+k]*a[k*n+j];
            }
        for (j=0;j<=m-1;j++)
        for (i=0;i<=n-1;i++)
          { p=i*m+j;
            if (i!=k)
              b[p]=b[p]-a[i*n+k]*b[k*m+j];
          }
      }
    for (k=n-1;k>=0;k--)
      if (js[k]!=k)
        for (j=0;j<=m-1;j++)
          { p=k*m+j; q=js[k]*m+j;
            t=b[p]; b[p]=b[q]; b[q]=t;
          }
    free(js);
    return(1);
  }

 /*
 Gauss-Jordan complete pivoting Elimination 

  AX=B

  form convert of the function 'agjdn'

  in: double **mat_a, : matrix A;
	  double **mat_b, : right vector;
	  int n, : rank of the linear function;
	  int m, : number of const right vector;
  out:double **mat_a,: destroyed.
	  double **mat_b, : root of the linear function;
 */
int matrix_gaussj(double **mat_a,double **mat_b,int n,int m)
{
	int i,j,k;
	double *a,*b;

	a=(double *)malloc(n*n*sizeof(double));
	b=(double *)malloc(n*m*sizeof(double));

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			k=i*n+j;
			a[k] = mat_a[i][j];
		}


		if(m == 1)
		{
	      for(i=0;i<n;i++)
		    for(j=0;j<m;j++)
			{
		 	 k=i;
			 b[k] = mat_b[j][i];
			}
		}
		else
		{
	     for(i=0;i<n;i++)
		   for(j=0;j<m;j++)
		   {
			 k=i*m+j;
			 b[k] = mat_b[i][j];
		   }
		}


    agjdn(a,b,n,m);


	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
		{
			k=i*n+j;
			mat_a[i][j] = a[k];
		}



		if(m == 1)
		{
	      for(i=0;i<n;i++)
		    for(j=0;j<m;j++)
			{
		 	 k=i;
			 mat_b[j][i] = b[k];
			}
		}
		else
		{
	      for(i=0;i<n;i++)
		    for(j=0;j<m;j++)
			{
		    	k=i*m+j;
		        mat_b[i][j] = b[k];
			}
		}

	free(a);
	free(b);

	return 1;
	
}




 double matrix_abs_max(double **mat_a,int n)
 {
  int i,j;
  double big;

  big = fabs(mat_a[0][0]);
  
  for(i=0;i<n;i++)
	  for(j=0;j<n;j++)
		  if(big < fabs(mat_a[i][j]))
		  {
	        big = fabs(mat_a[i][j]);
		  }

	return big;
 }





 /*
-------------------------------------------------------------------------------
 v001,22:56 2006-12-31;

 DESCRIPTION: GET vector dot multi. VALUE

 pre: double *vect1,*vect2:(INPUT);
	  int n : (input) 
 post: RETURN vector dot multi value.

 JCP 84(7),1986
---------------------------------------------------------------------------------
 */
  double vector_dot_multi(double const *vect1,double const *vect2,int const n)
  {
   int i;
   double sum=0.0;

   for(i=0;i<n;i++)
	   {
	    sum += vect1[i]*vect2[i];
	   }
   return sum;
  }














