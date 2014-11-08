

#ifndef _MATRIXOP_H_
#define _MATRIXOP_H_




/*
 2D to 1D
*/
 int num_sht(int i,int j);

/*
 4D to 1D pos.
 */
 int num4_sht(int i,int j,int m,int n);





/*
 in: int n
 out:double **matrix1

 cite:give zero value to the whole matrix;
*/
 void matrix_zero(double **matrix1,int n);

/*
 in:double **matrix_a,double **matrix_b,int n,int sign
 out:double **matrix_c,

 cite:  if variable sign < 0  then - ; else +;
 */
 void matrix_ab(double **matrix_a,double **matrix_b,double **matrix_c,int n,int sign);
 

/*
 in:double **matrix_a,double **matrix_b,int n,int sign
 out:double **matrix_a,

 cite:  if variable sign < 0  then - ; else +;
		in this function matrix_a is added or subed to matrix_a,not another matrix,return matrix_a.
 */
 void matrix_self_ab(double **matrix_a,double **matrix_b,int n,int sign);








 /*
  in: double **matrix_a,double **matrix_b,int n;
  out:double **matrix_c;
  cite: multi. two matrix
 */
 void matrix_multi(double **matrix_a,double **matrix_b,double **matrix_c,int n);


/*
 name: matrix_multi_const
 in:double **b,int n,double lambda
 out:double **a,
 cite: aij = lambda * bij ;
  */
 void matrix_multi_const(double **a,double **b,int n,double lambda);



 /*
  in:double **matrix_a,int n
  out:double **matrix_a,
  cite: turn of a matrix;
 */
  void matrix_turn(double **matrix_a,int n);

 /*
  in:double **matrix_a,int n
  out:double **matrix_b,
  cite: turn of a matrix ,but matrix_a is not changed;
 */

 void matrix_turn_copy(double **matrix_a,double **matrix_b,int n);

/*
 in:double **matrix_a,double **matrix_b,int n
 out:return added value.
 */
 double matrix_elem_multi_add(double **matrix_a,double **matrix_b,int n);

/*
 in:double **matrix_a,int n
 out: trace of the matrix
 */
 double matrix_trace(double **matrix_a,int n);


 /*
  in:double **mat_a, double **mat_b, int n
  out: return tr(mat_a * mat_b);
 */

 double matrix_m2_trace(double **mat_a, double **mat_b, int n);


/*
 in:double **matrix_ina,int n
 out:double **matrix_b,
 cite: copy matrix from a to b
 */
 void matrix_copy(double **matrix_ina,double **matrix_b,int n);



/*
  in:double *vec_ina, double *vec_inb, int n
  out: return dot value;
  used for dot multi of vector,double *vec_ina, double *vec_inb  can be different of the same;
 */
 double vector_dot(double *vec_ina, double *vec_inb, int n);






/*
 in:double *vec_in, double **mat_in,int n,int vec_left_flag
 out:double *vec_out,
if vec_left_flag == 0,mean mat_in * vec_in
else mean vec_in * mat_in

 */
 void vec_mat_multi(double *vec_in, double **mat_in,  double *vec_out, int n,int vec_left_flag);



/*
  in: double *vec_ina, double *vec_inb, double *vec_out, int n, int flag
  out:double *vec_out,

cite:
  if(flag==0)the vec_out = vec_ina - vec_inb;
  else vec_out = vec_ina + vec_inb;
 */
 void vec_add_sub(double *vec_ina, double *vec_inb, double *vec_out, int n, int flag);









/*
 this is copyed from a lib,which calculate matrix int the form of 1D array,
 i will changed it int 2D array,int the following function.

  Gauss-Jordan method 

  in: double a[],int n
  out: double a[];
  purpose:get the inverse of a matrix.
*/
 static int brinv(double a[],int n);

 int matrix_inv(double **mat,int n);
/* inverse of the mat is invmat.  */
 int matrix_copy_inv(double **mat,double **invmat,int n);






/*
Gauss-Jordan complete pivoting Elimination 

  AX=B

in:  double a[], matrix A;
	 double b[],matrix B,right vector;
	 n: rank of the linear function;
	 m: number of const right vector;

out:double a[], matrix A: destroyed.
    double b[],matrix B : root of the linear function;

*/
 static int agjdn(double a[],double b[],int n,int m);

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
 int matrix_gaussj(double **mat_a,double **mat_b,int n,int m);


/*
 in: double **mat_a,int n
 out: return the largest element of the matrix input.

cited: abs(element)
 
 */
 double matrix_abs_max(double **mat_a,int n);


/*
------------------------------------------------------------------------------
  DESCRIPTION: GET vector dot multi. VALUE
---------------------------------------------------------------------------------

 */
 double vector_dot_multi(double const *vect1,double const *vect2,int const n);



#endif  /* _MATRIXOP_H_ */