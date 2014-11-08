
#ifndef VEC_H_
#define VEC_H_

#define XX	0			/* Defines for indexing in	*/
#define	YY	1			/* vectors			*/
#define ZZ	2
#define DIM   	3
#define invsqrt(x) (1.0f/sqrt(x))

typedef double       	dvec[DIM];


 void unitv(const dvec src,dvec dest);
 
 void oprod(const dvec a,const dvec b,dvec c);
 
 double iprod(const dvec a,const dvec b);
 
 double norm2(const dvec a);
 
 double norm(const dvec a);

 void gen_vec(double *p0, double *p1, dvec vec);

 void copy_vec(const dvec src, dvec dest);



#endif /* VEC_H_ */