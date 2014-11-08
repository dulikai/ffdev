
#include <math.h>
#include "vec.h"


 void unitv(const dvec src,dvec dest)
{
  double linv;
  
  linv=invsqrt(norm2(src));
  dest[XX]=linv*src[XX];
  dest[YY]=linv*src[YY];
  dest[ZZ]=linv*src[ZZ];
}


 void oprod(const dvec a,const dvec b,dvec c)
{
  c[XX]=a[YY]*b[ZZ]-a[ZZ]*b[YY];
  c[YY]=a[ZZ]*b[XX]-a[XX]*b[ZZ];
  c[ZZ]=a[XX]*b[YY]-a[YY]*b[XX];
}

double iprod(const dvec a,const dvec b)
{
  return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}



double norm2(const dvec a)
{
  return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

double norm(const dvec a)
{
  return sqrt(a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]);
}

void gen_vec(double *p0, double *p1, dvec vec)
{
  vec[XX] = p1[XX] - p0[XX];
  vec[YY] = p1[YY] - p0[YY];
  vec[ZZ] = p1[ZZ] - p0[ZZ];

  return;
}

 void copy_vec(const dvec src, dvec dest)
 {
  dest[XX]=src[XX];
  dest[YY]=src[YY];
  dest[ZZ]=src[ZZ];

  return;
 }

