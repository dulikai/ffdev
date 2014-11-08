
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include "vec.h"
#include "tools/matrixop.h"
#include "tools/rdiotools.h"
#include "tools/createarr2.h"

#define NMAXPOINT 1000

typedef struct{
	double origin[3];
	double ex[3];
	double ey[3];
	double ez[3];

	double **rot;
}COORdinate;


 typedef struct{
         double value[3];
		 int id;
 }COORD;

 typedef struct{
	 COORD src[NMAXPOINT];
	 COORD dest[NMAXPOINT];
	 int nsrc, ndest;

	 COORdinate *coord;

	 COORD refsrc[3];
	 COORD refdest[3];
 
 }SysCoord;




 void gen_rot_mat_3d(double **rot, double *ex, double *ey, double *ez)
 {
   int i;

   for(i = 0; i < 3; i++)
   {
     rot[i][0] = ex[i];
	 rot[i][1] = ey[i];
	 rot[i][2] = ez[i];
   }

   return;
 }


  void mat_vec_mul_3d(double **mat, double *vec, int m, int n, double *dest)
  {
    int i, j;
	double s;

	for(i = 0; i < m; i++)
	{
      s = 0.0;
	  for(j = 0; j < n; j++)
	  {
	    s += mat[i][j] * vec[j];
	  }
	  dest[i] = s;
	}

	return;
  }

  void sub_vec(double *a, double *b, double *c, int n)
  {
    int i;

	for(i = 0; i < n; i++)
	{
	  c[i] = a[i] - b[i];
	}

	return;
  }

  void gen_new_point(double **rot, double *trans, double *old_point, double *new_point)
  {
    double dest[3];

    mat_vec_mul_3d(rot, old_point, 3, 3, dest);

	sub_vec(dest, trans, new_point, 3);

	return;
    
  }

  void gen_xyz_orign(SysCoord *sys, double *vec01, double *vec02)
  {
    double h;
    double vec0102c[3];
	double new_p2[3];
	double **rot, *p2;

	rot = sys->coord->rot;
    p2 = sys->refsrc[2].value;

    oprod(vec01, vec02, vec0102c);

	h = norm(vec0102c)/norm(vec01); 
	printf("\nheight of the triangle is %.6lf\n\n", h);
	new_p2[0] = 0.0;
    new_p2[1] = h;
    new_p2[2] = 0.0;

    gen_new_point(rot, new_p2, p2, sys->coord->origin);

	return;  
  }

  void gen_xyz_unit_vec(SysCoord *sys, double *vec01, double *vec02)
  {
    double ex[3], ey[3], ez[3];
    double vec0102c[3], vec31c[3];
    double **rot;

    rot = (double **)create_array2(3, 3, sizeof(double));

    //x
    unitv(vec01, ex);

    //z  vec0102c = vec01*vec02, cross multiple...
    oprod(vec01, vec02, vec0102c);
    unitv(vec0102c, ez);

    //y ,vec31c = e3*e1, cross multiple...
    oprod(ez, ex, vec31c);
    unitv(vec31c, ey);

    copy_vec(ex, sys->coord->ex);
	copy_vec(ey, sys->coord->ey);
	copy_vec(ez, sys->coord->ez);

	//use new xyz unit vector to generate rotation matrix....
    gen_rot_mat_3d(rot, ex, ey, ez);
	matrix_inv(rot, 3);
	sys->coord->rot = rot;

    
	return;
  }




  void gen_coordinate(SysCoord *sys)
  {
    COORdinate *coord;
	double *p0, *p1, *p2;
	double vec01[3], vec02[3];

	p0 = sys->refsrc[0].value;
	p1 = sys->refsrc[1].value;
	p2 = sys->refsrc[2].value;

	coord = sys->coord;

	//generate vector....
    gen_vec(p0, p1, vec01);
    gen_vec(p0, p2, vec02);

	//new xyz unit vector...
    gen_xyz_unit_vec(sys, vec01, vec02);

    //get the position of the new origin....
    gen_xyz_orign(sys, vec01, vec02);
	
	return;

  }



 void read_sys_coord(char *file, SysCoord *sys)
 {
   FILE *fp;
   char line[200];
   int i;
   int icoord = 0;

   fp = fopen(file, "r");
   if(fp == NULL)
   {
     printf("file:%s does't exit...\n");
	 system("PAUSE");
	 exit(0);
   }

   JmpOrNot(fp, '#'); 

   for(i = 0; i < 3; i++)
   {
     fscanf(fp, "%lf%lf%lf",
		 &sys->refsrc[i].value[0], &sys->refsrc[i].value[1], &sys->refsrc[i].value[2]); 
   }

   goto_target_line(fp, '#', 0L);   
   JmpOrNot(fp, '#'); 

   while(fgets(line, sizeof(line), fp) != NULL)
   {
      printf("%s", line);
      if(is_blank_line(line))
	  {
        break;
	  }
      else
      {
         sscanf(line, "%lf%lf%lf", 
			 &sys->src[icoord].value[0], &sys->src[icoord].value[1], &sys->src[icoord].value[2]);
         icoord++;
		 sys->src[icoord].id = icoord;
	  }   
	  
   }

   sys->coord = (COORdinate*)malloc(sizeof(COORdinate));

   sys->nsrc = icoord;

   return;
 }


 void cacl_new_points(SysCoord *sys)
 {
   int i;
   int nsrc;

   nsrc = sys->nsrc;

   for(i = 0; i < 3; i++)
   {
     gen_new_point(sys->coord->rot, sys->coord->origin, sys->refsrc[i].value, sys->refdest[i].value);
   }

   for(i = 0; i < nsrc; i++)
   {
     gen_new_point(sys->coord->rot, sys->coord->origin, sys->src[i].value, sys->dest[i].value);

	 sys->dest[i].id = sys->src[i].id;
   }
   sys->ndest = sys->nsrc;

   return;
 }


 void print_sys_point(char *file, SysCoord *sys)
 {
   FILE *fp;
   int i, j;
   int ndest;

   ndest = sys->ndest;

   fp = fopen(file, "w");

   fprintf(fp, "\n新坐标系的原点在旧坐标系的坐标为：\n");
   for(i = 0; i < 3; i++)
	   fprintf(fp, "%10.6lf ", sys->coord->origin[i]);


   fprintf(fp, "\n\n旋转矩阵为\n\n");
   for(i = 0; i < 3; i++)
   {
	   for(j = 0; j < 3; j++)
	   {
		   fprintf(fp, "%10.6lf ", sys->coord->rot[i][j]);
	   }
	   fprintf(fp, "\n");
   }

   fprintf(fp, "\n\nthe three point you enter, now converted in the now coordinate..\n\n");

   for(i = 0; i < 3; i++)
   {
      fprintf(fp, "%10.6lf %10.6lf %10.6lf  ==>  ", 
		  sys->refsrc[i].value[0], sys->refsrc[i].value[1], sys->refsrc[i].value[2]);

      fprintf(fp, "%10.6lf %10.6lf %10.6lf\n", 
		  sys->refdest[i].value[0], sys->refdest[i].value[1], sys->refdest[i].value[2]);
   }

   fprintf(fp, "\n\n所有需要转换的坐标，在这里输出\n\n");
   for(i = 0; i < ndest; i++)
   {
      fprintf(fp, "%10.6lf %10.6lf %10.6lf  ==>  ", 
		  sys->src[i].value[0], sys->src[i].value[1], sys->src[i].value[2]);

      fprintf(fp, "%10.6lf %10.6lf %10.6lf\n", 
		  sys->dest[i].value[0], sys->dest[i].value[1], sys->dest[i].value[2]);
   }

   return;
 }


 int main(void)
 {
   char *input = "input.txt";
   char *output = "output.txt";
   SysCoord sys;

   read_sys_coord(input, &sys);
   gen_coordinate(&sys);
   cacl_new_points(&sys);
   print_sys_point(output, &sys);

   system("PAUSE");

   return 0;

 }