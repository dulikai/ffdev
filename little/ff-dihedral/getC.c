
/*
二面角参数的转换 
*/

#include<stdio.h>

int main(void)
{
    double v0[100], v1[100], v2[100], v3[100], v4[100],v5[100],v6[100];
    
    int i;
    
    for(i = 0; i < 100; i++)
    {
    v0[i] = v5[i] = v6[i] = 0;          
          }
    

    
    FILE *fp, *fpw;
    int nline;
    
    fp = fopen("abc.txt", "r");
    
    fscanf(fp, "%d", &nline);
    
    for(i = 0; i < nline; i++)
    {
          fscanf(fp, "%lf%lf%lf%lf", &v1[i], &v2[i], &v3[i], &v4[i]);
          printf("%lf %lf %lf %lf\n", v1[i], v2[i], v3[i], v4[i]);
          }
    fclose(fp);
    double c0[100], c1[100], c2[100], c3[100], c4[100], c5[100];
    fpw = fopen("d.txt", "w"); 
    for(i = 0; i < nline; i++)
    {
        c0[i] = v0[i] + v2[i] + 0.5*(v1[i]+v3[i]);
        c1[i] = 0.5*(3*v3[i] - v1[i]);
        c2[i] = -v2[i];
        c3[i] = -2*v3[i];
        c4[i] = 0;
        c5[i] = 0;
        }     
        
        for(i = 0; i < nline; i++)
        {
        fprintf(fpw, "%.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf\n", c0[i], c1[i], c2[i], c3[i], c4[i], c5[i]);
              }
              
              fclose(fpw);
              getch();
              return 0;    
    }


  void get_dihedral_para(double v0, double v1, double v2, double v3, double v4, double c[5])
  {
        c[0] = v0 + v2 + 0.5*(v1+v3);
        c[1] = 0.5*(3*v3 - v1);
        c[2] = -v2;
        c[3] = -2*v3;
        c[4] = 0;
        c[5] = 0;
        
        return;
       }
