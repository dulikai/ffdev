
#ifndef _CREATE_ARRA2_H
#define _CREATE_ARRA2_H

void** create_array2(int m,int n,int size);
void free_array2(void** pp,int m);
void ****create_array4(int m,int n,int p,int q,int size);
void free_array4(double ****p4,int m,int n,int p,int q);
double ***Array3D(int columns, int rows, int floors);
void FreeArray3D(double ***uu,int columns, int rows, int floors);

#endif  /* _CREATE_ARRA2_H */ 