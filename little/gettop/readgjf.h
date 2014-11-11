
#ifndef _READGJF_H_
#define _READGJF_H_


//#include "graphic.h"

#ifndef VertexNum
#define VertexNum 100
#endif
#ifndef INFINITY
#define INFINITY 10000
#endif


//typedef struct{}index

typedef struct{
        int id;
		char element_name[10];
        char ffname[10];
        char unique_name[15];
        int inode;        
        }vertex;
        
typedef vertex VerType;
typedef int EdgeType;

typedef struct{
        VerType vex[VertexNum];
        EdgeType edge[VertexNum][VertexNum];
		char name[50];
        int n, e;
        }structGraphic;
        
        

 void createGraphic(structGraphic *g, char *file);

/*
 fill graphic g, with original value.
*/
 void init_graphic(structGraphic *g);

 /*
  read gjf Cartesian file.
 */
 void read_gjf_cartesian(char *file, structGraphic *g);


 int read_gjf_match_ff(char *gjffile, char *suffix, structGraphic *g);

/*
  print the element in graphic g...
*/ 
 void print_graphic(char *file, structGraphic *g);
 
 

#endif  /* _READGJF_H_ */
