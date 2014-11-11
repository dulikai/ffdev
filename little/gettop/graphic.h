
/*
  Name: 
  Copyright: 
  Author: 
  Date: 16-12-08 09:48
  Description: 
*/

#ifndef _GRAPHIC_H_
#define _GRAPHIC_H_


#include "readgjf.h"

#ifndef VertexNum
#define VertexNum 100
#endif
#ifndef INFINITY
#define INFINITY 10000
#endif


typedef vertex com_vertex;

typedef struct{
        int i;
        int j;
        }com_bond;
        
typedef struct{
        int i;
        int j;
        int k;
        }com_angle;
        
typedef struct{
        int i;
        int j;
        int k;
        int l;
        }com_dihedral;
        
typedef struct{
        int i;
        int l;
        }com_pair;
                                
typedef struct{
        com_vertex vex[VertexNum];
        com_bond bond[VertexNum];
        com_angle angle[VertexNum];
        com_dihedral dihedral[VertexNum];
        com_pair pair[VertexNum];
        
        int nvex, nbond, nangle, ndihedral, npair;
        
        }analyzeGraphic;
        
        
        



 static void gen_vex(structGraphic *g, analyzeGraphic *com);
 static void gen_bond(structGraphic *g, analyzeGraphic *com);    
 static void gen_angle(structGraphic *g, analyzeGraphic *com);
 static void gen_dihedral(structGraphic *g, analyzeGraphic *com);    
 static void gen_pair(structGraphic *g, analyzeGraphic *com);
 static int is_redundance_pair(analyzeGraphic *com, int npair, int u, int v);
 static void Dijkstra_distance(structGraphic *g, int v0, int *dist);
 static void assign_vex_name(structGraphic *g);


 void gen_com(structGraphic *g, analyzeGraphic *com);
 void print_com(char *file, analyzeGraphic *com);
 int is_self_iter(int list[], int n);

 int get_inode_form_id(analyzeGraphic *com, int id);
 void get_element_from_com(analyzeGraphic *com, int inode, char *name);


#endif /* _GRAPHIC_H_ */
