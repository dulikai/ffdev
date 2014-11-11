
#ifndef _READFF_H_
#define _READFF_H_

#include <stdio.h>
#include <stdlib.h>
#include "utils/mytype.h"

#include "conndb.h"
/*
structure of ff database.
*/

#define FFtypeNum 100
#define FFbondNum 200
#define FFangleNum 300
#define FFdihNum 1000
#define FFimNum 50

#define MaxPara 10

#define INODE 0
#define JNODE 1
#define KNODE 2
#define LNODE 3

#define FF_bond_PARA_RE 0
#define FF_bond_PARA_KA 1

#define FF_angle_PARA_TH 0
#define FF_angle_PARA_KA 1

#define FF_dihedral_PARA_C0 0
#define FF_dihedral_PARA_C1 1
#define FF_dihedral_PARA_C2 2
#define FF_dihedral_PARA_C3 3
#define FF_dihedral_PARA_C4 4
#define FF_dihedral_PARA_C5 5

#define FF_improper_PARA_C0 0
#define FF_improper_PARA_C1 1
#define FF_improper_PARA_C2 2
#define FF_improper_PARA_C3 3
#define FF_improper_PARA_C4 4
#define FF_improper_PARA_C5 5


typedef struct{
        char ffname[10];
		int inode;
}ffndx;
        
typedef struct{
        int inode;
        char name[15];
        char element[5];
        int systype;
        }ffmap;     
        
typedef struct{
        int inode;
        decimal mass;
        decimal charge;
        decimal sig, eps;
        char field[20];
        }ffatomtype; 

typedef struct{
        int id;
        int idnst;
        int isource;
        }ffalias;     
        
typedef struct{
        int id;
        int inode, jnode;
        decimal re, ka;
        char type[10];
        }ffbond;            
        
typedef struct{
        int id;
        int inode, jnode, knode;
        decimal th, ka;
        char type[10];
        }ffangle;       


/*
typedef struct{
        int id;
		int inode, jnode, knode, lnode;
		decimal para[10];
		int npara;
		char type[10];
		int multi;
}ffdihedral; 
*/


       
typedef struct{
        int id;
        int inode, jnode, knode, lnode;
        decimal v1, v2, v3, v4;
        decimal para[10];
        int npara;
        
        char type[10];  
        }ffdihedral;          
      
typedef ffdihedral ffimproper;               



typedef struct{
        ffndx ndx[FFtypeNum];
        ffmap map[FFtypeNum];
        ffatomtype atomtype[FFtypeNum];
        ffalias alias[FFtypeNum];
        ffbond bond[FFbondNum];
        ffangle angle[FFangleNum];
        ffdihedral dihedral[FFdihNum];
        ffimproper improper[FFimNum];
    	char ffpname[100];
     	int nndx, natomtype, nmap, nalias, nbond, nangle, ndihedral, nimproper;
        }FFP;




/*
routines...
*/



  FFP *init_ffp();

  void read_ff_db(char *file, FFP *f);

  void print_ff_db(char *file, FFP *f);

#endif /* _READFF_H_ */
