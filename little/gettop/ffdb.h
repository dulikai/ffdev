
#ifndef _FFDB_H_
#define _FFDB_H_

#include "readff.h"

typedef ffndx affndx;
typedef ffatomtype affatomtype;
typedef ffalias afftrans;


typedef struct{
        int id;
        int node[2];
        decimal para[MaxPara];
		int npara;
        char type[10];
}affbond;

typedef struct{
        int id;
        int node[3];
        decimal para[MaxPara];
		int npara;
        char type[10];
}affangle;

typedef struct{
        int id;
        int node[4];
        decimal para[MaxPara];
        int npara;        
        char type[10];  
        int multi;      
}affdihedral;

typedef affdihedral affimproper;
/*
typedef struct{
	int funct;

}affpair;
*/
typedef struct{
	char description[500];
	char type[100];
	affndx ndx[FFtypeNum];
	affatomtype atomtype[FFtypeNum];
	afftrans trans[FFtypeNum];
    affbond bond[FFbondNum];
    affangle angle[FFangleNum];
    affdihedral dihedral[FFdihNum];
    affimproper improper[FFimNum];
    char affpname[100];
    int nndx, natomtype, nmap, ntrans, nbond, nangle, ndihedral, nimproper;
}aFFP;



  void convert_ffp2affp(FFP *f, aFFP *af);
  void print_aff(char *file, aFFP *af);






#endif
