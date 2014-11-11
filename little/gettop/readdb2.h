

#ifndef _READDB2_H_
#define _READDB2_H_

#include "utils/mytype.h"

#define AA_NUM 30
#define AA_NDX_NUM 100
#define AA_ATOMCNST_NUM 50
#define AA_ATOMTYPE_NUM 50
#define AA_TRANS 100

#define AA_BONDCNST_NUM 100
#define AA_BOND_NUM 50

#define AA_ANGLECNST_NUM 100
#define AA_ANGLE_NUM 50

#define AA_DIHEDRAL_NUM 1000



typedef struct{
	char ffname[20];
	int inode;
}aa_pubndx;

typedef struct{
	int id;
	char wffname[20];
	decimal sig, eps;
	decimal mass;
	char source[50];
}aa_atomcnst;

typedef struct{
	int id;
	char ffname[20];
	decimal charge;
}aa_atomtype;

typedef struct{
	int id;
	char aaname[20];
	int natomtype;
	aa_atomtype atomtype[AA_ATOMTYPE_NUM];
}aa_atomvar;


typedef struct{
     int id;
     char dnst[20];
     char source[20];
}aa_trans;

typedef struct{
	int id;
    int inode, jnode;
	decimal kr;
	char source[50];
}aa_bondcnst;

typedef struct{
	int id;
	char catom[20];
	char datom[20];
	int inode, jnode;
	decimal r0;
}aa_bond;

typedef struct{
	int id;
	char aaname[20];
	int nbond;
	aa_bond bond[AA_BOND_NUM];
}aa_bondvar;

typedef struct{
	int id;
    int inode, jnode, knode;
	decimal th;
	char source[50];
}aa_angle;


typedef struct{
	int id;
    int inode, jnode, knode;
	decimal kth;
	char source[50];
}aa_anglecnst;

typedef struct{
	int id;
	int nangle;
	char aaname[20];
	aa_angle angle[AA_ANGLE_NUM];
}aa_anglevar;


typedef struct{
	int id;
    int inode, jnode, knode, lnode;
	decimal gam, kphi;
	int multi;
	char source[50];
}aa_dihedral;

typedef struct{
	char aaname[20];
	int id;
}AAseries;


typedef struct{
    AAseries name[AA_NUM];
	aa_pubndx ndx[AA_NDX_NUM];
	aa_atomcnst atomcnst[AA_ATOMCNST_NUM];
	aa_atomvar atomvar[AA_NUM];
	aa_trans trans[AA_TRANS];
	aa_bondcnst bondcnst[AA_BONDCNST_NUM];
	aa_bondvar bondvar[AA_NUM];
	aa_anglecnst anglecnst[AA_ANGLECNST_NUM];
	aa_anglevar anglevar[AA_NUM];
	aa_dihedral dihedral[AA_DIHEDRAL_NUM];
	int nname, nndx;
	int natomcnst, natomvar, ntrans, nbondcnst, nbondvar, nanglecnst, nanglevar, ndihedral;
}AAFFP;



 void read_aaff_db(char *file, AAFFP *aa);




#endif  /* _READDB2_H_ */