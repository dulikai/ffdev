

#ifndef _HANDSHAKE_H_
#define _HANDSHAKE_H_

#include "graphic.h"
#include "ffdb.h"

#define ITPATOMNUM 100
#define ITPBONDNUM 200
#define ITPANGLENUM 500
#define ITPDIHEDRALNUM 1000
#define ITPIMPROPERNUM 100
#define ITPPAIRNUM 1000

#define ITPTYPENUM 50


//#define OPLS_init 1000

/*
structure of data in itp file.
*/
typedef struct{
        char cite[500];
        char name[50];
        int nrexcl;

		char atomtype_cite[500];
		char atomtype_name_prefix[20];
		int atomtype_name_start;
		char atomtype_ptype[10];

		char atom_cite[500];
		int atom_resnr;
		int atom_cgnr;
		int atom_quote_qtot;

		char bond_cite[500];
		int bond_funct;

		char pair_cite[500];
		int pair_funct;

		char angle_cite[500];
		int angle_funct;

		char dihedral_cite[500];
		int dihedral_funct;
		int dihedral_multi;

}ItpControl;

typedef struct{
        char cite[500];
        char name[50];
        int nrexcl;

		char atomtype_cite[500];
		//char atomtype_name_prefix[20];
		//int atomtype_name_start;
		//char atomtype_ptype[10];

		char atom_cite[500];
		//int atom_resnr;
		//int atom_cgnr;
		//int atom_quote_qtot;

		char bond_cite[500];
		//int bond_funct;

		char pair_cite[500];
		//int pair_funct;

		char angle_cite[500];
		//int angle_funct;

		char dihedral_cite[500];
		//int dihedral_funct;

        }moleculetype;        
        
typedef struct{
        //char cite[500];
        int nr;
		int inode;  
        char type[20];
        int resnr;
        char residue[20];
        char atom[15];
        int cgnr;
        decimal charge, mass;
		decimal qtot;
        char quote[50];
        }atoms;      
        
typedef struct{
        //char cite[500];
        int id;
		int inode;
		char name[20];
        char bond_type[20];
		int ncharge;
		char ptype[10];
        decimal charge, mass;
		decimal sig, eps;
        }atomtypes;      


typedef struct{
        int ai;
        int aj;
        int funct;
        decimal para[10];
		int npara;
        char quote[50];
        }bonds;         
        
typedef struct{
        int ai;
        int aj;
        int funct;
        //decimal para[10];
		//int npara;
        char quote[50];
        }pairs;
        
typedef struct{
        int ai;
        int aj;
        int ak;
        int funct;
        decimal para[10];
		int npara;
        char quote[50];
        }angles;                 

typedef struct{
        int ai;
        int aj;
        int ak;
        int al;
        int funct;
		int multi;
        decimal para[10];
		int npara;
        char quote[50];
        }dihedrals;
        
typedef struct{
        int ai, aj, ak, al;
        int funct;
        decimal para[10];
		int npara;
        char quote[50];
        }impropers;        


typedef struct{
        moleculetype mo;
        atomtypes atomtype[ITPTYPENUM];
        atoms atom[ITPATOMNUM];
        bonds bond[ITPBONDNUM];
        pairs pair[ITPPAIRNUM];
        angles angle[ITPANGLENUM];
        dihedrals dihedral[ITPDIHEDRALNUM];
        impropers improper[ITPIMPROPERNUM];
		int natom, nbond, npair, nangle, ndihedral, nimproper, natomtype;
        }NOBITPFILE;
        


 void con_com2ff(analyzeGraphic *com, aFFP *af);

 void init_itp_behavior(char *file, ItpControl *itpc);
 void print_itp_behavior(char *file, ItpControl *itpc);

 void gen_itp(analyzeGraphic *com, aFFP *af, ItpControl *itpc, NOBITPFILE *itp);





#endif  /* _HANDSHAKE_H_ */
