

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "utils/mytype.h"
#include "tools/rdiotools.h"
#include "tools/opstr.h"


#include "ffdb.h"


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

 void convert_ffp2affp(FFP *f, aFFP *af)
 {
   int i, multi;
   double c[10];

   /* these two section is ignored now */
   strcpy(af->description, 
"#read from\n# il.ff (kJ/mol, A, deg)\n# Molecular force field for ionic liquids\n\
# DL_POLY format except dihedral functions that may have 4 cosine terms.\
# A. Padua, J. Canongia Lopes, K. Shimizu, 12-05-2008\
# questions to: agilio.padua@univ-bpclermont.fr\
# latest version donwloadable at: http://therm10.univ-bpclermont.fr/~apadua/\
"
);
   strcpy(af->type, "FF for IL:unit[amu, e, kj/mol, nm]");

   /* node index table */
   af->nndx = f->nndx;
   for(i = 0; i < f->nndx; i++)
   {
     af->ndx[i].inode = f->ndx[i].inode;
	 strcpy(af->ndx[i].ffname, f->ndx[i].ffname);
   }

   /* the definition of atomtype */
   af->natomtype = f->natomtype;
   for(i = 0; i < f->natomtype; i++)
   {
	 af->atomtype[i].inode = f->atomtype[i].inode;

     af->atomtype[i].charge = f->atomtype[i].charge;
	 af->atomtype[i].eps = f->atomtype[i].eps;
	 af->atomtype[i].mass = f->atomtype[i].mass;

	 //unit switch happened
	 af->atomtype[i].sig = f->atomtype[i].sig * 0.1;  //unit: A-->nm
	 strcpy(af->atomtype[i].field, f->atomtype[i].field);
   }

   /* for translation section */
   af->ntrans = f->nalias;
   for(i = 0; i < f->nalias; i++)
   {
     af->trans[i].id = f->alias[i].id;

	 af->trans[i].idnst = f->alias[i].idnst;
	 af->trans[i].isource = f->alias[i].isource;
   }


 /* for pair section */

 /* for bond parameter */
   af->nbond = f->nbond;
   for(i = 0; i < f->nbond; i++)
   {
     af->bond[i].id = f->bond[i].id;

	 af->bond[i].node[INODE] = f->bond[i].inode;
	 af->bond[i].node[JNODE] = f->bond[i].jnode;
	 af->bond[i].npara = 2;

	 //unit switch happened
	 af->bond[i].para[FF_bond_PARA_RE] = f->bond[i].re * 0.1;  //A-->nm
	 af->bond[i].para[FF_bond_PARA_KA] = fabs (f->bond[i].ka) * 100;  //kj*mol-1*A-2 --> kj*mol-1*A-2; also neglect 'negative means a rigid constraint'

	 strcpy(af->bond[i].type, f->bond[i].type);
   }

   /* for angle parameter */
   af->nangle = f->nangle;
   for(i = 0; i < f->nangle; i++)
   {
     af->angle[i].id = f->angle[i].id;

	 af->angle[i].node[INODE] = f->angle[i].inode;
	 af->angle[i].node[JNODE] = f->angle[i].jnode;
	 af->angle[i].node[KNODE] = f->angle[i].knode;

	 af->angle[i].npara = 2;

	 af->angle[i].para[FF_angle_PARA_TH] = f->angle[i].th;
	 af->angle[i].para[FF_angle_PARA_KA] = f->angle[i].ka;

	 strcpy(af->angle[i].type, f->angle[i].type);
   }

   /* for dihedral parameter */
   af->ndihedral = f->ndihedral;
   for(i = 0; i < f->ndihedral; i++)
   {
     af->dihedral[i].id = f->dihedral[i].id;

	 af->dihedral[i].node[INODE] = f->dihedral[i].inode;
	 af->dihedral[i].node[JNODE] = f->dihedral[i].jnode;
	 af->dihedral[i].node[KNODE] = f->dihedral[i].knode;
	 af->dihedral[i].node[LNODE] = f->dihedral[i].lnode;

	 af->dihedral[i].npara = 6;

	 get_dihedral_para(0.0, f->dihedral[i].v1, f->dihedral[i].v2, f->dihedral[i].v3, f->dihedral[i].v4, c);

	 af->dihedral[i].para[FF_dihedral_PARA_C0] = c[0];
	 af->dihedral[i].para[FF_dihedral_PARA_C1] = c[1];
	 af->dihedral[i].para[FF_dihedral_PARA_C2] = c[2];
	 af->dihedral[i].para[FF_dihedral_PARA_C3] = c[3];
	 af->dihedral[i].para[FF_dihedral_PARA_C4] = c[4];
	 af->dihedral[i].para[FF_dihedral_PARA_C5] = c[5];

	 strcpy(af->dihedral[i].type, f->dihedral[i].type);

     sscanf(f->dihedral[i].type, "cos%d", &multi);
	 af->dihedral[i].multi = multi;
   }

   af->nimproper = f->nimproper;
   for(i = 0; i < f->nimproper; i++)
   {
     af->improper[i].id = f->improper[i].id;

	 af->improper[i].node[INODE] = f->improper[i].inode;
	 af->improper[i].node[JNODE] = f->improper[i].jnode;
	 af->improper[i].node[KNODE] = f->improper[i].knode;
	 af->improper[i].node[LNODE] = f->improper[i].lnode;

	 af->improper[i].npara = 5;

	 get_dihedral_para(0.0, f->improper[i].v1, f->improper[i].v2, f->improper[i].v3, f->improper[i].v4, c);

	 af->improper[i].para[FF_improper_PARA_C0] = c[0];
	 af->improper[i].para[FF_improper_PARA_C1] = c[1];
	 af->improper[i].para[FF_improper_PARA_C2] = c[2];
	 af->improper[i].para[FF_improper_PARA_C3] = c[3];
	 af->improper[i].para[FF_improper_PARA_C4] = c[4];
	 af->improper[i].para[FF_improper_PARA_C5] = c[5];

	 strcpy(af->improper[i].type, f->improper[i].type);
   }

   return;

 }



