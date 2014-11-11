

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "utils/mytype.h"
#include "tools/rdiotools.h"
#include "tools/opstr.h"

#include "ffdb.h"
#include "ffdb2.h"



#define MYLOOP(type, entity) n = aa->n##type; \
 for(i = 0; i < n; i++) \
 { \
 if(strstr(aa->type[i].aaname, entity)) \
 { p = (void*) (&aa->type[i]); break;} \
 } \
if(i == n) \
 { \
  printf(" %s entity not exist in " #type" table...\ntime:%s %s\n", entity, __DATE__, __TIME__); \
  exit(1); \
 }


 void * select_aa_entity(AAFFP *aa, char *entity, int type)
 {
   int i, n;
   void *p;

   
   switch(type)
   {
   case AA_ATOMS_TABLE:
	   n = aa->natomvar;
	   for(i = 0; i < n; i++){
		   if(strstr(aa->atomvar[i].aaname, entity))
		   {
             p = (void*) (&aa->atomvar[i]);
		     break;
		   }
	   }
	   if(i == n)
	   {
	     printf(" %s entity not exist in atoms table...\n", entity);
		 exit(1);
	   }
	   
	   break;

   case AA_BONDS_TABLE:
	   MYLOOP(bondvar, entity);
	   break;

   case AA_ANGLES_TABLE:
	   MYLOOP(anglevar, entity);
	   break;

   default:
	   printf("invalid table come...\n");
	   p = NULL;
	   break;

   }
   return p;
 }
#undef MYLOOP


 void gen_aa2af_ndx(AAFFP *aa, aFFP *af, char *entity)
 {
   int i, nndx, j, natomtype;
   int flag = 1;
   aa_atomvar *var;

	 /* node index table */
   nndx = aa->nndx;
   for(i = 0; i < aa->nndx; i++)
   {
     af->ndx[i].inode = aa->ndx[i].inode;
	 strcpy(af->ndx[i].ffname, aa->ndx[i].ffname);
   }
   var = (aa_atomvar *)select_aa_entity(aa, entity, AA_ATOMS_TABLE);

   natomtype = var->natomtype;

   for(i = 0; i < natomtype; i++)
   {
     flag = 1;
     for(j = 0; j < nndx; j++)
	 {
	   if(!strcmp(var->atomtype[i].ffname, af->ndx[j].ffname))
	   {
	     flag = 0;
		 break;
	   }
	 }

	 if(flag)
	 {
	   strcpy(af->ndx[nndx].ffname, var->atomtype[i].ffname);
	   af->ndx[nndx].inode = nndx++;
	 }
   }
   af->nndx = nndx;

   return;   
 }

 void translate_aa_atomtype(char *atom, char *translation, AAFFP *aa)
 {
   int i, n;

   n = aa->ntrans;

   strcpy(translation, atom);

   for(i = 0; i < n; i++)
   {
     if(!strcmp(atom, aa->trans[i].source))
	 {
	   strcpy(translation, aa->trans[i].dnst);
	   break;
	 }
   }
   return;
 }


 void gen_aa2af_atomtype(AAFFP *aa, aFFP *af, char *entity)
 {
   int i, nndx, j, natomtype;
   char ffname[20];
   aa_atomvar *var;


	 /* node atom table */
   nndx = af->nndx;
   var = (aa_atomvar *)select_aa_entity(aa, entity, AA_ATOMS_TABLE);

   natomtype = var->natomtype;

   for(i = 0; i < natomtype; i++)
     for(j = 0; j < nndx; j++)
	 {
	   if(!strcmp(var->atomtype[i].ffname, af->ndx[j].ffname))
	   {
	     af->atomtype[i].inode = af->ndx[j].inode;
		 break;
	   }
	 }
	 af->natomtype = natomtype;
   
   //SWITCH HAPPED
	 for(i = 0; i < var->natomtype; i++)
	 {
       for(j = 0; j < aa->natomcnst; j++)
	   { 
        translate_aa_atomtype(var->atomtype[i].ffname, ffname, aa);
		if(!strcmp(ffname, aa->atomcnst[j].wffname))
		{
		   af->atomtype[i].sig = aa->atomcnst[j].sig  *  0.1;  //CONVERT angstroms TO nm
		   af->atomtype[i].eps = aa->atomcnst[j].eps;
		   af->atomtype[i].mass = aa->atomcnst[j].mass;
		   break;
		}

	   }
			af->atomtype[i].charge = var->atomtype[i].charge;
		    strcpy(af->atomtype[i].field, "AA");
	 }


   return;   
 
 }


   /* for translation section */
 
#define For(n)    for(i = 0; i < n; i++){
#define Endfor }

 static int get_inode_from_af2(aFFP *af, char *ffname)
 {
   int i;

   For(af->nndx)
	   if(!strcmp(af->ndx[i].ffname, ffname))
	   {
		   break;
	   }
   Endfor
   if(i > af->nndx)
   {
     printf("atomtype [%s] can't find inode...\n", ffname);
	 exit(0);
   }
   return i;
 }


 void gen_aa2af_translation(AAFFP *aa, aFFP *af)
 {
   int i;

   For(aa->ntrans)
   {
     af->trans[i].id = i;
	 af->trans[i].idnst = get_inode_from_af2(af, aa->trans[i].dnst);
	 af->trans[i].isource = get_inode_from_af2(af, aa->trans[i].source);
   }
   Endfor

   af->ntrans = aa->ntrans;

   return;
   
 }
#undef For
#undef Endfor


  static int compare_head(int a[], int b[], int n)
  {
    int i, flag;
	//int flag;

	flag = 1;

    for(i = 0; i < n; i++)
	{
	  if(a[i] != b[i])
	  {
	    flag = 0;
		break;
	  }
	}
	return flag;      
  }
       
       
  static int compare_tail(int *a, int *b, int n)
  {
    int i;
	int flag;

	flag = 1;

    for(i = 0; i < n; i++)
	{
	  if(a[i] != b[n-i-1])
	  {
	    flag = 0;
		break;
	  }
	}
	return flag;      
  }

 void gen_aa2af_bonds(AAFFP *aa, aFFP *af, char *entity)
 {
   int i, j;
   int node[2];
   aa_bondvar *var;

   var = (aa_bondvar *)select_aa_entity(aa, entity, AA_BONDS_TABLE);
   af->nbond = var->nbond;
   for(i = 0; i < var->nbond; i++)
   {
     af->bond[i].id = var->bond[i].id;

	 af->bond[i].node[INODE] = var->bond[i].inode;
	 af->bond[i].node[JNODE] = var->bond[i].jnode;
	 af->bond[i].npara = 2;

	 //unit switch happened
	 af->bond[i].para[AA_BOND_PARA_R0] = var->bond[i].r0 * 0.1;  //A-->nm

	 strcpy(af->bond[i].type, "");

	 //nodea[0] = var->bond[i].inode;
	 //nodea[1] = var->bond[i].jnode;


	 for(j = 0; j < aa->nbondcnst; j++)
	 {
       node[0] = aa->bondcnst[j].inode;
       node[1] = aa->bondcnst[j].jnode;

	   if(compare_head(af->bond[i].node, node, 2) || compare_tail(af->bond[i].node, node, 2))
	   {
		  //kj*mol-1*A-2 --> kj*mol-1*A-2; also neglect 'negative means a rigid constraint'
	      af->bond[i].para[1] = fabs (aa->bondcnst[j].kr) * 200;	  
		  
		  break;
	   }
	 }
   }
   return;
 }



 void gen_aa2af_angles(AAFFP *aa, aFFP *af, char *entity)
 {
   int i, j;
   int node[3];
   aa_anglevar *var;

   var = (aa_anglevar *)select_aa_entity(aa, entity, AA_ANGLES_TABLE);
   af->nangle = var->nangle;
   for(i = 0; i < var->nangle; i++)
   {
     af->angle[i].id = var->angle[i].id;

	 af->angle[i].node[INODE] = var->angle[i].inode;
	 af->angle[i].node[JNODE] = var->angle[i].jnode;
	 af->angle[i].node[KNODE] = var->angle[i].knode;
	 af->angle[i].npara = 2;

	 //unit switch happened
	 af->angle[i].para[AA_ANGLE_TH] = var->angle[i].th;

	 strcpy(af->angle[i].type, "");

	 //nodea[0] = var->angle[i].inode;
	 //nodea[1] = var->angle[i].jnode;


	 for(j = 0; j < aa->nanglecnst; j++)
	 {
       node[0] = aa->anglecnst[j].inode;
       node[1] = aa->anglecnst[j].jnode;
       node[2] = aa->anglecnst[j].knode;

	   if(compare_head(af->angle[i].node, node, 3) || compare_tail(af->angle[i].node, node, 3))
	   {
		  //kj*mol-1*rad-2 --> kj*mol-1*A-2; also neglect 'negative means a rigid constraint'
	      af->angle[i].para[AA_ANGLE_KTH] = fabs (aa->anglecnst[j].kth) * 2;
		  break;
		  //k*x^2 or 0.5*k*x^2
	   }
	 }
   }
   return;
 }

 void get_aa_dihedral_para(decimal gam, decimal kphi, decimal *c)
 {
   c[0] = gam;
   c[1] = kphi;

   return;
 }

 void gen_aa2af_dihedrals(AAFFP *aa, aFFP *af)
 {
   int i;
   decimal c[10];

   af->ndihedral = aa->ndihedral;

   for(i = 0; i < aa->ndihedral; i++)
   {
     af->dihedral[i].id = aa->dihedral[i].id;

	 af->dihedral[i].node[INODE] = aa->dihedral[i].inode;
	 af->dihedral[i].node[JNODE] = aa->dihedral[i].jnode;
	 af->dihedral[i].node[KNODE] = aa->dihedral[i].knode;
	 af->dihedral[i].node[LNODE] = aa->dihedral[i].lnode;

	 af->dihedral[i].npara = 2;

	 get_aa_dihedral_para(aa->dihedral[i].gam, aa->dihedral[i].kphi, c);

	 af->dihedral[i].para[AA_dihedral_PARA_C0] = c[0];
	 af->dihedral[i].para[AA_dihedral_PARA_C1] = c[1];
	 //af->dihedral[i].para[FF_dihedral_PARA_C2] = c[2];
	 //af->dihedral[i].para[FF_dihedral_PARA_C3] = c[3];
	 //af->dihedral[i].para[FF_dihedral_PARA_C4] = c[4];
	 //af->dihedral[i].para[FF_dihedral_PARA_C5] = c[5];

	 strcpy(af->dihedral[i].type, "AA");

     //sscanf(f->dihedral[i].type, "cos%d", &multi);
	 af->dihedral[i].multi = aa->dihedral[i].multi;
   }
 }



 void convert_aaffp2afp(AAFFP *aa, aFFP *af, char *entity)
 {

   /* these two section is ignored now */
   strcpy(af->description, 
"\
#amino acid parameter file \
#Molecular force field for ionic liquids \
#Guohui Zhoua,b, Xiaomin Liua,b, Suojiang Zhanga, Guangren Yua,b, Hongyan He \
#convert from the support information of this article : \
#A Force Field for Molecular Simulation of Tetrabutylphosphonium Amino Acid Ionic \
#J. Phys. Chem. B 2007, 111, 7078-7084 \
#added by dulikai \
#date 2008.12.25  \
"
);
   strcpy(af->type, "FF for AA:unit[angstroms,amu,kjmol-1,e,deg]");

   /* node index table */

   gen_aa2af_ndx(aa, af, entity);

   /* atoms table */
   gen_aa2af_atomtype(aa, af, entity);

   /*  */
   gen_aa2af_translation(aa, af);

 /* for bond parameter */
   gen_aa2af_bonds(aa, af, entity);


  /* for angle parameter */
   gen_aa2af_angles(aa, af, entity);


   /* for dihedral parameter */

   gen_aa2af_dihedrals(aa, af);
	   
	   
   return;

 }



