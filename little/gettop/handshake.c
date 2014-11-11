
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools/rdiotools.h"

#include "utils/utils.h"
#include "graphic.h"
#include "opffdb.h"
#include "handshake.h"

 void init_itp_behavior(char *file, ItpControl *itpc)
 {
   FILE *fp;

   fp = fopen(file, "r");

   //deal with atomtypes section
   JmpOrNot(fp, '#');
   fgets(itpc->atomtype_cite, 500, fp);
   //atomtypes section : type part
   JmpOrNot(fp, '$');
   fscanf(fp, "%s%d", itpc->atomtype_name_prefix, &itpc->atomtype_name_start);
   //atomtypes section : ptype
   JmpOrNot(fp, '$');
   fscanf(fp, "%s", itpc->atomtype_ptype);


   //deal with mol section
   JmpOrNot(fp, '#');
   fgets(itpc->cite, 500, fp);
   fscanf(fp, "%s%d", itpc->name, &itpc->nrexcl);

   //deal with atoms section
   JmpOrNot(fp, '#');
   fgets(itpc->atom_cite, 500, fp);
   //atoms section : resnr
   JmpOrNot(fp, '$');
   fscanf(fp, "%d", &itpc->atom_resnr);
   //atoms section : cgnr
   JmpOrNot(fp, '$');
   fscanf(fp, "%d", &itpc->atom_cgnr);
   //atoms section : qtot
   JmpOrNot(fp, '$');
   fscanf(fp, "%d", &itpc->atom_quote_qtot);

   //deal with bonds section
   JmpOrNot(fp, '#');
   fgets(itpc->bond_cite, 500, fp);
   //bonds section : funct part
   JmpOrNot(fp, '$');
   fscanf(fp, "%d", &itpc->bond_funct);

   //deal with pairs section
   JmpOrNot(fp, '#');
   fgets(itpc->pair_cite, 500, fp);
   //pairs section: funct part
   JmpOrNot(fp, '$');
   fscanf(fp, "%d", &itpc->pair_funct);

   //deal with angles section
   JmpOrNot(fp, '#');
   fgets(itpc->angle_cite, 500, fp);
   //angles section: funct part
   JmpOrNot(fp, '$');
   fscanf(fp, "%d", &itpc->angle_funct);

   //deal with dihedrals section
   JmpOrNot(fp, '#');
   fgets(itpc->dihedral_cite, 500, fp);
   //dihedrals section: funct part
   JmpOrNot(fp, '$');
   fscanf(fp, "%d", &itpc->dihedral_funct);


   JmpOrNot(fp, '#');
   JmpOrNot(fp, '$');
   fscanf(fp, "%d", &itpc->dihedral_multi);

   fclose(fp);

   return;
 }

 void print_itp_behavior(char *file, ItpControl *itpc)
 {
   FILE *fp;

   fp = fopen(file, "w");

   fputs("#atomtypes", fp);
   fputs(itpc->atomtype_cite, fp);

   fprintf(fp, "$type\n%s  %d\n", itpc->atomtype_name_prefix, itpc->atomtype_name_start);
   //atomtypes section : ptype

   fprintf(fp, "$ptype\n%s\n", itpc->atomtype_ptype);


   //

   fclose(fp);

   return;

 }

 void gen_itp_www(ItpControl *itpc, NOBITPFILE *itp)
 {
   strcpy(itp->mo.cite, itpc->cite);
   strcpy(itp->mo.name, itpc->name);
   itp->mo.nrexcl = itpc->nrexcl;

   strcpy(itp->mo.atomtype_cite, itpc->atomtype_cite);

   // section...atom bond angle
   strcpy(itp->mo.atom_cite, itpc->atom_cite);
   strcpy(itp->mo.pair_cite, itpc->pair_cite);
   strcpy(itp->mo.bond_cite, itpc->bond_cite);
   strcpy(itp->mo.angle_cite, itpc->angle_cite);
   strcpy(itp->mo.dihedral_cite, itpc->dihedral_cite);


   return;
 }




 void con_com2ff(analyzeGraphic *com, aFFP *af)
 {
   int nndx, n;
   int i, j;

   nndx = af->nndx;
   n = com->nvex;

   for(i = 0; i < n; i++)
   {
     for(j = 0; j < nndx; j++)
	 {
	   if(!strcmp(com->vex[i].ffname, af->ndx[j].ffname))
	   {
	     com->vex[i].inode = af->ndx[j].inode;
		 break;
	   }
	 }
   }

   return;
 }

//atomtype_init;;;Î´Íê³É¡£¡£¡£
 int is_in_array(int *a, int len, int check)
 {
   int i;
   int flag = 0;

   for(i = 0; i < len; i++)
   {
     if(check == a[i])
	 {
	   flag = 1;
	   break;
	 }
   }
   return flag;
 }

 int gen_node_from_com(analyzeGraphic *com, int *a)
 {
   int i;
   int k, n;

   n = com->nvex;
   k = 0;

   for(i = 0; i < n; i++)
   {
      if(!is_in_array(a, k, com->vex[i].inode))
	  {
	     a[k] = com->vex[i].inode;
		 k++;
	  }     
   }
   return k;
 }


 void itp_atomtypes_section_init(analyzeGraphic *com, aFFP *af, ItpControl *itpc, NOBITPFILE *itp)
 {
   int i, n;
   int a[100], len;
   char element_name[10];
   affatomtype *atomtype;

   n = com->nvex;

   //generate from com...
   len = gen_node_from_com(com, a);
   itp->natomtype = len;
   
   for(i = 0; i < len; i++)
   {
     //automation  label....  id
     itp->atomtype[i].id = i;

	 //generate from com ..&ag...,  inode
	 itp->atomtype[i].inode = a[i];
     //bondtype
	 get_atomtype_name(af, a[i], itp->atomtype[i].bond_type);

	 //from control itp file.... opls_*  name, ptype
	 sprintf(itp->atomtype[i].name, "%s_%d", 
		 itpc->atomtype_name_prefix, itpc->atomtype_name_start+i+1);
	 strcpy(itp->atomtype[i].ptype, itpc->atomtype_ptype);

	 //from ag.... nuclear charge....
     get_element_from_com(com, a[i], element_name);
     itp->atomtype[i].ncharge = get_element_order(element_name);

	 //from db : charge, mass, sigma, eps
	 atomtype = (affatomtype *)ff_db_select(af, ATOMS_TABLE, &a[i]);

	 itp->atomtype[i].charge = atomtype->charge;
	 itp->atomtype[i].mass = atomtype->mass;
	 itp->atomtype[i].sig = atomtype->sig;
	 itp->atomtype[i].eps = atomtype->eps;
   } 

   return;
 }
 
 void get_typename_from_itp(NOBITPFILE *itp, int inode, char *type)
 {
   int i, ntype;

   ntype = itp->natomtype;

   for(i = 0; i < ntype; i++)
   {
      if(inode == itp->atomtype[i].inode)
	  {
	    strcpy(type, itp->atomtype[i].name);
		break;
	  }
   }
   
   return;
 }


 void itp_atoms_section_init(analyzeGraphic *com, aFFP *af, ItpControl *itpc, NOBITPFILE *itp)
 {
   int natom, i;
   char type[20];
   int inode ;
   affatomtype *atomtype;
   decimal qtot = 0.0;

   natom = com->nvex;
   itp->natom = natom;
   for(i = 0; i < natom; i++)
   {
      //determined by com... atom name, & nr(means id of atom)
      strcpy(itp->atom[i].atom, com->vex[i].unique_name);
	  itp->atom[i].nr = com->vex[i].id;

      //i don't think it's useful, but maybe...
	  inode = com->vex[i].inode;
	  itp->atom[i].inode = inode;

      //determined by db==af para... charge & mass of the atoms.
	  atomtype = (affatomtype*)ff_db_select(af, ATOMS_TABLE, &inode);
	  itp->atom[i].charge = atomtype->charge;
	  itp->atom[i].mass = atomtype->mass;
	  if(itpc->atom_quote_qtot == 1)
	  {
	     qtot += atomtype->charge;
	  }
	  itp->atom[i].qtot = qtot;

	  sprintf(itp->atom[i].quote, "qtot %lf", qtot);

	  //determined by itpcontrol behavior ...
	  //sprintf(type, "%s_%d", itpc->atomtype_name_prefix, itpc->atomtype_name_start+i+1);
      get_typename_from_itp(itp, inode, type);
	  strcpy(itp->atom[i].type, type);

	  strcpy(itp->atom[i].residue, itpc->name);
	  itp->atom[i].resnr = itpc->atom_resnr;
	  itp->atom[i].cgnr = itpc->atom_cgnr;
   }

   return ;

 }

 


 void itp_bonds_section_init(analyzeGraphic *com, aFFP *af, ItpControl *itpc, NOBITPFILE *itp)
 {
   int nbond, i, j, npara;
   int node[2], tnode[2];
   affbond *bond;
   char quote[100];
   //affndx *atom1;

   nbond = com->nbond;
   itp->nbond = nbond;

   for(i = 0; i < nbond; i++)
   {
	   //determined by com .. connection... ai,aj
     itp->bond[i].ai = com->bond[i].i;
	 itp->bond[i].aj = com->bond[i].j;

	 //determined by itpc... control items  funct
	 itp->bond[i].funct = itpc->bond_funct;

	 //determined by db af....  para, & npara
	 //node form com... 
	 node[0] = get_inode_form_id(com, com->bond[i].i);
	 node[1] = get_inode_form_id(com, com->bond[i].j);
     //select form db...
	 translate_db_inode(tnode, af, node, 2);
	 bond = (affbond*)ff_db_select(af, BONDS_TABLE, tnode);
	 if(bond != NULL)
	 {
	    itp->bond[i].npara = bond->npara;
	    npara = 2;
	    for(j = 0; j < npara; j++)
		{
	       itp->bond[i].para[j] = bond->para[j];
		}

	    //quote for human readable...link from db af
	    get_ffatom_names(af, node, 2, quote, "-");
	    sprintf(itp->bond[i].quote, "%s", quote);   
	 }
	 else
	 {
	    itp->bond[i].npara = 2;
	    npara = 2;
	    for(j = 0; j < npara; j++)
		{
	       itp->bond[i].para[j] = -1;
		}

	    //quote for human readable...link from db af
	    //get_ffatom_names(af, node, 2, quote, "-");
		strcpy(quote, "");
	    sprintf(itp->bond[i].quote, "%s", quote);   
		
		//warning...
		printf("warning: bond between %d %d parameter can't be found in db\n well, \
			ignore the warning or exit[y/n]\n", 
			itp->bond[i].ai, itp->bond[i].aj);
	    to_exist("enter 'y' to continue or press other key to exit\n", BOND_EXIT, 'y');
	 }
   }

   return;
 }



 void itp_pairs_section_init(analyzeGraphic *com, aFFP *af, ItpControl *itpc, NOBITPFILE *itp)
 {
   int i, npair;
   int node[2];
   char quote[100];
   //affpair *pair;

   npair = com->npair;
   itp->npair = npair;

   for(i = 0; i < npair; i++)
   {
	   //determined by com .. connection... ai,aj
      itp->pair[i].ai = com->pair[i].i;
	  itp->pair[i].aj = com->pair[i].l;

	 //determined by itpc... control items  funct
	 itp->pair[i].funct = itpc->pair_funct;

	 //determined by db af....  para, & npara
	 //node form com... 
	 node[0] = get_inode_form_id(com, com->pair[i].i);
	 node[1] = get_inode_form_id(com, com->pair[i].l);
     //select form db...
	 //pair = (affpair*)ff_db_select(af, 'PAIRS_TABLE', node);
	 //itp->pair[i].npara = pair->npara;
/*	 for(j = 0; j < npara; j++)
	 {
	   itp->pair[i].para[j] = pair->para[j];
	 }
*/
	 //quote for human readable...link from db af
	 get_ffatom_names(af, node, 2, quote, "--");
	 sprintf(itp->pair[i].quote, "%s", quote);   

   }

   return ;
 }


 void itp_angles_section_init(analyzeGraphic *com, aFFP *af, ItpControl *itpc, NOBITPFILE *itp)
 {
   int i, j, nangle, npara;
   int node[3], tnode[3];
   affangle *angle;
   char quote[100];

   nangle = com->nangle;
   itp->nangle = nangle;

   for(i = 0; i < nangle; i++)
   {
	   //determined by com .. connection... ai,aj
      itp->angle[i].ai = com->angle[i].i;
	  itp->angle[i].aj = com->angle[i].j;
	  itp->angle[i].ak = com->angle[i].k;

	 //determined by itpc... control items  funct
	 itp->angle[i].funct = itpc->angle_funct;

	 //determined by db af....  para, & npara
	 //node form com... 
	 node[0] = get_inode_form_id(com, com->angle[i].i);
	 node[1] = get_inode_form_id(com, com->angle[i].j);
	 node[2] = get_inode_form_id(com, com->angle[i].k);

     //select form db...
	 translate_db_inode(tnode, af, node, 3);
	 angle = (affangle*)ff_db_select(af, ANGLES_TABLE, tnode);
	 if(angle != NULL)
	 {
	    itp->angle[i].npara = angle->npara;
	    npara = angle->npara; 
	    for(j = 0; j < npara; j++)
		{
	       itp->angle[i].para[j] = angle->para[j];
		}

	    //quote for human readable...link from db af
	    get_ffatom_names(af, node, 3, quote, "-");
	    sprintf(itp->angle[i].quote, "%s", quote);   

	 }
     else
	 {
	    itp->angle[i].npara = 2;
	    npara = 3; 
	    for(j = 0; j < npara; j++)
		{
	       itp->angle[i].para[j] = -1;
		}

	    //quote for human readable...link from db af
	    //get_ffatom_names(af, node, 3, quote, "-");
		strcpy(quote, "");
	    sprintf(itp->angle[i].quote, "%s", quote); 

		//warning...
		printf("warning: angle between %d %d %d parameter can't be found in db\n well, \
			ignore the warning or exit\n", 
			itp->angle[i].ai, itp->angle[i].aj, itp->angle[i].ak);
		
		to_exist("enter 'y' to continue or press other key to exit\n", ANGLE_EXIT, 'y');

	 }
   }

   return ;

 }


 void itp_dihedrals_section_init(analyzeGraphic *com, aFFP *af, ItpControl *itpc, NOBITPFILE *itp)
 {
   int i, j, ndihedral, npara;
   int node[4], tnode[4];
   affdihedral *dihedral;
   char quote[100];

   ndihedral = com->ndihedral;
   itp->ndihedral = ndihedral;

   for(i = 0; i < ndihedral; i++)
   {
	   //determined by com .. connection... ai,aj
      itp->dihedral[i].ai = com->dihedral[i].i;
	  itp->dihedral[i].aj = com->dihedral[i].j;
	  itp->dihedral[i].ak = com->dihedral[i].k;
	  itp->dihedral[i].al = com->dihedral[i].l;

	 //determined by itpc... control items  funct
	 itp->dihedral[i].funct = itpc->dihedral_funct;
	 //itp->dihedral[i].multi = itpc->dihedral_funct;

	 //determined by db af....  para, & npara
	 //node form com... 
	 node[0] = get_inode_form_id(com, com->dihedral[i].i);
	 node[1] = get_inode_form_id(com, com->dihedral[i].j);
	 node[2] = get_inode_form_id(com, com->dihedral[i].k);
	 node[3] = get_inode_form_id(com, com->dihedral[i].l);

     //select form db...
	 translate_db_inode(tnode, af, node, 4);
	 dihedral = (affdihedral*)ff_db_select(af, DIHEDRALS_TABLE, tnode);
	 if(dihedral != NULL)
	 {
	    itp->dihedral[i].npara = dihedral->npara;
	    npara = dihedral->npara;
	    for(j = 0; j < npara; j++)
		{
	       itp->dihedral[i].para[j] = dihedral->para[j];
		}

	    //quote for human readable...link from db af
	    get_ffatom_names(af, node, 4, quote, "-");
	    sprintf(itp->dihedral[i].quote, "%s", quote); 

		if(itpc->dihedral_multi > 0)
		{
		   itp->dihedral[i].multi = dihedral->multi;
		}
		else
		{
		   itp->dihedral[i].multi = -1;
		}
		

	 }
	 else
	 {
	    itp->dihedral[i].npara = 5;
	    npara = 5;
	    for(j = 0; j < npara; j++)
		{
	       itp->dihedral[i].para[j] = -1;
		}

	    //quote for human readable...link from db af
	    //get_ffatom_names(af, node, 4, quote, "-");
		strcpy(quote, "");
	    sprintf(itp->dihedral[i].quote, "%s", quote);  


		printf("warning: DIHEDRAL between ATOMID %d %d %d %d parameter can't be found in db\n well, \
			ignore the warning or exit[y/n]:enter 'y' to continue\n", 
			itp->dihedral[i].ai, itp->dihedral[i].aj, itp->dihedral[i].ak, itp->dihedral[i].al);

		    to_exist("enter 'y' to continue or press other key to exit\n", DIHEDRAL_EXIT, 'y');

	 }
   }

   return ;
 }



 void itp_impropers_section_init(analyzeGraphic *com, aFFP *af, ItpControl *itpc, NOBITPFILE *itp)
 {
   return;
 }


 void gen_itp(analyzeGraphic *com, aFFP *af, ItpControl *itpc, NOBITPFILE *itp)
 {

	gen_itp_www(itpc, itp);
    itp_atomtypes_section_init(com, af, itpc, itp);
    itp_atoms_section_init(com, af, itpc, itp);
    itp_bonds_section_init(com, af, itpc, itp);
    itp_pairs_section_init(com, af, itpc, itp);
    itp_angles_section_init(com, af, itpc, itp);
    itp_dihedrals_section_init(com, af, itpc, itp);

    itp_impropers_section_init(com, af, itpc, itp);

    return;
 }

