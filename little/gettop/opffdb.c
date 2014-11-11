
#include<stdio.h>
#include<string.h>
#include<stdlib.h>

//#include "readff.h"
#include "opffdb.h"

 void *ff_db_select(aFFP *af, int select_table, int *node)
 {
   int mytable;
   void *p;

   mytable = select_table;

   switch(mytable)
   {
   case ATOMS_TABLE:
	   p = ff_db_select_atoms(af, node[0]);
	   break;

   case TRANSLATIONS_TABLE:
	   p = ff_db_select_translations(af, node[0]);
	   break;

   case BONDS_TABLE:
	   p = ff_db_select_bonds(af, node);
	   break;

   case ANGLES_TABLE:
	   p = ff_db_select_angles(af, node);
	   break;

   case DIHEDRALS_TABLE:
	   p = ff_db_select_dihedrals(af, node);
	   break;

   case IMPROPERS_TABLE:
	   p = ff_db_select_impropers(af, node);
	   break;

   case INDEX_TABLE:
	   p = ff_db_select_index(af, node[0]);

   default:
	   printf("select failed, action ignored. press any key to continue...\n");
	   system("PAUSE");   
   }

   return p;
 }



 void *ff_db_select_atoms(aFFP *af, int inode)
 {
   int i;

   for(i = 0; i < af->natomtype; i++)
   {
     if(af->atomtype[i].inode == inode)
		 break;
   }

   return (void*)((&af->atomtype[i]));
 }


#define iseq(a, b) if(a == b)


 void *ff_db_select_translations(aFFP *af, int inode)
 {
   int i, n;
   void *p = NULL;

   n = af->ntrans;

   for(i = 0; i < n; i++)
   {
     iseq(inode, af->trans[i].isource)
	 {
	   break;
	 }
   }

   if(i < n)
      p = (void*)(&(af->trans[i]));

   return p;
 }



  int compare_head(int *a, int *b, int n)
  {
    int i;
	int flag;

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
       
       
  int compare_tail(int *a, int *b, int n)
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
       

#define Loop(i, n) for(i = 0; i < n; i++){
#define endLoop }
 void *ff_db_select_bonds(aFFP *af, int *node)
 {
   int i, n;
   void *p = NULL;

   n = af->nbond;

   for(i = 0; i < n; i++)
   {
	   if(compare_head(af->bond[i].node, node, 2) || compare_tail(af->bond[i].node, node, 2))
		   break;
   }
  
   if(i < n)
	   p = (void*)(&(af->bond[i]));   

   return p;

 }


 void *ff_db_select_angles(aFFP *af, int *node)
 {
   int i, n;
   void *p = NULL;

   n = af->nangle;

   Loop(i, n)
	   if(compare_head(af->angle[i].node, node, 3) || compare_tail(af->angle[i].node, node, 3))
		   break;
   endLoop
  
   if(i < n)
	   p = (void*)(&(af->angle[i]));   

   return p;

 }


#define Node(x) af->x[i].node
 void *ff_db_select_dihedrals(aFFP *af, int *node)
 {
   int i, n;
   void *p = NULL;

   n = af->ndihedral;

   Loop(i, n)
	   if(compare_head(Node(dihedral), node, 4) || compare_tail(Node(dihedral), node, 4))
		   break;
   endLoop
  
   if(i < n)
	   p = (void*)(&(af->dihedral[i]));   

   return p;

 }
#undef Node


 void *ff_db_select_impropers(aFFP *af, int *node)
 {
   return NULL;
 }


 void *ff_db_select_index(aFFP *af, int inode)
 {
   return (void*)((&af->ndx[inode]));
 }




 static void get_ff_atom_name(aFFP *af, int inode, char *name)
 {
   int n;

   n = af->nndx;

   strcpy(name, af->ndx[inode].ffname);

   return;
 }


 void get_ffatom_names(aFFP *af, int *node, int n, char *line, char *link)
 {
   char name[10], tmp[15];
   //char *line;
   int i;

   //line = (char *)malloc(100*sizeof(char));
   strcpy(line, "");

   for(i = 0; i < n; i++)
   {
     get_ff_atom_name(af, node[i], name);
	 if(i < n-1)
	   sprintf(tmp, "%s%s", name, link);
	 else
	   sprintf(tmp, "%s", name);

	 strcat(line, tmp);
   }

   return ;
 }

 void get_atomtype_name(aFFP *af, int inode, char *name)
 {
   int i;

   for(i = 0; i < af->nndx; i++)
   {
     if(inode == af->ndx[i].inode)
	 {
	   strcpy(name, af->ndx[i].ffname);
	   break;
	 }
   }

   return;
 }

/*
 in: af, node, n
 out: tnode
 */
 void translate_db_inode(int *tnode, aFFP *af, int *node, int n)
 {
   int i, j, is_trans;


   for(i = 0; i < n; i++)
   {
	  is_trans = 0;
      for(j = 0; j < af->ntrans; j++)
	  {
	    if(af->trans[j].isource == node[i])
		{
		  is_trans = 1;
		  break;
		}
	  }
	  tnode[i] = is_trans?af->trans[j].idnst:node[i];
   }

   return;
 }
