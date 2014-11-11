
#include "ffpara.h"
#include "ffdb2.h"

 void gen_ffpara(FILEControl *fc, aFFP *af)
 {
    FFP *f;
	AAFFP aa;

	if(!stricmp(fc->dbtype, "A"))
	{
       f = init_ffp();
       read_ff_db(fc->dbfile, f);
       //print_ff_db("outff.txt", f); 
   
       convert_ffp2affp(f, af);
	}
	else if(!stricmp(fc->dbtype, "B"))
	{
       read_aaff_db(fc->dbfile, &aa);

       convert_aaffp2afp(&aa, af, fc->dbpara);
	}
	else
	{
	  printf("error occured : db type unrecognize, sorry, ...\nyour dbtype: %s\n", fc->dbtype);
	  exit(1);
	}
   //print_aff("aff.txt", af);

    return;
             
 }



 void print_aff(char *file, aFFP *af)
 {
   FILE *fp;
   int i;

   fp = fopen(file, "w");

   fputs(af->description, fp);

   fputs("\n\n", fp);
   fputs(af->type, fp);

  
   /* node index table */
   fprintf(fp, "\n\nnode index table: number %2d\n", af->nndx);
   for(i = 0; i < af->nndx; i++)
   {
     fprintf(fp, "inode %d : %s\n", af->ndx[i].inode, af->ndx[i].ffname);
   }


   /*
   the definition of atomtype
   */
   fprintf(fp, "\n\nthe definition of atomtype\n");
   fprintf(fp, "inode mass charge eps sig\n");
   for(i = 0; i < af->natomtype; i++)
   {
     fprintf(fp, "%d %10.6lf %10.6lf %10.6lf %10.6lf %s\n", af->atomtype[i].inode, af->atomtype[i].mass, 
		 af->atomtype[i].charge,af->atomtype[i].eps,
		 af->atomtype[i].sig, af->atomtype[i].field);
   }


   /* for translation section */
   fprintf(fp, "\n\nfor translation section\n");
   fprintf(fp, "id idnst isource\n");

   for(i = 0; i < af->ntrans; i++)
   {
	 fprintf(fp, "%d %d %d\n", af->trans[i].id, 
		 af->trans[i].idnst, af->trans[i].isource); 
   }

 /* for bond parameter */
   fprintf(fp, "\n\nfor bond parameter\n");
   fprintf(fp, "id inode jnode para1  para2 type\n");
   for(i = 0; i < af->nbond; i++)
   {
	 fprintf(fp, "%2d %2d %2d %10.6le %10.6le %s\n", af->bond[i].id, 
		 af->bond[i].node[INODE], af->bond[i].node[JNODE],
		 af->bond[i].para[FF_bond_PARA_RE], af->bond[i].para[FF_bond_PARA_KA],
		 af->bond[i].type);
   }

   /* for angle parameter */
   fprintf(fp, "\n\nfor angle parameter\n");
   fprintf(fp, "id inode jnode jnode para1  para2 type\n");
   for(i = 0; i < af->nangle; i++)
   {
	 fprintf(fp, "%2d  %2d %2d %2d %10.6le %10.6le %s\n", af->angle[i].id, 
		 af->angle[i].node[INODE], af->angle[i].node[JNODE], af->angle[i].node[JNODE],
		 af->angle[i].para[0], af->angle[i].para[1],
		 af->angle[i].type);
   }

   /* for dihedral parameter */
   fprintf(fp, "\n\nfor dihedral parameter\n");
   fprintf(fp, "id inode jnode jnode para1  para2 type\n");
   for(i = 0; i < af->ndihedral; i++)
   {
	 fprintf(fp, "%2d  %2d %2d %2d %2d %10.6le %10.6le %s\n", af->dihedral[i].id, 
		 af->dihedral[i].node[INODE], af->dihedral[i].node[JNODE], 
		 af->dihedral[i].node[KNODE],af->dihedral[i].node[LNODE],
		 af->dihedral[i].para[0], af->dihedral[i].para[1],
		 af->dihedral[i].type);
   }


   fclose(fp);

   return;
 }
