
#include "print.h"

 void print_declare(FILE *fp)
 {
   fprintf(fp, ";FIRST: generate automation\n");
   fprintf(fp, ";SECOND: improper not included default\n");
   fprintf(fp, ";THIRD: no responsibile to the correct of this top file\n");
   fprintf(fp, ";AT LAST: MAYBE IT'S USEFUL\n");
 }

 void print_itp(char *file, NOBITPFILE *itp)
 {
   FILE *fp;
   int i, j, n;


   fp = fopen(file, "w");

   print_declare(fp);
   fprintf(fp, "[ moleculetype ]\n");
   fprintf(fp, "%s", itp->mo.cite);
   fprintf(fp, "%s       %d\n", itp->mo.name, itp->mo.nrexcl);

   fprintf(fp, "\n");

   fprintf(fp, "[ atoms ]\n");
   fprintf(fp, "%s", itp->mo.atom_cite);

   //atoms section
   n = itp->natom;
   for(i = 0; i < n; i++)
   {
      fprintf(fp, "%5d  %12s  %2d %6s %6s %5d %12.5lf %12.5lf", itp->atom[i].nr, itp->atom[i].type, 
		  itp->atom[i].resnr, itp->atom[i].residue, 
		  itp->atom[i].atom, itp->atom[i].cgnr,
		  itp->atom[i].charge, itp->atom[i].mass);

	      fprintf(fp, "   ; qtot %.5lf\n", itp->atom[i].qtot);
   }

   //bonds section
   fprintf(fp, "\n[ bonds ]\n");
   fprintf(fp, "%s", itp->mo.bond_cite);
   n = itp->nbond;
   for(i = 0; i < n; i++)
   {
      fprintf(fp, "%4d  %4d  %2d       ", itp->bond[i].ai, itp->bond[i].aj, 
		  itp->bond[i].funct);
	  //printf("%5d  %5d  %2d ", itp->bond[i].ai, itp->bond[i].aj, itp->bond[i].funct);
	  for(j = 0; j < itp->bond[i].npara; j++)
	  {
	    fprintf(fp, "%.6le  ", itp->bond[i].para[j]);
		//printf("%.6le ", itp->bond[i].para[j]);
	  }

	  fprintf(fp, "   ; %s\n", itp->bond[i].quote);
	  //printf("   ; %s\n", itp->bond[i].quote);
   }

//pair section
   fprintf(fp, "\n[ pairs ]\n");
   fprintf(fp, "%s", itp->mo.pair_cite);
   n = itp->npair;
   for(i = 0; i < n; i++)
   {
      fprintf(fp, "%4d  %4d  %2d       ", itp->pair[i].ai, itp->pair[i].aj, 
		  itp->pair[i].funct);
	  //printf("%5d  %5d  %2d ", itp->pair[i].ai, itp->pair[i].aj, itp->pair[i].funct);
	 
	  /*for(j = 0; j < itp->pair[i].npara; j++)
	  {
	    fprintf(fp, "%.6le  ", itp->pair[i].para[j]);
		printf("%.6le ", itp->pair[i].para[j]);
	  }*/

	  fprintf(fp, "                           ; %s\n", itp->pair[i].quote);
	  //printf("   ; %s\n", itp->pair[i].quote);
   }


   //angles section 重复真讨厌！下次一定用模板，我保证。。。^_^ ^_^
   fprintf(fp, "\n[ angles ]\n");
   fprintf(fp, "%s", itp->mo.angle_cite);
   n = itp->nangle;
   for(i = 0; i < n; i++)
   {
      fprintf(fp, "%4d  %4d  %4d  %2d       ", itp->angle[i].ai, itp->angle[i].aj, itp->angle[i].ak,
		  itp->angle[i].funct);
	  //printf("%5d  %5d  %2d ", itp->angle[i].ai, itp->angle[i].aj, itp->angle[i].funct);
	  for(j = 0; j < itp->angle[i].npara; j++)
	  {
	    fprintf(fp, "%.6le  ", itp->angle[i].para[j]);
		//printf("%.6le ", itp->angle[i].para[j]);
	  }

	  fprintf(fp, "   ; %s\n", itp->angle[i].quote);
	  //printf("   ; %s\n", itp->angle[i].quote);
   }


   //dihedrals section 重复真讨厌！下次一定用模板，我保证。。。^_^ ^_^
   fprintf(fp, "\n[ dihedrals ]\n");
   fprintf(fp, "%s", itp->mo.dihedral_cite);
   n = itp->ndihedral;
   for(i = 0; i < n; i++)
   {
      fprintf(fp, "%4d  %4d  %4d  %4d  %2d       ", itp->dihedral[i].ai, itp->dihedral[i].aj, 
		  itp->dihedral[i].ak, itp->dihedral[i].al,
		  itp->dihedral[i].funct);
	  //printf("%5d  %5d  %2d ", itp->dihedral[i].ai, itp->dihedral[i].aj, itp->dihedral[i].funct);
	  for(j = 0; j < itp->dihedral[i].npara; j++)
	  {
	    fprintf(fp, "%10.6lf  ", itp->dihedral[i].para[j]);
		//printf("%.6le ", itp->dihedral[i].para[j]);
	  }
	  if(itp->dihedral[i].multi > 0)
	      fprintf(fp, "%5d ", itp->dihedral[i].multi);

	  fprintf(fp, "   ; %s\n", itp->dihedral[i].quote);
	  //printf("   ; %s\n", itp->dihedral[i].quote);
   }

   fprintf(fp, "\n; warning也许存在非正常二面角，but我没有生成improper dihedrals, 不想写了，自己加上吧...\n");
   fprintf(fp, ";warning improper dihedrals don't included, I don't want to do that, do it yourself\n");
   fprintf(fp, ";I can't promise this file is correct, thank you.\n;^_^\n");



   fclose(fp);

 }

  void print_itp_nob(char *file, NOBITPFILE *itp)
  {
   FILE *fp;
   int i, n;

   fp = fopen(file, "w");

   //print_declare(fp);

   n = itp->natomtype;

   fprintf(fp, "[ atomtypes ]\n");
   fprintf(fp, "%s", itp->mo.atomtype_cite);
   //fprintf(fp, "; %s\n", itp->mo.name);

   for(i = 0; i < n; i++)
   {
     fprintf(fp, "%10s %5s %2d %10.5lf %10.5lf %5s %13.5le %13.5le\n", itp->atomtype[i].name, itp->atomtype[i].bond_type,
		 itp->atomtype[i].ncharge,itp->atomtype[i].mass, itp->atomtype[i].charge,
		 itp->atomtype[i].ptype, itp->atomtype[i].sig, itp->atomtype[i].eps);
   }

   fclose(fp);
  }


 void print_top(FILEControl *fc, NOBITPFILE *itp)
 {
      
    print_itp(fc->outitpfile, itp);

    print_itp_nob(fc->itp_nobfile, itp);
    
    return;
}
    
