
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "utils/mytype.h"
#include "tools/rdiotools.h"
#include "tools/opstr.h"
#include "readff.h"


 int gen_ff_inode(FFP *f, char *ffname)
 {
   int i, n, inode;

   inode = -1; 
   n = f->nndx;

   for(i = 0; i < n; i++)
   {
     if(!strcmp(f->ndx[i].ffname, ffname))
	 {
       inode = f->ndx[i].inode;
	   break;
	 }
   }

   return inode;
 }

 
 FFP *init_ffp()
 {
   FFP *f;

   f = (FFP*)malloc(sizeof(FFP));
/*
   f->ndx = (ffndx*)malloc(50*sizeof(ffndx));
   f->atomtype = (ffatomtype*)malloc(50*sizeof(ffatomtype));
   f->alias = (ffalias*)malloc(10*sizeof(ffalias));
   f->bond = (ffbond*)malloc(100*sizeof(ffbond));
   f->angle = (ffangle *)malloc(100*sizeof(ffangle));
   f->dihedral = (ffdihedral *)malloc(200*sizeof(ffdihedral));
   f->improper = (ffimproper *)malloc(50*sizeof(ffimproper));
*/
//可以将所有量初始化为零。。。，嫌麻烦，没有做

   return f;
 }

 void read_ff_db_atoms(FILE *fp, FFP *f)
 {
    int iatom;
    decimal mass, charge, sig, eps;
    char ffname[10], field[20];
	char str[50], title[50], line[200];

	iatom = 0;

    fgets(str, 50, fp);
	sscanf(str, "%s", title);
	//printf("'%s'\n", title);
    if(strcmp(title, "ATOMS"))
    {
      printf("illegal section ATOM...\n");
	  exit(0);
    }
   JmpOrNot(fp, '#');

   /*
   read ATOMS section...
   */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
	  printf("%s", line);
      if(is_blank_line(line))
          break;
	  else if(strchr(line, '#'))
		   continue;
      else
      {
          sscanf(line, "%s%lf%lf%s%lf%lf", ffname, &mass, &charge, field, &sig, &eps);
          
		  f->ndx[iatom].inode = iatom;
          strcpy(f->ndx[iatom].ffname, ffname);

          f->atomtype[iatom].inode = iatom;
		  f->atomtype[iatom].mass = mass;
		  f->atomtype[iatom].charge = charge;
		  f->atomtype[iatom].sig = sig;
		  f->atomtype[iatom].eps = eps;
		  strcpy(f->atomtype[iatom].field, field);         
          //printf("%s\n", name);
          iatom++;
      }   
	  
    }
    f->natomtype = iatom;
	f->nndx = iatom;

	return ;
 
 }




 void read_ff_db_translation(FILE *fp, FFP *f)
 {
    int itran;
	char dnst[10], source[10];
	char str[50], title[50], line[200];
	itran = 0;

    fgets(str, 50, fp);
	sscanf(str, "%s", title);
    if(strcmp(title, "TRANSLATION"))
    {
      printf("illegal section TRANSLATION...\n");
	  exit(0);
    } 
   
   /*
   read TRANSLATION section...
   */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
      if(is_blank_line(line))
          break;
	  else if(strchr(line, '#'))
		   continue;
      else
      {
          sscanf(line, "%s%s", dnst, source);
          
		  f->alias[itran].id = itran;
          f->alias[itran].idnst = gen_ff_inode(f, dnst);
          f->alias[itran].isource = gen_ff_inode(f, source);

          //printf("%s\n", name);
          itran++;
      }                     
    }
    f->nalias = itran;
	
	return;
 }


  void read_ff_db_bonds(FILE *fp, FFP *f)
  {
     int ibond;
	 char str[50], title[50], line[200];
	 char atom_a[10], atom_b[10], type[10];
	 decimal re, ka;
	 int inode, jnode;

	 ibond = 0;

     fgets(str, 50, fp);
	 sscanf(str, "%s", title);
     if(strcmp(title, "BONDS"))
	 {
       printf("illegal section BONDS...\n");
	   exit(1);
	 } 

     /*
      read BONDS section...
    */
     while(fgets(line, sizeof(line), fp) != NULL)
	 {
        if(is_blank_line(line))
           break;
	    else if(strchr(line, '#'))
		   continue;
        else
		{
           sscanf(line, "%s%s%s%lf%lf", atom_a, atom_b, type, &re, &ka);
		   inode = gen_ff_inode(f, atom_a);
		   jnode = gen_ff_inode(f, atom_b);

		   f->bond[ibond].id = ibond;
		   f->bond[ibond].inode = inode;
		   f->bond[ibond].jnode = jnode;
		   f->bond[ibond].re = re;
		   f->bond[ibond].ka = ka;

		   strcpy(f->bond[ibond].type, type);          

           //printf("%s\n", name);
           ibond++;
		}                     
	 }

     f->nbond = ibond;

	 return ;
  }

 void read_ff_db_angles(FILE *fp, FFP *f)
 {
   int iangle;
   int inode, jnode, knode;
   char str[50], title[50], line[200];
   char atom_a[10], atom_b[10], atom_c[10], type[10];
   decimal th, ka;

   iangle = 0;

   fgets(str, 50, fp);
   sscanf(str, "%s", title);
   if(strcmp(title, "ANGLES"))
   {
       printf("illegal section ANGLES...\n");
	   exit(1);
   } 

     /*
      read ANGLES section...
    */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
      if(is_blank_line(line))
         break;
	  else if(strchr(line, '#'))
		   continue;
      else
	  {
        sscanf(line, "%s%s%s%s%lf%lf", atom_a, atom_b, atom_c, type, &th, &ka);
		inode = gen_ff_inode(f, atom_a);
		jnode = gen_ff_inode(f, atom_b);
		knode = gen_ff_inode(f, atom_c);

		f->angle[iangle].id = iangle;

		f->angle[iangle].inode = inode;
		f->angle[iangle].jnode = jnode;
		f->angle[iangle].knode = knode;

		f->angle[iangle].th = th;
		f->angle[iangle].ka = ka;

		strcpy(f->angle[iangle].type, type);          

           //printf("%s\n", name);
        iangle++;
	  }                     
   }

   f->nangle = iangle;

    return ;

 }



 void read_ff_db_dihedral(FILE *fp, FFP *f)
 {
   int idihedral;
   int inode, jnode, knode, lnode;
   char str[50], title[50], line[200];
   char atom_a[10], atom_b[10], atom_c[10], atom_d[10], type[10];
   decimal v1, v2, v3, v4;

   idihedral = 0;

   fgets(str, 50, fp);
   sscanf(str, "%s", title);
   if(strcmp(title, "DIHEDRALS"))
   {
       printf("illegal section DIHEDRALS...\n");
	   exit(1);
   } 

     /*
      read ANGLES section...
    */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
      if(is_blank_line(line))
         break;
	  else if(strchr(line, '#'))
		   continue;
      else
	  {
        sscanf(line, "%s%s%s%s%s%lf%lf%lf%lf", atom_a, atom_b, atom_c, atom_d, type, &v1, &v2, &v3, &v4);
		inode = gen_ff_inode(f, atom_a);
		jnode = gen_ff_inode(f, atom_b);
		knode = gen_ff_inode(f, atom_c);
		lnode = gen_ff_inode(f, atom_d);

		f->dihedral[idihedral].id = idihedral;

		f->dihedral[idihedral].inode = inode;
		f->dihedral[idihedral].jnode = jnode;
		f->dihedral[idihedral].knode = knode;
		f->dihedral[idihedral].lnode = lnode;


		f->dihedral[idihedral].v1 = v1;
		f->dihedral[idihedral].v2 = v2;
		f->dihedral[idihedral].v3 = v3;
		f->dihedral[idihedral].v4 = v4;

		strcpy(f->dihedral[idihedral].type, type);          

           //printf("%s\n", name);
        idihedral++;
	  }                     
   }

   f->ndihedral = idihedral;

    return ;
 }


 void read_ff_db_improper(FILE *fp, FFP *f)
 {
   int iimproper;
   int inode, jnode, knode, lnode;
   char str[50], title[50], line[200];
   char atom_a[10], atom_b[10], atom_c[10], atom_d[10], type[10];
   decimal v1, v2, v3, v4;

   iimproper = 0;

   fgets(str, 50, fp);
   sscanf(str, "%s", title);
   if(strcmp(title, "IMPROPERS"))
   {
       printf("illegal section IMPROPERS DIHEDRALS..\n");
	   exit(1);
   } 

     /*
      read ANGLES section...
    */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
      if(is_blank_line(line))
         break;
	  else if(strchr(line, '#'))
		   continue;
      else
	  {
        sscanf(line, "%s%s%s%s%s%lf%lf%lf%lf", atom_a, atom_b, atom_c, atom_d, type, &v1, &v2, &v3, &v4);
		inode = gen_ff_inode(f, atom_a);
		jnode = gen_ff_inode(f, atom_b);
		knode = gen_ff_inode(f, atom_c);
		lnode = gen_ff_inode(f, atom_d);

		f->improper[iimproper].id = iimproper;

		f->improper[iimproper].inode = inode;
		f->improper[iimproper].jnode = jnode;
		f->improper[iimproper].knode = knode;
		f->improper[iimproper].lnode = lnode;


		f->improper[iimproper].v1 = v1;
		f->improper[iimproper].v2 = v2;
		f->improper[iimproper].v3 = v3;
		f->improper[iimproper].v4 = v4;

		strcpy(f->improper[iimproper].type, type);          

           //printf("%s\n", name);
        iimproper++;
	  }
   }                     

   f->nimproper = iimproper;

   return ;
 }




 void read_ff_db(char *file, FFP *f)
 {
   FILE *fp;

   fp = connect_ff_db(file);

   JmpOrNot(fp, '#'); //jump over all line contain # or space character.
   //JmpSpace(fp);

   /* read atoms... */
   read_ff_db_atoms(fp, f);

   JmpOrNot(fp, '#');
   JmpSpace(fp);

   /* read translation ... */
   read_ff_db_translation(fp, f);

   JmpOrNot(fp, '#');

   /* read bonds ... */
   read_ff_db_bonds(fp, f);

   JmpOrNot(fp, '#');
   /* read angles... */
   read_ff_db_angles(fp, f);

   JmpOrNot(fp, '#');
   /* read dihedral... */
   read_ff_db_dihedral(fp, f);

   JmpOrNot(fp, '#');
   /* read improper... */
   read_ff_db_improper(fp, f);

   close_ff_db(fp);

   return;
 }


 void print_ff_db(char *file, FFP *f)
 {
   int i;
   FILE *fp;

   fp = fopen(file, "w");

   fprintf(fp, "variable in index table\n");
   for(i = 0; i < f->nndx; i++)
   {
     fprintf(fp, "%d:%s\n", f->ndx[i].inode, f->ndx[i].ffname);
   }
   
   fprintf(fp, "\nvariable of atomtype\n");
   for(i = 0; i < f->natomtype; i++)
   {
     fprintf(fp, "%d: mass=%.6lf, charge=%.6lf, eps=%.6lf, sig=%.6lf\n", f->atomtype[i].inode, 
		 f->atomtype[i].charge, f->atomtype[i].eps, f->atomtype[i].sig, f->atomtype[i].mass);
   }

   fprintf(fp, "\nvariable of bonds\n");
   for(i = 0; i < f->nbond; i++)
   {
     fprintf(fp, "%d: %d %d  re = %lf, ka=%lf\n", f->bond[i].id, f->bond[i].inode, 
		 f->bond[i].jnode, f->bond[i].ka, f->bond[i].re);
   }

   fprintf(fp, "\nvariable of angles\n");
   for(i = 0; i < f->nangle; i++)
   {
     fprintf(fp, "%d: %d %d %d  th = %lf, ka=%lf\n", f->angle[i].id, f->angle[i].inode, 
		 f->angle[i].jnode, f->angle[i].knode, f->angle[i].th, f->angle[i].ka);
   }

   fprintf(fp, "\nvariable of dihedrals\n");
   for(i = 0; i < f->ndihedral; i++)
   {
     fprintf(fp, "%d: %d %d %d %d v1 = %lf, v2=%lf, v3=%lf, v4=%lf\n", f->dihedral[i].id, f->dihedral[i].inode, 
		 f->dihedral[i].jnode, f->dihedral[i].knode, f->dihedral[i].lnode, 
		 f->dihedral[i].v1, f->dihedral[i].v2, f->dihedral[i].v3, f->dihedral[i].v4);
   }

   fprintf(fp, "\nvariable of impropers\n");
   for(i = 0; i < f->nimproper; i++)
   {
     fprintf(fp, "%d: %d %d %d %d v1 = %lf, v2=%lf, v3=%lf, v4=%lf\n", f->improper[i].id, f->improper[i].inode, 
		 f->improper[i].jnode, f->improper[i].knode, f->improper[i].lnode, 
		 f->improper[i].v1, f->improper[i].v2, f->improper[i].v3, f->improper[i].v4);
   }

   fclose(fp);
 }


 //未完成工作，对ffp结构的初始化，将未用值解析为0



