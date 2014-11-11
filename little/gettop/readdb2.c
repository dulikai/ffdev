

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools/rdiotools.h"
#include "tools/opstr.h"
#include "tools/tokenst.h"
#include "conndb.h"
#include "readdb2.h"

#define For(i, n) for(i=0;i<n;i++){
#define End }




 int gen_aaff_inode(AAFFP *aa, char *ffname)
 {
   int i, n, inode;

   inode = -1; 
   n = aa->nndx;

   for(i = 0; i < n; i++)
   {
     if(!strcmp(aa->ndx[i].ffname, ffname))
	 {
       inode = aa->ndx[i].inode;
	   break;
	 }
   }

   return inode;
 }

 void read_aaff_db_aaname(FILE *fp, AAFFP *aa)
 {
       //AAseries name[AA_NUM];
    int i, nname, ncount;
	char str[50], title[50], line[200];

	nname = 0;
    fgets(str, 50, fp);
	sscanf(str, "%s", title);
	//printf("'%s'\n", title);
    if(strcmp(title, "ENTITY"))
    {
      printf("illegal section ENTITY...\n");
	  exit(1);
    }
   JmpOrNot(fp, '#');

   /*
   read ENTITY section...
   */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
	  //printf("%s", line);
      if(is_blank_line(line))
          break;
	  else if(strchr(line, '#'))
		   continue;
      else
      {
		  ncount = count_tokens(ltrim(line), " "); 
		  for(i = 0; i < ncount; i++)
		  {
		    strcpy(aa->name[i].aaname, gettoken(line," ",i+1));
			aa->name[i].id = nname;
			nname++;
		  }
      }   
	  
    }
    //aa->natomtype = iatom;
	aa->nname = nname;
 
 }


 int is_aaff_db_entity(AAFFP *aa, char *name)
 {
   int flag = 0;
   int i;

   For(i, aa->nname)
	   if(!stricmp(aa->name[i].aaname, name))
	   {
		   flag = 1; break;
	   }
   End
   
   return flag;

 }

 void read_aaff_db_ndxs(FILE *fp, AAFFP *aa)
 {
    int iatom;
    decimal mass, sig, eps;
    char ffname[10] ;
	char source[100];
	char str[50], title[50], line[512];


	iatom = 0;

    fgets(str, 50, fp);
	sscanf(str, "%s", title);
	//printf("'%s'\n", title);
    if(strcmp(title, "ATOMS"))
    {
      printf("illegal section ATOM...\n");
	  exit(1);
    }
   JmpOrNot(fp, '#');

   /*
   read ATOMS section...
   */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
	  //printf("%s", line);
      if(is_blank_line(line))
          break;
	  else if(strchr(line, '#'))
		   continue;
      else
      {
		  //fgets(scanline, sizeof(scanline), fp);
          sscanf(line, "%s%lf%lf%lf%s", ffname, &sig, &eps, &mass, source);
		  printf("%4s %-.5lf %-.5lf %-.5lf %10s\n", ffname, sig, eps, mass, source);
          
		  aa->ndx[iatom].inode = iatom;
          strcpy(aa->ndx[iatom].ffname, ffname);

          strcpy(aa->atomcnst[iatom].wffname, ffname);
		  aa->atomcnst[iatom].id = iatom;
		  aa->atomcnst[iatom].mass = mass;
		  aa->atomcnst[iatom].sig = sig;
		  aa->atomcnst[iatom].eps = eps;
		  strcpy(aa->atomcnst[iatom].source, source);
          iatom++;
      }   
	  
    }
    aa->natomcnst = iatom;
	aa->nndx = iatom;

	return ; 
 }

 void read_aaff_db_atomvar(FILE *fp, AAFFP *aa)
 {
   char ffname[20];

    int iatom, iaa;
	int nblankline;
    decimal charge ;
	char title[50], line[512];


	iatom = 0;
	iaa = 0;
	nblankline = 0;
   JmpOrNot(fp, '#');

   /*
   read ATOMS charge section...
   */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
	  printf("%s", line);
      if(is_blank_line(line))
	  {
	     nblankline++;

		 if(nblankline == 1)
		 {
		    aa->atomvar[iaa].natomtype = iatom;
			iatom = 0;
			iaa++;
            continue;
		 }
	     else
		 {
		     break;
		 }

	  }
	  else if(strchr(line, '$'))
	  {
         sscanf(line, "$%s", title);
		 if(!is_aaff_db_entity(aa, title))
		 {
           printf("illegal section ATOM: illegal entity occur...\n");
	       exit(1);
		 }
		 strcpy(aa->atomvar[iaa].aaname, title);
	     aa->atomvar[iaa].id = iaa;	 
         nblankline = 0;
	  }
	  else if(strchr(line, '#'))
		   continue;
      else
      {
		  //fgets(scanline, sizeof(scanline), fp);
          sscanf(line, "%s%lf", ffname, &charge);
		  //printf("%4s %-.5lf \n", ffname,  charge);
          
		  aa->atomvar[iaa].atomtype[iatom].charge = charge;
		  aa->atomvar[iaa].atomtype[iatom].id = iatom;
		  aa->ndx[iatom].inode = iatom;
          strcpy(aa->atomvar[iaa].atomtype[iatom].ffname, ffname);

          iatom++;
      } 	  
    }

   aa->natomvar = iaa;

	return ; 

 }

 void read_aaff_db_atoms(FILE *fp, AAFFP *aa)
 {
    read_aaff_db_aaname(fp, aa);
    read_aaff_db_ndxs(fp, aa);
    read_aaff_db_atomvar(fp, aa);

	return;
 }


 void read_aaff_db_translation(FILE *fp, AAFFP *aa)
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
	   printf("%s", line);
      if(is_blank_line(line))
          break;
	  else if(strchr(line, '#'))
		   continue;
      else
      {
          sscanf(line, "%s%s", dnst, source);
          
		  aa->trans[itran].id = itran;
		  strcpy(aa->trans[itran].dnst, dnst);
		  strcpy(aa->trans[itran].source, source);

          //f->alias[itran].idnst = gen_ff_inode(f, dnst);
          //f->alias[itran].isource = gen_ff_inode(f, source);

          itran++;
      }                     
    }
    aa->ntrans = itran;
	
	return;
 }

 void read_aaff_db_bondcnst(FILE *fp, AAFFP *aa)
 {
     int ibondcnst;
	 char str[50], title[50], line[512];
	 char *atom_a, *atom_b, atom[20], source[100];
	 decimal kr;
	 int inode, jnode;

	 ibondcnst = 0;

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
	    printf("%s", line);
        if(is_blank_line(line))
           break;
	    else if(strchr(line, '#'))
		   continue;
        else
		{
           sscanf(line, "%s%lf%s", atom, &kr, source);
		   atom_a = strtok(atom, "-");
		   atom_b = strtok(NULL, "-");
		   inode = gen_aaff_inode(aa, atom_a);
		   jnode = gen_aaff_inode(aa, atom_b);

		   aa->atomcnst[ibondcnst].id = ibondcnst;
		   aa->bondcnst[ibondcnst].inode = inode;
		   aa->bondcnst[ibondcnst].jnode = jnode;
		   aa->bondcnst[ibondcnst].kr = kr;

		   strcpy(aa->bondcnst[ibondcnst].source, source);          

           //printf("%s\n", name);
           ibondcnst++;
		}                     
	 }
     aa->nbondcnst = ibondcnst;

	 return ;
 }

 void read_aaff_db_bondvar(FILE *fp, AAFFP *aa)
 {
    int ibond, iaa;
	int nblankline;
    decimal r0;
	int inode, jnode;
	char *atom_a, *atom_b, atom[20] ;
	char title[50], line[512];

	ibond = 0;
	iaa = 0;
	nblankline = 0;
    JmpOrNot(fp, '#');

   /*
   read ATOMS charge section...
   */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
	  printf("%s", line);
      if(is_blank_line(line))
	  {
	     nblankline++;

		 if(nblankline == 1)
		 {
		    aa->bondvar[iaa].nbond = ibond;
			ibond = 0;
			iaa++;
            continue;
		 }
	     else
		 {
		     break;
		 }

	  }
	  else if(strchr(line, '$'))
	  {
         sscanf(line, "$%s", title);
		 if(!is_aaff_db_entity(aa, title))
		 {
           printf("illegal section BONDS: illegal entity occur...\n");
	       exit(1);
		 }
		 strcpy(aa->bondvar[iaa].aaname, title);
	     aa->bondvar[iaa].id = iaa;	 
         nblankline = 0;
	  }
	  else if(strchr(line, '#'))
		   continue;
      else
      {
		  //fgets(scanline, sizeof(scanline), fp);
          sscanf(line, "%s%lf", atom, &r0);
		  atom_a = strtok(atom, "-");
		  atom_b = strtok(NULL, "-");
		  //printf("%4s %-.5lf \n", ffname,  charge);
          inode = gen_aaff_inode(aa, atom_a);
		  jnode = gen_aaff_inode(aa, atom_b);

		  aa->bondvar[iaa].bond[ibond].r0 = r0;
		  aa->bondvar[iaa].bond[ibond].id = ibond;
		  aa->bondvar[iaa].bond[ibond].inode = inode;
		  aa->bondvar[iaa].bond[ibond].jnode = jnode;
          strcpy(aa->bondvar[iaa].bond[ibond].catom, atom_a);
          strcpy(aa->bondvar[iaa].bond[ibond].datom, atom_b);

          ibond++;
      } 	  
    }
   aa->nbondvar = iaa;

	return ; 
 }


 void read_aaff_db_bonds(FILE *fp, AAFFP *aa)
 {
     read_aaff_db_bondcnst(fp, aa);
     read_aaff_db_bondvar(fp, aa); 

	 return;
 }

 void read_aaff_db_angleconst(FILE *fp, AAFFP *aa)
 {
   int ianglecnst;
   int inode, jnode, knode;
   char str[50], title[50], line[200];
   char *atom_a, *atom_b, *atom_c, atom[30] ;
   decimal kth ;

   ianglecnst = 0;

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
	  printf("%s", line);
      if(is_blank_line(line))
         break;
	  else if(strchr(line, '#'))
		   continue;
      else
	  {
        sscanf(line, "%s%lf", atom, &kth);
		atom_a = strtok(atom, "-");
		atom_b = strtok(NULL, "-");
		atom_c = strtok(NULL, "-");

		inode = gen_aaff_inode(aa, atom_a);
		jnode = gen_aaff_inode(aa, atom_b);
		knode = gen_aaff_inode(aa, atom_c);

		aa->anglecnst[ianglecnst].id = ianglecnst;
      //aa->anglecnst[ianglecnst]->id
		aa->anglecnst[ianglecnst].inode = inode;
		aa->anglecnst[ianglecnst].jnode = jnode;
		aa->anglecnst[ianglecnst].knode = knode;

		aa->anglecnst[ianglecnst].kth = kth;

		//strcpy(f->anglecnst[ianglecnst].type, type);          

           //printf("%s\n", name);
        ianglecnst++;
	  }                     
   }

   aa->nanglecnst = ianglecnst;

    return ;

 }

 void read_aaff_db_anglevar(FILE *fp, AAFFP *aa)
 {
    int iangle, iaa;
	int nblankline;
    decimal th;
	int inode, jnode, knode;
    char *atom_a, *atom_b, *atom_c, atom[30] ;
	char title[50], line[512];

	iangle = 0;
	iaa = 0;
	nblankline = 0;
    JmpOrNot(fp, '#');

   /*
   read ATOMS charge section...
   */
   while(fgets(line, sizeof(line), fp) != NULL)
   {
	  printf("%s", line);
      if(is_blank_line(line))
	  {
	     nblankline++;

		 if(nblankline == 1)
		 {
		    aa->anglevar[iaa].nangle = iangle;
			iangle = 0;
			iaa++;
            continue;
		 }
	     else
		 {
		     break;
		 }

	  }
	  else if(strchr(line, '$'))
	  {
         sscanf(line, "$%s", title);
		 if(!is_aaff_db_entity(aa, title))
		 {
           printf("illegal section angleS: illegal entity occur...\n");
	       exit(1);
		 }
		 strcpy(aa->anglevar[iaa].aaname, title);
	     aa->anglevar[iaa].id = iaa;	 
         nblankline = 0;
	  }
	  else if(strchr(line, '#'))
		   continue;
      else
      {
		  //fgets(scanline, sizeof(scanline), fp);
          sscanf(line, "%s%lf", atom, &th);

		  atom_a = strtok(atom, "-");
		  atom_b = strtok(NULL, "-");
		  atom_c = strtok(NULL, "-");
		  //printf("%4s %-.5lf \n", ffname,  charge);
          inode = gen_aaff_inode(aa, atom_a);
		  jnode = gen_aaff_inode(aa, atom_b);
		  knode = gen_aaff_inode(aa, atom_c);

		  aa->anglevar[iaa].angle[iangle].th = th;
		  aa->anglevar[iaa].angle[iangle].id = iangle;
		  aa->anglevar[iaa].angle[iangle].inode = inode;
		  aa->anglevar[iaa].angle[iangle].jnode = jnode;
		  aa->anglevar[iaa].angle[iangle].knode = knode;

          //strcpy(aa->anglevar[iaa].atomtype[iangle].ffname, ffname);

          iangle++;
      } 	  
    }
   aa->nanglevar = iaa;

	return ; 
 }


 void read_aaff_db_angles(FILE *fp, AAFFP *aa)
 {
   read_aaff_db_angleconst(fp, aa);
   read_aaff_db_anglevar(fp, aa);

   return;
 }

 void read_aaff_db_dihedral(FILE *fp, AAFFP *aa)
 {
   int idihedral;
   int inode, jnode, knode, lnode;
   int multi;
   char str[50], title[50], line[200];
   char *atom_a, *atom_b, *atom_c, *atom_d, atom[40];
   char source[100];
   decimal gam, kphi;

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
        printf("%s", line); 
	   if(is_blank_line(line))
         break;
	  else if(strchr(line, '#'))
		   continue;
      else
	  {
        sscanf(line, "%s%lf%lf%d%s", atom, &gam, &kphi, &multi, source);
	    atom_a = strtok(atom, "-");
	    atom_b = strtok(NULL, "-");
	    atom_c = strtok(NULL, "-");
	    atom_d = strtok(NULL, "-");

		inode = gen_aaff_inode(aa, atom_a);
		jnode = gen_aaff_inode(aa, atom_b);
		knode = gen_aaff_inode(aa, atom_c);
		lnode = gen_aaff_inode(aa, atom_d);

		aa->dihedral[idihedral].id = idihedral;

		aa->dihedral[idihedral].inode = inode;
		aa->dihedral[idihedral].jnode = jnode;
		aa->dihedral[idihedral].knode = knode;
		aa->dihedral[idihedral].lnode = lnode;


		aa->dihedral[idihedral].gam = gam;
		aa->dihedral[idihedral].kphi = kphi;
		aa->dihedral[idihedral].multi = multi;

		strcpy(aa->dihedral[idihedral].source, source);          

           //printf("%s\n", name);
        idihedral++;
	  }                     
   }

   aa->ndihedral = idihedral;

    return ;
 }


 void read_aaff_db(char *file, AAFFP *aa)
 {
    FILE *fp;

    fp = connect_ff_db(file);

    JmpOrNot(fp, '#'); //jump over all line contain # or space character.
    //JmpSpace(fp);

    /* read atoms... */
    read_aaff_db_atoms(fp, aa);

    JmpOrNot(fp, '#');

   /* read translation ... */
    read_aaff_db_translation(fp, aa);

    JmpOrNot(fp, '#');

   /* read bonds ... */
    read_aaff_db_bonds(fp, aa);

	goto_target_line(fp, '!', 1L);
    JmpOrNot(fp, '#');
   /* read angles... */
    read_aaff_db_angles(fp, aa);

   goto_target_line(fp, '!', 1L);
   JmpOrNot(fp, '#');
   /* read dihedral... */
   read_aaff_db_dihedral(fp, aa);

   //JmpOrNot(fp, '#');
   /* read improper... */
   //read_aaff_db_improper(fp, f);

   close_ff_db(fp);

   return; 
 }

#undef For
#undef End



 void print_aaffp(char *file, AAFFP *aa)
 {
   FILE *fp;
   int i,j;

   fp = fopen(file, "w");

   fprintf(fp, "\nindex table\n");

   for(i = 0; i < aa->nndx; i++)
   {
     fprintf(fp, "inode: %d --> ffname: %s\n", aa->ndx[i].inode, aa->ndx[i].ffname);
   }

   fprintf(fp, "\natomtype table : vdW parameter\n");
   for(i = 0; i < aa->natomcnst; i++)
   {
     fprintf(fp, "id:%3d mass:%10.6lf sig:%10.6lf eps:%10.6lf\n", aa->atomcnst[i].id, 
		 aa->atomcnst[i].mass,aa->atomcnst[i].sig, aa->atomcnst[i].eps);
   }

   fprintf(fp, "atomtype table: charge\n\n\n");
   for(i = 0; i < aa->natomvar; i++)
   {
      fprintf(fp, "\n[AA:%s]==> number: %3d\n", aa->atomvar[i].aaname, aa->atomvar[i].natomtype);
	  for(j = 0; j < aa->atomvar[i].natomtype; j++)
	  {
	    fprintf(fp, "%4s : charge %10.6lf\n", aa->atomvar[i].atomtype[j].ffname,
			aa->atomvar[i].atomtype[j].charge);
	  }
   }

   fprintf(fp, "bonds table: Kr\n\n");

   for(i = 0; i < aa->nbondcnst; i++)
   {
     fprintf(fp, "inode | jnode %2d|%2d  kr %10.6lf\n", aa->bondcnst[i].inode,
		 aa->bondcnst[i].jnode, aa->bondcnst[i].kr);
   }

   fprintf(fp, "bonds table: degree\n\n");
   for(i = 0; i < aa->nbondvar; i++)
   {
      fprintf(fp, "\n[AA:%s]==> number: %3d\n", aa->bondvar[i].aaname, aa->bondvar[i].nbond);
	  for(j = 0; j < aa->bondvar[i].nbond; j++)
	  {
	    fprintf(fp, "%4s-%4s[%2d %2d]: r0 %10.6lf\n", aa->bondvar[i].bond[j].catom,
		aa->bondvar[i].bond[j].datom, 
		aa->bondvar[i].bond[j].inode, aa->bondvar[i].bond[j].jnode,
		aa->bondvar[i].bond[j].r0);
	  }
   }


   fprintf(fp, "angles table: Kr\n\n");

   for(i = 0; i < aa->nanglecnst; i++)
   {
     fprintf(fp, "%2d %2d %2d  kth %10.6lf\n", 
		 aa->anglecnst[i].inode, aa->anglecnst[i].jnode,aa->anglecnst[i].knode, 
		 aa->anglecnst[i].kth);
   }

   fprintf(fp, "angles table: degree\n\n");
   for(i = 0; i < aa->nanglevar; i++)
   {
      fprintf(fp, "\n[AA:%s]==> number: %3d\n", aa->anglevar[i].aaname, aa->anglevar[i].nangle);
	  for(j = 0; j < aa->anglevar[i].nangle; j++)
	  {
	    fprintf(fp, "[%2d %2d %2d]: th %10.6lf\n", 
		aa->anglevar[i].angle[j].inode, aa->anglevar[i].angle[j].jnode, aa->anglevar[i].angle[j].knode,
		aa->anglevar[i].angle[j].th);
	  }
   }

   
   fprintf(fp, "dihedrals table: all\n\n");

   for(i = 0; i < aa->ndihedral; i++)
   {
     fprintf(fp, "%2d %2d %2d %2d gam %10.6lf  kphi %10.6lf multi %2d\n",
		 aa->dihedral[i].inode, aa->dihedral[i].jnode, 
		 aa->dihedral[i].inode, aa->dihedral[i].jnode,
		 aa->dihedral[i].gam, aa->dihedral[i].kphi,aa->dihedral[i].multi );
   }


   
   fclose(fp);

	 return;
 }


