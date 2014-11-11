
#include<stdio.h>
#include<string.h>

#include<stdlib.h>

#include "tools/rdiotools.h"

#include "fileinit.h"

 void init_filecontrol(char *file, FILEControl *fc)
 {
   FILE *fp;

   if(strlen(file) > 0)   
   {
       fp = fopen(file, "r");
   }
   else
   {
       fp = fopen("control.txt", "r");
       }
       
   if(fp == NULL)
   {
     printf("fatal error: match file of gjf not exist! \nfilename: %s\n", file);
	 exit(1);
   }
   
      JmpOrNot(fp, '#');
      fscanf(fp, "%s%s", fc->structfile, fc->structtype);
      
      JmpOrNot(fp, '#');
      fscanf(fp, "%s%s%s", fc->dbfile, fc->dbtype, fc->dbpara);      

      JmpOrNot(fp, '#');
      fscanf(fp, "%s", fc->templatefile);   
      
      JmpOrNot(fp, '#');
      fscanf(fp, "%s%s", fc->outitpfile, fc->itp_nobfile); 
      
      fclose(fp);    
           
      }        
