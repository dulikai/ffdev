

#include<stdio.h>
#include "utils/mytype.h"
#include "tools/rdiotools.h"
#include "tools/tokenst.h"
#include "tools/opstr.h"

#include "readgjf.h"


#ifndef INIT_GRAPHIC
#define INIT_GRAPHIC

void createGraphic(structGraphic *g, char *file)
{
     return;
     }    

//#if 0
/*
 fill graphic g, with original value.
*/
 void init_graphic(structGraphic *g)
 {
   int i, j;
   
   for(i = 0; i < VertexNum; i++)
   {
      for(j = 0; j < VertexNum; j++)
        {
            g->edge[i][j] = INFINITY;
        }  

	  g->vex[i].id = -1;
	  strcpy(g->vex[i].ffname, "");
	  strcpy(g->vex[i].unique_name, "");
	  g->vex[i].inode = -1;
   }
  }

//#endif
#endif  /* INIT_GRAPHIC */


 int read_gjf_match_ff(char *gjffile, char *suffix, structGraphic *g)
 {
   FILE *fp;
   char matchfile[50];
   char line[200], name[10], file[50];
   char *p;
   int i;

   i = 0;

   strcpy(file, gjffile);
   p = strtok(file, ".");
   sprintf(matchfile, "%s_%s.txt", p, suffix);

   fp = fopen(matchfile, "r");
   if(fp == NULL)
   {
     printf("fatal error: match file of gjf not exist! \nfilename: %s\n", matchfile);
	 exit(1);
   }

   //JmpOrNot(fp, '#');
   //fscanf(fp, "%s", name);
   //strcpy(g->name, name);
   JmpOrNot(fp, '#');

   while(fgets(line, sizeof(line), fp) != NULL)
   {
      if(is_blank_line(line))
          break;
      else
      {
          sscanf(line, "%s", name);
         
          strcpy(g->vex[i].ffname, name);
          //printf("%s\n", name);
          i++;
      }                     
    }

   return i;
 }


 /*
  read gjf Cartesian file.
 */
void read_gjf_cartesian(char *file, structGraphic *g)
{
   FILE *fp;
   int flag;
   int i, k,  nedge;
   int ncount, key[20];
   decimal value[20];
   char line[300], name[15];
   
   flag = 0;   i = k = 0;   nedge = 0;
   
   init_graphic(g);

   fp = fopen(file, "r");
     
   JmpOrNot(fp, '%');
   JmpOrNot(fp, '#');
   JmpSpace(fp);
   toss(fp, 3);
   
   while(fgets(line, sizeof(line), fp) != NULL)
   {
      if(is_blank_line(line))
          break;
      else
      {
          sscanf(line, "%s", name);
          g->vex[i].id = i+1;
         
          strcpy(g->vex[i].element_name, name);
          //printf("%s\n", name);
          i++;
      }                     
    }
    
    g->n = i;
    
    while(fgets(line, sizeof(line), fp) != NULL)
    {
         if(is_blank_line(line))
             break;
                      
         ncount = count_tokens(ltrim(line), " ");  
		 //printf("%d\n\n", ncount);
         nedge += ncount/2; 
		 k = 0;
         for(i = 0; i < ncount; i++)
         {
             if(i == 0 || i%2 == 1)             
                 key[k++] = atoi(gettoken(line," ",i+1)); 
             else
                 value[0] = (decimal)atof(gettoken(line," ",i+1));              
         }                
          
         //g->edge[key[0]-1][key[0]-1] = INFINITY;
		 //printf("%d\n", k);
         for(i = 1; i < k; i++)
         {
            g->edge[key[0]-1][key[i]-1] = 1;
            g->edge[key[i]-1][key[0]-1] = 1;
			
         }  
                
     }
     
     g->e = nedge;
	 //printf("%d\n", nedge);

    if(read_gjf_match_ff(file, "match", g) != g->n)
	{
	  printf("ERROR OCCURED... FILE:%s LINE:%s\n", __FILE__, __LINE__);
	  exit(1);
	}

  return;

 }


/*
  print the element in graphic g...
*/
 void print_graphic(char *file, structGraphic *g)
 {
   int n, e;
   int i,j;
   FILE *fp;

   fp = fopen(file, "w");

   n = g->n;
   e = g->e;

   fprintf(fp, "number of atom : %d, number of edge: %d\n", n, e);
   fprintf(fp, "\n[the NAME OF the SYSTEM: ]\n");

   for(i = 0; i < n; i++)
   {
     fprintf(fp, "id|inode: %3d|%3d  element name: %s  FF type name: %s\n", 
		 g->vex[i].id, g->vex[i].inode, g->vex[i].element_name, g->vex[i].ffname);
   }

   for(i = 0; i < n; i++)
   {
	   fprintf(fp, "\n%d ", i+1);
	   for(j = 0; j < n; j++)
	   {
		   if(g->edge[i][j] < INFINITY && i < j)
		      fprintf(fp, "%d ", j+1);
	   }
   }

   fprintf(fp, "\n\nMAP\n    ");
   for(i = 0; i < n; i++)
	   fprintf(fp, "%5d ", i+1);

   for(i = 0; i < n; i++)
   {
	   fprintf(fp, "\n%3d ", i+1);
	   for(j = 0; j < n; j++)
	   {
		   if(g->edge[i][j] < INFINITY)
		      fprintf(fp, "%5d ", 1);
		   else
		      fprintf(fp, "%5d ", 0);
	   }
   }

   fclose(fp);
 }

