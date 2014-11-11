

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include "graphic.h"

     
void gen_com(structGraphic *g, analyzeGraphic *com)
{
     assign_vex_name(g);
     gen_vex(g, com);
     gen_bond(g, com);
     gen_angle(g, com);
     gen_dihedral(g, com);
     gen_pair(g, com);

     return;
}     

 void get_element_from_com(analyzeGraphic *com, int inode, char *name)
 {
   int i, n;

   n = com->nvex;

   for(i = 0; i < n; i++)
   {
     if(inode == com->vex[i].inode)
	 {
	   strcpy(name, com->vex[i].element_name);
	   break;
	 }
   }

   return ;
 }

 static void gen_vex(structGraphic *g, analyzeGraphic *com)
{
     int i, n;
     n = g->n;
     
     for(i = 0; i < n; i++)
     {
           com->vex[i].id = g->vex[i].id;
		   strcpy(com->vex[i].element_name, g->vex[i].element_name);
           strcpy(com->vex[i].ffname, g->vex[i].ffname);
           strcpy(com->vex[i].unique_name, g->vex[i].unique_name);
           com->vex[i].inode = g->vex[i].inode;
           }
     com->nvex = n;
     
     return;      
     }
     
     
     
 static void gen_bond(structGraphic *g, analyzeGraphic *com)    
{
     int i, j, n;
     int nbond = 0;
     
     n = g->n;
     for(i = 0; i < n; i++)
        for(j = i+1; j < n; j++)
           {
                if(g->edge[i][j] < INFINITY)
                   {
                      com->bond[nbond].i = g->vex[i].id;
                      com->bond[nbond].j = g->vex[j].id; 
                      nbond++;          
                      }
                }
      com->nbond = nbond;
      
      return;           
     }
     
 static void gen_angle(structGraphic *g, analyzeGraphic *com)   
{
     int i, j, k , n;
     int nangle = 0;

     n = g->n;
     for(i = 0; i < n; i++)
         for(j = 0; j < n; j++)
            for(k = i+1; k < n; k++)
              {
                  if(g->edge[i][j] < INFINITY && g->edge[j][k] < INFINITY)
                    {
                      com->angle[nangle].i = g->vex[i].id;
                      com->angle[nangle].j = g->vex[j].id;
                      com->angle[nangle].k = g->vex[k].id;
                      nangle++;                                    
				  }
			}
     com->nangle = nangle;
     
     return;     
 }  
     
     
 static void gen_dihedral(structGraphic *g, analyzeGraphic *com)    
{
     int i, j, k, l, n;
     int ndihedral = 0;
	 int list[4];
	 
     
     n = g->n;
     for(i = 0; i < n; i++)
         for(j = 0; j < n; j++)
            for(k = 0; k < n; k++)
               for(l = i+1; l < n; l++)
              {
				  list[0] = i; list[1] = j; list[2] = k; list[3] = l;
                  if(g->edge[i][j] < INFINITY && g->edge[j][k] < INFINITY && g->edge[k][l] < INFINITY && !is_self_iter(list, 4))
                    {
                      com->dihedral[ndihedral].i = g->vex[i].id;
                      com->dihedral[ndihedral].j = g->vex[j].id;
                      com->dihedral[ndihedral].k = g->vex[k].id;
                      com->dihedral[ndihedral].l = g->vex[l].id;
                      ndihedral++;                                    
				  }
			   }
     com->ndihedral = ndihedral;
     
     return;     
     
 } 

 static void gen_pair(structGraphic *g, analyzeGraphic *com)
{
     int i, u, v;
	 int n;
     int ndihedral, npair = 0;
     int dist[VertexNum];
     
     n = g->n;
     ndihedral = com->ndihedral;
     
     //floyd_shortpath();
     
     for(i = 0; i < ndihedral; i++)
     {
         u = com->dihedral[i].i;
         v = com->dihedral[i].l;
         
         Dijkstra_distance(g, u-1, dist);
		/* for(j = 0; j < n; j++)
		 {
		   printf("%d %d: %d\n", u, v, dist[j]);
		 }
*/
         if(dist[v-1] >= 3 && !is_redundance_pair(com, npair, u, v))
         {
            com->pair[npair].i = u;
            com->pair[npair].l = v;    
            npair++;    
		 } 
          /*else if(dist[v-1] == 3 && !is_redundance_pair(com, npair, u, v))      
          {
            com->pair[npair].i = u;
            com->pair[npair].l = v;    
            npair++;                   
		  } */
           else
           {
               ;
		   }   
	 }
           
           com->npair = npair;
           
           return;     
 }


int is_self_iter(int list[], int n)
{
  int i, j;
  int flag = 0;

  for(i = 0; i < n; i++)
	  for(j = i+1; j < n; j++)
	  {
         if(list[i] == list[j])
		 {
			 flag = 1;
			 break;
		 }
	  }

	  return flag;
}

 static int is_redundance_pair(analyzeGraphic *com, int npair, int u, int v)
{
    int i;
    int is_redundance = 0;
    
    for(i = 0; i < npair; i++)
    {
          if(com->pair[i].i == u && com->pair[i].l == v)
          {
              is_redundance = 1;              
          }
          
    }
          
     return is_redundance;     
 }


 /*
 返回值dist满足如下条件
 索引编号从1到n，0为空，以符合编号规则。
 */
 static void Dijkstra_distance(structGraphic *g, int v0, int *dist)
{
     int i, v, w;
     int n;
     int final[VertexNum];
	 int min;
     
     n = g->n;

	 //v0 = v0 -1;
     
     for(v = 0; v < n; v++)
     {
        final[v]= 0;
        dist[v] = g->edge[v0][v];
     }
     dist[v0] = 0;
     final[v0] = 1;
     
     for(i = 1; i < n; i++)
     {
        min = INFINITY;
        for(w = 0; w < n; w++)
        {
              if(!final[w])
                  if(dist[w] < min)
                  {
                      v = w;
                      min = dist[w];         
                  }                           
          } 
          final[v] = 1;
          
          for(w = 0; w < n; w++)
             if(!final[w] && (min+g->edge[v][w] < dist[w]))
             {
                dist[w] = min + g->edge[v][w];
             }  
       }
       
       return;      
     }


 static void assign_vex_name(structGraphic *g)
 {
   int i, j, k, n;;

   n = g->n;

   for(i = 0; i < n; i++)
   {
	   k = 0;
	   if(strlen(g->vex[i].unique_name) > 0)
	   {
		  continue;
	   }
	   for(j = i; j < n; j++)
	   {
          if(strcmp(g->vex[i].ffname, g->vex[j].ffname) == 0)
		  {
		    k++;
			sprintf(g->vex[j].unique_name, "%s%d", g->vex[i].ffname, k);
		  }
		  else
		  {
		     continue;
		  }
	   }
   }

   return;
 }


  int get_inode_form_id(analyzeGraphic *com, int id)
  {
    int i, n, inode;
	
	n = com->nvex;

	for(i = 0; i < n; i++)
	{
	  if(id == com->vex[i].id)
	  {
	    break;
	  }
	}
	inode = com->vex[i].inode;

	return inode;
  }


 void print_com(char *file, analyzeGraphic *com)
 {
	int i;
	FILE *fp;

	fp = fopen(file, "w");

	fprintf(fp, "vertex value\n");
	for(i = 0; i < com->nvex; i++)
	{
	   fprintf(fp, "id|inode: %3d|%3d  element name: %s  ffname: %s unique name: %s\n", 
		   com->vex[i].id, com->vex[i].inode, com->vex[i].element_name,
		   com->vex[i].ffname, com->vex[i].unique_name);
	}

	fprintf(fp, "bond value\n");
	for(i = 0; i < com->nbond; i++)
	{
	   fprintf(fp, "id %2d: %d %d\n", i+1,
		   com->bond[i].i, com->bond[i].j);
	}

	fprintf(fp, "angle value\n");
	for(i = 0; i < com->nangle; i++)
	{
	   fprintf(fp, "id %2d: %d %d %d\n", i+1,
		   com->angle[i].i, com->angle[i].j, com->angle[i].k);
	}

	fprintf(fp, "dihedral value\n");
	for(i = 0; i < com->ndihedral; i++)
	{
	   fprintf(fp, "id %2d: %d %d %d %d\n", i+1,
		   com->dihedral[i].i, com->dihedral[i].j, com->dihedral[i].k, com->dihedral[i].l);
	}

	fprintf(fp, "pair value\n");
	for(i = 0; i < com->npair; i++)
	{
	   fprintf(fp, "id %2d: %d %d\n", i+1,
		   com->pair[i].i, com->pair[i].l);
	}

	fclose(fp);


 return;
 }
