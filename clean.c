#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void clean2D(Domain *D)
{
    Particle **particle;
    particle=D->particle;
    LoadList *LL,*tmpLL;
    LaserList *L, *tmpL;
    int i,j,istart,iend,jstart,jend,s;
    ptclList *p,*tmp;
 
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    //remove particles
    for(i=0; i<=iend; i++)
      for(j=jstart-1; j<=jend; j++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;
          while(p)
          {	
            tmp=p->next;
            particle[i][j].head[s]->pt=tmp; 
            p->next=NULL;
            free(p);
            p=particle[i][j].head[s]->pt;
          }
          free(particle[i][j].head[s]);
        }
        free(particle[i][j].head);
      }		
      
    LL=D->loadList;
    while(LL)
    {	
      tmpLL=LL->next;
      D->loadList=tmpLL; 
      LL->next=NULL;
      free(LL);
      LL=D->loadList;
    }
    free(D->loadList);

    L=D->laserList;
    while(L)
    {	
      tmpL=L->next;
      D->laserList=tmpL; 
      L->next=NULL;
      free(L);
      L=D->laserList;
    }
    free(D->laserList);

    //remove trans
    free(D->btJ);
    free(D->upJ);


    //remove field
    if(D->fieldType==1)
    {
      FieldDSX **field;
      field=D->fieldDSX;

      for(i=0; i<D->nxSub+5; i++)
        free(field[i]);
      free(field);
    }

    //remove probe
    if(D->probeNum>0)
    {
      for(i=0; i<D->probeNum; i++)
       free(D->probe[i]);
      free(D->probe);
      free(D->probeX);
      free(D->probeY);
    }

/*
    //remove boost field
    Boost **boost;
    boost=D->boost;

    for(i=0; i<D->nxSub+5; i++)
      free(boost[i]);
    free(boost);
*/
}

