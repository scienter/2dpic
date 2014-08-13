#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void rearrangeParticles2D(Domain *D)
{
    Particle **particle;
    particle=D->particle;

    int i,j,s,intX,intY,cnt,deleteFlag=0;
    int istart,iend,jstart,jend;
    float x,y;
    ptclList *p,*New,*prev,*tmp;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    if(D->fieldType==1)
    {
      for(i=istart; i<=iend; i++)
        for(j=jstart; j<jend; j++)
          for(s=0; s<D->nSpecies; s++)
          {
            cnt=1;
            p=particle[i][j].head[s]->pt;
            while(p)
            {
              x=p->x;   y=p->y; 
              if(x>=1.0)  {
                intX=(int)x;
                x-=intX;
                deleteFlag=1;
              }
              else if(x<0) {              
                intX=(int)(x-1);
                x-=intX;
                deleteFlag=1;
              } 
              else   intX=0;
              if(y>=1.0)  {
                intY=(int)y;
                y-=intY;
                deleteFlag=1;
              }
              else if(y<0) {
                intY=(int)(y-1);
                y-=intY;
                deleteFlag=1;
              } 
              else   intY=0;            
  
              if(deleteFlag==1)
              {
                New = (ptclList *)malloc(sizeof(ptclList)); 
                New->next = particle[i+intX][j+intY].head[s]->pt;
                particle[i+intX][j+intY].head[s]->pt = New;
                New->x=x;    New->oldX=p->oldX;
                New->y=y;    New->oldY=p->oldY;
                New->p1=p->p1;  New->p2=p->p2;  New->p3=p->p3;
                New->index=p->index;
  
                if(cnt==1)
                {
                  tmp=p->next;
                  particle[i][j].head[s]->pt=tmp; 
                  p->next=NULL;
                  free(p);
                  p=particle[i][j].head[s]->pt; 
                  cnt=1;
                }
                else
                {
                  prev->next=p->next;
                  p->next=NULL;
                  free(p);
                  p=prev->next; 
                  cnt++;
                }
              }		//End of if(deleteFlag==1)
              else
              {
                prev=p;
                p=p->next;
                cnt++;
              }

              deleteFlag=0;
            }
          }		//End of for(s)
    }          //End of fieldType=1
}

