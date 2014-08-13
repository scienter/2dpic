#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

void movingDomain2D(Domain *D)
{
    int i,j,istart,iend,jstart,jend,s;
    float x;
    ptclList *p,*tmp,*New;      
    Particle **particle;
    particle=D->particle;
   
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
/* 
    //delete first column   
    for(j=jstart-1; j<=jend; j++)
    {
      for(s=0; s<D->nSpecies; s++)
      {
        p=particle[istart][j].head[s]->pt;
        while(p)
        {
          tmp=p->next;
          particle[istart][j].head[s]->pt=tmp;
          p->next=NULL;
          free(p);
          p=particle[istart][j].head[s]->pt;
        }
      }
      
    }
 
*/
    if(D->fieldType==1)
    {
      FieldDSX **field;
      field=D->fieldDSX;

      for(i=istart-1; i<iend; i++)
        for(j=jstart-1; j<jend; j++)
        {
          field[i][j].E1=field[i+1][j].E1;
          field[i][j].B1=field[i+1][j].B1;
          field[i][j].Pr=field[i+1][j].Pr;
          field[i][j].Pl=field[i+1][j].Pl;
          field[i][j].Sr=field[i+1][j].Sr;
          field[i][j].Sl=field[i+1][j].Sl;
          field[i][j].E1C=field[i+1][j].E1C;
          field[i][j].B1C=field[i+1][j].B1C;
          field[i][j].PrC=field[i+1][j].PrC;
          field[i][j].PlC=field[i+1][j].PlC;
          field[i][j].SrC=field[i+1][j].SrC;
          field[i][j].SlC=field[i+1][j].SlC;
          field[i][j].J1Old=field[i+1][j].J1Old;
          field[i][j].J2Old=field[i+1][j].J2Old;
          field[i][j].J3Old=field[i+1][j].J3Old;
          field[i][j].J1=field[i+1][j].J1;
          field[i][j].J2=field[i+1][j].J2;
          field[i][j].J3=field[i+1][j].J3;
        }     
 
      for(i=istart; i<=iend; i++)
        for(j=jstart; j<jend; j++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            p=particle[i][j].head[s]->pt;
            while(p)
            {
              p->x-=1.0;  
              p->oldX-=1.0;
              p=p->next;
            }
          } 
        }	
    }          //End of fieldType=1

    D->minXSub+=1;
    D->maxXSub+=1;
}  
     
