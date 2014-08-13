#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

void probe(Domain *D,int iteration)
{
    int i,j,istart,iend,jstart,jend,n,iter;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
   
    iter=iteration;

    if(D->fieldType==1)
    {
       FieldDSX **field;
       field=D->fieldDSX;

       for(n=0; n<D->probeNum; n++)
       {
         if(D->probeX[n]>=D->minXSub && D->probeX[n]<D->maxXSub && 
            D->probeY[n]>=D->minYSub && D->probeY[n]<D->maxYSub)
         {
            i=D->probeX[n]-D->minXSub+istart;
            j=D->probeY[n]-D->minYSub+jstart;
            D->probe[n][iter].E1=field[i][j].E1;    
            D->probe[n][iter].Pr=field[i][j].Pr;
            D->probe[n][iter].Pl=field[i][j].Pl;
            D->probe[n][iter].B1=field[i][j].B1;    
            D->probe[n][iter].Sr=field[i][j].Sr;
            D->probe[n][iter].Sl=field[i][j].Sl;      
         }
       }
    }

}

