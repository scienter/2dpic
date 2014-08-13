#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void boostShot(Domain *D,int iteration,int labSaveStep)
{
    int i,j,labStep,istart,iend,jstart,jend;
    char name[100];
    float x,y,e1,e2,e3,b1,b2,b3;
    float xx,yy,E1,E2,E3,B1,B2,B3;
    float tmp,factor,factor1;
    int myrank, nprocs;
    Boost **boost;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    if(D->fieldType==1)
    {   
      FieldDSX **field;
      field=D->fieldDSX;
      boost=D->boost;

      factor=D->gamma*(1+D->beta);
      factor1=1+D->beta;
      i=(int)(labSaveStep/D->gamma/D->beta/factor-iteration/D->beta-D->minXSub+istart);
      if(i>=istart && i<iend)
        for(j=jstart; j<jend; j++) 
        {       
           //boost frame data
           x=(i-istart+D->minXSub)*D->dx*D->lambda;
           y=(j-jstart+D->minYSub)*D->dy*D->lambda;
           e1=field[i][j].E1;
           e2=field[i][j].Pr+field[i][j].Pl;
           e3=field[i][j].Sr+field[i][j].Sl;
           b1=field[i][j].B1;
           b2=field[i][j].Sl-field[i][j].Sr;
           b3=field[i][j].Pr-field[i][j].Pl;

           //lab frame data
           xx=D->gamma*(x+D->beta*iteration*D->dt*D->lambda);
           yy=y;
           E1=e1/factor;
           E2=(e2+D->beta*b3)/factor1;
           E3=(e3-D->beta*b2)/factor1;
           B1=b1/factor;
           B2=(b2-D->beta*b3)/factor1;
           B3=(b3+D->beta*e2)/factor1;

           boost[i][j].x=xx;   
           boost[i][j].y=yy;   
           boost[i][j].E1=E1;
           boost[i][j].B1=B1;
           boost[i][j].Pr=0.5*(E2+B3);
           boost[i][j].Pl=0.5*(E2-B3);
           boost[i][j].Sr=0.5*(E3-B2);
           boost[i][j].Sl=0.5*(E3+B2);
        }	//End of for(j)

    }			//End of fieldType=1
}

