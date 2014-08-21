#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>


void prePML(Domain *D)
{
    int i,j,istart,iend,jstart,jend,nxSub,nySub;  
    float dx,dy,dt,preB1C;
    float nowPr,nowSr,prevPr,prevSr;
    float nowPrC,nowSrC,prevPrC,prevSrC;
    FieldDSX **field;
    field=D->fieldDSX;
    UpPML **upPML;
    upPML=D->UpPML;
    DnPML **dnPML;
    dnPML=D->DnPML;

    dx=D->dx;
    dy=D->dy;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    //up
    for(i=istart; i<iend; i++)
      for(j=1; j<10; j++)
      {
        upPML[i][j].Dx=field[i][jend-1].PrC-field[i][jend-1].PlC;
        upPML[i][0].Hx=0.5*(field[i][jend].B1C+field[i][jend-1].B1C);
      }

}


void absorpbing(Domain *D)
{
    int i,j,istart,iend,jstart,jend,nxSub,nySub;  
    float dx,dy,dt,preB1C;
    float nowPr,nowSr,prevPr,prevSr;
    float nowPrC,nowSrC,prevPrC,prevSrC;
    float pmlValue,pmlCell,s,def;
    int myrank, nTasks;
    double sigma();

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    FieldDSX **field;
    field=D->fieldDSX;
    UpPML **upPML;
    upPML=D->UpPML;
    DnPML **dnPML;
    dnPML=D->DnPML;

    dx=D->dx;
    dy=D->dy;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    pmlCell=20;
    pmlValue=0.5;
    //up
    if(myrank==nTasks)
    {
      def=jend-pmlCell;
      for(i=istart; i<iend; i++)
        for(j=jend-1-pmlCell; j<jend; j++)
        {
          s=pmlValue*(j-def-0.5)*(j-def-0.5);
          field[i][j].E1 -= s*field[i][j].E1;

          s=pmlValue*(j-def)*(j-def);
          field[i][j].Pr -= s*field[i][j].Pr;
          field[i][j].Pl -= s*field[i][j].Pl;
        }
    }

}

void absorpbingC(Domain *D)
{
    int i,j,istart,iend,jstart,jend,nxSub,nySub;  
    float dx,dy,dt,preB1C;
    float nowPr,nowSr,prevPr,prevSr;
    float nowPrC,nowSrC,prevPrC,prevSrC;
    float pmlValue,pmlCell,s,def;
    int myrank, nTasks;
    double sigma();

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    FieldDSX **field;
    field=D->fieldDSX;
    UpPML **upPML;
    upPML=D->UpPML;
    DnPML **dnPML;
    dnPML=D->DnPML;

    dx=D->dx;
    dy=D->dy;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    pmlCell=20;
    pmlValue=0.5;
    //up
    if(myrank==nTasks)
    {
      def=jend-pmlCell;
      for(i=istart; i<iend; i++)
        for(j=jend-1-pmlCell; j<jend; j++)
        {
          s=pmlValue*(j-def-0.5)*(j-def-0.5);
          field[i][j].E1C -= s*field[i][j].E1C;

          s=pmlValue*(j-def)*(j-def);
          field[i][j].PrC -= s*field[i][j].PrC;
          field[i][j].PlC -= s*field[i][j].PlC;
        }
    }

}
double sigma(double pmlValue, double i, double pmlCell)
{
   float r,tmp;
 
   tmp=i/pmlCell;
   r = pmlValue*tmp*tmp;

   return r;
}
