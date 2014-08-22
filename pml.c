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
    float pmlValue,s;
    int myrank, nTasks,pmlCell,def;
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

    pmlCell=D->pmlCell;
    pmlValue=1.0/(float)pmlCell/(float)pmlCell;
//    pmlValue=1.0/(float)pmlCell/(float)pmlCell/(float)pmlCell;

    //up
    if(myrank==nTasks-1)
    {
      def=jend-pmlCell;
      for(i=istart-1; i<=iend; i++)
        for(j=def; j<=jend; j++)
        {
          s=pmlValue*(j-def-0.5)*(j-def-0.5);
//          s=pmlValue*(j-def-0.5)*(j-def-0.5)*(j-def-0.5);
          field[i][j].E1 -= s*field[i][j].E1;

          s=pmlValue*(j-def)*(j-def);
//          s=pmlValue*(j-def)*(j-def)*(j-def);
          field[i][j].Pr -= s*field[i][j].Pr;
          field[i][j].Pl -= s*field[i][j].Pl;
        }
    }

    //down
    if(myrank==0)
    {
      def=jstart-1+pmlCell;
      for(i=istart-1; i<=iend; i++)
        for(j=jstart-1; j<=def; j++)
        {
          s=pmlValue*(def-j-0.5)*(def-j-0.5);
//          s=pmlValue*(j-def-0.5)*(j-def-0.5)*(j-def-0.5);
          field[i][j].E1 -= s*field[i][j].E1;

          s=pmlValue*(def-j)*(def-j);
//          s=pmlValue*(j-def)*(j-def)*(j-def);
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
    float pmlValue,s;
    int myrank, nTasks,pmlCell,def;
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

    pmlCell=D->pmlCell;
    pmlValue=1.0/(float)pmlCell/(float)pmlCell;
    def=jend-pmlCell;
//    pmlValue=1.0/(float)pmlCell/(float)pmlCell/(float)pmlCell;
    //up
    if(myrank==nTasks-1)
    {
      for(i=istart-1; i<=iend; i++)
        for(j=def; j<=jend; j++)
        {
          s=pmlValue*(j-def-0.5)*(j-def-0.5);
//          field[i][j].E1C = (1.0-s)*field[i][j].E1C;
          field[i][j].E1C -= s*field[i][j].E1C;

          s=pmlValue*(j-def)*(j-def);
          field[i][j].PrC -= s*field[i][j].PrC;
          field[i][j].PlC -= s*field[i][j].PlC;
//          field[i][j].PrC = (1.0-s)*field[i][j].PrC;
//          field[i][j].PlC = (1.0-s)*field[i][j].PlC;
        }
    }

    //down
    if(myrank==0)
    {
      def=jstart-1+pmlCell;
      for(i=istart-1; i<=iend; i++)
        for(j=jstart-1; j<=def; j++)
        {
          s=pmlValue*(def-j-0.5)*(def-j-0.5);
          field[i][j].E1C -= s*field[i][j].E1C;
//          field[i][j].E1C = (1.0-s)*field[i][j].E1C;

          s=pmlValue*(def-j)*(def-j);
          field[i][j].PrC -= s*field[i][j].PrC;
          field[i][j].PlC -= s*field[i][j].PlC;
//          field[i][j].PrC = (1.0-s)*field[i][j].PrC;
//          field[i][j].PlC = (1.0-s)*field[i][j].PlC;
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
