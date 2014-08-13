#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"

void solveField2DC_DSX(Domain *D)
{
    int i,j,istart,iend,jstart,jend,nxSub,nySub;  
    float dx,dy,dt,preB1C;
    float nowPr,nowSr,prevPr,prevSr;
    float nowPrC,nowSrC,prevPrC,prevSrC;
    FieldDSX **field;
    field=D->fieldDSX;

    dx=D->dx;
    dy=D->dy;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    // PrC,PlC,E1C,SrC,SlC,B1C
    for(j=jstart; j<jend; j++)
    {
      nowPrC=field[0][j].PrC;
      nowSrC=field[0][j].SrC;
      for(i=istart; i<iend; i++)
      {
         field[i][j].E1C+=0.5*dt/dy*(field[i-1][j].Pr-field[i-1][j].Pl+field[i][j].Pr-field[i][j].Pl-field[i-1][j-1].Pr+field[i-1][j-1].Pl-field[i][j-1].Pr+field[i][j-1].Pl)-0.5*pi*dt*(field[i-1][j].J1Old+field[i-1][j].J1+field[i][j].J1Old+field[i][j].J1);
         field[i][j].B1C+=-0.5*dt/dy*(field[i-1][j].Sr+field[i-1][j].Sl+field[i][j].Sr+field[i][j].Sl-field[i-1][j-1].Sr-field[i-1][j-1].Sl-field[i][j-1].Sr-field[i][j-1].Sl);

         prevSrC=nowSrC;
         nowSrC=field[i][j].SrC;
         field[i][j].SrC=prevSrC-0.5*dt/dy*(field[i][j+1].B1-field[i][j].B1)-0.5*pi*dt*(field[i][j].J3+field[i][j].J3Old+field[i][j+1].J3+field[i][j+1].J3);
         field[i-1][j].SlC=field[i][j].SlC-0.5*dt/dy*(field[i][j+1].B1-field[i][j].B1)-0.5*pi*dt*(field[i][j].J3+field[i][j].J3Old+field[i][j+1].J3+field[i][j+1].J3Old);     
         prevPrC=nowPrC;
         nowPrC=field[i][j].PrC;
         field[i][j].PrC=prevPrC+0.5*dt/dy*(field[i][j+1].E1-field[i][j].E1)-0.5*pi*dt*(field[i][j].J2Old+field[i][j].J2);
         field[i-1][j].PlC=field[i][j].PlC-0.5*dt/dy*(field[i][j+1].E1-field[i][j].E1)-0.5*pi*dt*(field[i][j].J2Old+field[i][j].J2);
      }
    }
}

void solveField2D_DSX(Domain *D)
{
    int i,j,istart,iend,jstart,jend,nxSub,nySub;  
    float dx,dy,dt,preB1C;
    float nowPr,nowSr,prevPr,prevSr;
    float nowPrC,nowSrC,prevPrC,prevSrC;
    FieldDSX **field;
    field=D->fieldDSX;

    dx=D->dx;
    dy=D->dy;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    // Pr,Pl,E1,Sr,Sl,B1
    for(j=jstart; j<jend; j++)
    {
      nowPr=field[0][j].Pr;
      nowSr=field[0][j].Sr;
      for(i=istart; i<iend; i++)
      {
         field[i][j].E1+=0.5*dt/dy*(field[i-1][j].PrC-field[i-1][j].PlC+field[i][j].PrC-field[i][j].PlC-field[i-1][j-1].PrC+field[i-1][j-1].PlC-field[i][j-1].PrC+field[i][j-1].PlC)-pi*dt*(field[i-1][j].J1+field[i][j].J1);
         field[i][j].B1+=-0.5*dt/dy*(field[i-1][j].SrC+field[i-1][j].SlC+field[i][j].SrC+field[i][j].SlC-field[i-1][j-1].SrC-field[i-1][j-1].SlC-field[i][j-1].SrC-field[i][j-1].SlC);

         prevSr=nowSr;
         nowSr=field[i][j].Sr;
         field[i][j].Sr=prevSr-0.5*dt/dy*(field[i][j+1].B1C-field[i][j].B1C)-0.5*pi*dt*(field[i][j].J3+field[i][j+1].J3);
         field[i-1][j].Sl=field[i][j].Sl-0.5*dt/dy*(field[i][j+1].B1C-field[i][j].B1C)-0.5*pi*dt*(field[i][j].J3+field[i][j+1].J3);     
         prevPr=nowPr;
         nowPr=field[i][j].Pr;
         field[i][j].Pr=prevPr+0.5*dt/dy*(field[i][j+1].E1C-field[i][j].E1C)-pi*dt*field[i][j].J2;
         field[i-1][j].Pl=field[i][j].Pl-0.5*dt/dy*(field[i][j+1].E1C-field[i][j].E1C)-pi*dt*field[i][j].J2;
      }
    }
}

