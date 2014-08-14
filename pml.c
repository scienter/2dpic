#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"

void pmlC_DSX(Domain *D)
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
    {
      upPML[i][0].PrC=field[i][jend-1].PrC;
      upPML[i][0].PlC=field[i][jend-1].PlC;
      upPML[i][0].E1C=field[i][jend-1].E1C;
      upPML[i][0].SrC=field[i][jend-1].SrC;
      upPML[i][0].SlC=field[i][jend-1].SlC;
      upPML[i][0].B1C=field[i][jend-1].B1C;
      upPML[i][0].Pr=field[i][jend-1].Pr;
      upPML[i][0].Pl=field[i][jend-1].Pl;
      upPML[i][0].E1=field[i][jend-1].E1;
      upPML[i][0].Sr=field[i][jend-1].Sr;
      upPML[i][0].Sl=field[i][jend-1].Sl;
      upPML[i][0].B1=field[i][jend-1].B1;
      upPML[i][0].Bzy=field[i][jend-1].Pr-field[i][jend-1].Pl;
    }

    for(j=1; j<10; j++)
    {
      nowPrC=upPML[0][j].PrC;
      nowSrC=upPML[0][j].SrC;
      for(i=istart; i<iend; i++)
      {
         upPML[i][j].E1C+=0.5*dt/dy*(upPML[i-1][j].Pr-upPML[i-1][j].Pl+upPML[i][j].Pr-upPML[i][j].Pl-upPML[i-1][j-1].Pr+upPML[i-1][j-1].Pl-upPML[i][j-1].Pr+upPML[i][j-1].Pl);
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

