#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"

void boundary(Domain *D,External *Ext)
{
     int i,j,s,rank,remain,tmp,sub,nxSub,nySub,numdataUp,numdataBt;
     float min,max;
     int myrank, nTasks;
     MPI_Status status;

     MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

     if(D->divDirection==1)
     {
        D->nxSub=D->nx/nTasks;
        sub=D->nxSub;
        remain=D->nx%nTasks;
        min=max=0;
        rank=0;
        while(rank<nTasks)
        {
           if(rank<remain)   tmp=sub+1;
           else              tmp=sub;
           min=max;
           max=min+tmp;
           if(myrank==rank)
           {
              D->minXSub=min+D->minXSub;
              D->maxXSub=max+D->minXSub;
              D->nxSub=tmp;
           }
           rank++;
        }
        D->nySub=D->maxYSub=D->ny;
        D->minYSub=0;
     }
     else if(D->divDirection==2)
     {
        D->nySub=D->ny/nTasks;
        sub=D->nySub;
        remain=D->ny%nTasks;
        min=max=0;
        rank=0;
        while(rank<nTasks)
        {
           if(rank<remain)   tmp=sub+1;
           else              tmp=sub;
           min=max;
           max=min+tmp;
           if(myrank==rank)
           {
              D->minYSub=min;
              D->maxYSub=max;
              D->nySub=tmp;
           }
           rank++;
        }
        D->nxSub=D->nx;
        D->maxXSub=D->minXSub+D->nx;

        D->istart=2;
        D->iend=D->nxSub+2;
        D->jstart=2;
        D->jend=D->nySub+2;
     }
printf("nx=%d, nxSub=%d, minXSub=%d, maxXSub=%d\n",D->nx,D->nxSub,D->minXSub,D->maxXSub);
printf("ny=%d, nySub=%d, minYSub=%d, maxYSub=%d\n",D->ny,D->nySub,D->minYSub,D->maxYSub);

     if(D->fieldType==1)
     {
       D->fieldDSX = (FieldDSX **)malloc((D->nxSub+5)*sizeof(FieldDSX *));
       for(i=0; i<D->nxSub+5; i++) 
         D->fieldDSX[i] = (FieldDSX *)malloc((D->nySub+5)*sizeof(FieldDSX ));
     
       for(i=0; i<D->nxSub+5; i++)
         for(j=0; j<D->nySub+5; j++)
         {
           D->fieldDSX[i][j].Pr=0.0;//+Ext->Pr;
           D->fieldDSX[i][j].Pl=0.0;//+Ext->Pl;
           D->fieldDSX[i][j].E1=0.0;//+Ext->E1;
           D->fieldDSX[i][j].PrC=0.0;//+Ext->Pr;
           D->fieldDSX[i][j].PlC=0.0;//+Ext->Pl;
           D->fieldDSX[i][j].E1C=0.0;//+Ext->E1;
           D->fieldDSX[i][j].Sr=0.0;//+Ext->Sr;
           D->fieldDSX[i][j].Sl=0.0;//+Ext->Sl;
           D->fieldDSX[i][j].B1=0.0;//+Ext->B1;
           D->fieldDSX[i][j].SrC=0.0;//+Ext->Sr;
           D->fieldDSX[i][j].SlC=0.0;//+Ext->Sl;
           D->fieldDSX[i][j].B1C=0.0;//+Ext->B1;
           D->fieldDSX[i][j].J1=0.0;        
           D->fieldDSX[i][j].J2=0.0; 
           D->fieldDSX[i][j].J3=0.0; 
           D->fieldDSX[i][j].J1Old=0.0;        
           D->fieldDSX[i][j].J2Old=0.0; 
           D->fieldDSX[i][j].J3Old=0.0; 
         } 

       D->UpPML = (UpPML **)malloc((D->nxSub+5)*sizeof(UpPML *));
       for(i=0; i<D->nxSub+5; i++) 
         D->UpPML[i] = (UpPML *)malloc((10)*sizeof(UpPML ));     
       D->DnPML = (DnPML **)malloc((D->nxSub+5)*sizeof(DnPML *));
       for(i=0; i<D->nxSub+5; i++) 
         D->DnPML[i] = (DnPML *)malloc((10)*sizeof(DnPML ));
       for(i=0; i<D->nxSub+5; i++)
         for(j=0; j<10; j++)
         {
           D->UpPML[i][j].Bzx=0.0;
     }	//End of DSX
   
     D->particle = (Particle **)malloc((D->nxSub+3)*sizeof(Particle *));
     for(i=0; i<D->nxSub+3; i++) 
       D->particle[i] = (Particle *)malloc((D->nySub+3)*sizeof(Particle ));
    // setting up particle's pointer
     for(i=0; i<D->iend+1; i++)
       for(j=D->jstart-1; j<D->jend+1; j++)
       {
         D->particle[i][j].rho=0.0;  
         D->particle[i][j].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
         for(s=0; s<D->nSpecies; s++)
         {
           D->particle[i][j].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
           D->particle[i][j].head[s]->pt = NULL;
         }
       }
   
    // current J trasffering boundary
    numdataUp=(D->nx+5)*3*3;            
    numdataBt=(D->nx+5)*8;            
    D->upJ=(double *)malloc(numdataUp*sizeof(double ));
    D->btJ=(double *)malloc(numdataBt*sizeof(double ));


     D->boost = (Boost **)malloc((D->nxSub+5)*sizeof(Boost ));
     for(i=0; i<D->nxSub+5; i++) 
       D->boost[i] = (Boost *)malloc((D->nySub+5)*sizeof(Boost ));
     for(i=0; i<D->nxSub+5; i++)
       for(j=0; j<D->nySub+5; j++)
       {
         D->boost[i][j].E1=0.0;
         D->boost[i][j].B1=0.0;
         D->boost[i][j].Pr=0.0;
         D->boost[i][j].Pl=0.0;
         D->boost[i][j].Sr=0.0;
         D->boost[i][j].Sl=0.0;
       }

     D->probe = (Probe **)malloc(D->probeNum*sizeof(Probe *));
     for(i=0; i<D->probeNum; i++)
       D->probe[i] = (Probe *)malloc((D->maxStep+1)*sizeof(Probe ));
     for(i=0; i<D->probeNum; i++)
       for(j=0; j<=D->maxStep; j++)
       {
         D->probe[i][j].E1=0.0;
         D->probe[i][j].Pr=0.0;
         D->probe[i][j].Pl=0.0;
         D->probe[i][j].B1=0.0;
         D->probe[i][j].Sr=0.0;
         D->probe[i][j].Sl=0.0;
       }

}
                        
