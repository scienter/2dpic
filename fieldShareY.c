#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"
/*//
void MPI_TransferF_XplusFilter(Domain *D)
{
    int i,rank,numberData;
    int myrank, nTasks; 
    float *behindF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=5;
    behindF=(float *)malloc(numberData*sizeof(float )); 
    
    //Transferring even ~ odd cores             
    behindF[0]=D->field[D->nx+2].E1;
    behindF[1]=D->field[D->nx+2].Pr;
    behindF[2]=D->field[D->nx+2].Pl;
    behindF[3]=D->field[D->nx+2].Sr;
    behindF[4]=D->field[D->nx+2].Sl;
        
    if(myrank%2==1)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[0].E1=behindF[0];
       D->field[0].Pr=behindF[1];
       D->field[0].Pl=behindF[2];
       D->field[0].Sr=behindF[3];
       D->field[0].Sl=behindF[4];
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    behindF[0]=D->field[D->nx+2].E1;
    behindF[1]=D->field[D->nx+2].Pr;
    behindF[2]=D->field[D->nx+2].Pl;
    behindF[3]=D->field[D->nx+2].Sr;
    behindF[4]=D->field[D->nx+2].Sl;
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[0].E1=behindF[0];
       D->field[0].Pr=behindF[1];
       D->field[0].Pl=behindF[2];
       D->field[0].Sr=behindF[3];
       D->field[0].Sl=behindF[4];
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(behindF);
}
*/

void MPI_TransferJ_DSX_Yplus(Domain *D)
{
    int i,n,numberData,start,end,nx,nySub,ibegin,ibottom;
    int istart,iend,jstart,jend;
    int myrank, nTasks; 
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    FieldDSX **field;
    field=D->fieldDSX;

    nx=D->nx;
    nySub=D->nySub;
    numberData=3*3*(nx+5);
    ibegin=0;
    ibottom=nx+5;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    //Transferring even ~ odd cores 
    if(myrank%2==1)
    {
       MPI_Recv(D->upJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(n=0; n<3; n++)
       {
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J1+=D->upJ[i+start];
         start+=ibottom;
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J2+=D->upJ[i+start];
         start+=ibottom;
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J3+=D->upJ[i+start];
         start+=ibottom;
       }
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
    {
      start=0; 
      for(n=0; n<3; n++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J1;
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J2;
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J3;
        start+=ibottom;
      }
        
      MPI_Send(D->upJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(D->upJ,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(n=0; n<3; n++)
       {
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J1+=D->upJ[i+start];
         start+=ibottom;
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J2+=D->upJ[i+start];
         start+=ibottom;
         for(i=ibegin; i<ibottom; i++)
           field[i][jstart+n].J3+=D->upJ[i+start];
         start+=ibottom;
       }
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
    {
      start=0; 
      for(n=0; n<3; n++)
      {
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J1;
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J2;
        start+=ibottom;
        for(i=ibegin; i<ibottom; i++)
          D->upJ[start+i]=field[i][nySub+jstart+n].J3;
        start+=ibottom;
      }        
      MPI_Send(D->upJ,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD);             
    }     
    MPI_Barrier(MPI_COMM_WORLD);        
}


void MPI_TransferJ_Yminus(Domain *D)
{
    int i,n,numberData,start,end,nx,nySub,ibegin,ibottom;
    int istart,iend,jstart,jend;
    int myrank, nTasks; 
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    nx=D->nx;
    nySub=D->nySub;
    numberData=8*(nx+5);
    ibegin=0;
    ibottom=nx+5;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    //Transferring even ~ odd cores 
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(D->btJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J1+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J2+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J3+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J1+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J2+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J3+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+2].J3=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+2].J3Old=D->btJ[i+start];
    }
    else if(myrank%2==0 && myrank!=0)
    {
      start=0; 
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J1;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J2;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J1;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J2;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][2].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][2].J3Old;
        
      MPI_Send(D->btJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(D->btJ,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J1+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J2+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub].J3+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J1+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J2+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+1].J3+=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+2].J3=D->btJ[i+start];
       start+=ibottom;
       for(i=ibegin; i<ibottom; i++)
         field[i][nySub+2].J3Old=D->btJ[i+start];
    }
    else if(myrank%2==1)
    {
      start=0; 
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J1;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J2;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][0].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J1;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J2;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][1].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][2].J3;
      start+=ibottom;
      for(i=ibegin; i<ibottom; i++)
        D->btJ[start+i]=field[i][2].J3Old;
        
      MPI_Send(D->btJ,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD);             
    }
     
    MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_TransferF_Yminus(Domain *D)
{
    int i,numberData,start,end,nx,nySub,ibegin,ibottom;
    int istart,iend,jstart,jend;
    int myrank, nTasks; 
    float *btF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    nx=D->nx;
    nySub=D->nySub;
    ibegin=1;
    ibottom=nx+3;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;


    numberData=6*(nx+2);
    btF=(float *)malloc(numberData*sizeof(float ));             

    //Transferring even ~ odd cores 
    start=0; 
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].E1;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].B1;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].Pr;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].Pl;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].Sr;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].Sl;
        
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(btF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
          D->field[i][nySub+jstart].E1=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          D->field[i][nySub+istart].B1=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          D->field[i][nySub+istart].Pr=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          D->field[i][nySub+istart].Pl=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          D->field[i][nySub+istart].Sr=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          D->field[i][nySub+istart].Sl=btF[i+start];
    }
    else if(myrank%2==1)
       MPI_Send(btF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
/*
    //Transferring odd ~ even cores             
    start=0; 
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].E1;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].B1;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].Pr;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].Pl;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].Sr;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].Sl;
        
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(btF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+jstart].E1=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+jstart].B1=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+jstart].Pr=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+jstart].Pl=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+jstart].Sr=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+jstart].Sl=btF[i+start];
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(btF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             

    MPI_Barrier(MPI_COMM_WORLD);
*/    
    free(btF);
}

void MPI_TransferF_YminusC(Domain *D)
{
    int i,numberData,start,end,nx,nySub,ibegin,ibottom;
    int istart,iend,jstart,jend;
    int myrank, nTasks; 
    float *btF;
    MPI_Status status;         

    nx=D->nx;
    nySub=D->nySub;
    ibegin=1;
    ibottom=nx+3;   

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=2*(nx+2);
    btF=(float *)malloc(numberData*sizeof(float ));             

    //Transferring even ~ odd cores 
    start=0; 
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].E1C;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].B1C;
        
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(btF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+2].E1C=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+2].B1C=btF[i+start];
    }
    else if(myrank%2==1)
       MPI_Send(btF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].E1C;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       btF[start+i]=field[i][jstart].B1C;
        
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(btF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+jstart].E1C=btF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][nySub+jstart].B1C=btF[i+start];
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(btF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(btF);
}

void MPI_TransferF_Yplus(Domain *D)
{
    int i,numberData,start,end,nx,nySub,ibegin,ibottom;
    int istart,iend,jstart,jend;
    int myrank, nTasks; 
    float *upF;
    MPI_Status status;         
   
    nx=D->nx;
    nySub=D->nySub;
    ibegin=1;
    ibottom=nx+3;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=6*(nx+2);
    upF=(float *)malloc(numberData*sizeof(float ));             

    //Transferring even ~ odd cores 
    
    if(myrank%2==1)
    {
       MPI_Recv(upF,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].Pr=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].Pl=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].Sr=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].Sl=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].E1=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].B1=upF[i+start];
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
    {
      start=0; 
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].Pr;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].Pl;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].Sr;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].Sl;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].E1;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].B1;
       MPI_Send(upF,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);             
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(upF,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].Pr=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].Pl=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].Sr=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].Sl=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].E1=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].B1=upF[i+start];
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
    {
      start=0; 
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].Pr;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].Pl;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].Sr;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].Sl;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].E1;
      start+=nx+2;
      for(i=ibegin; i<ibottom; i++)
         upF[start+i]=field[i][nySub+1].B1;
       MPI_Send(upF,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);             
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(upF);
}

void MPI_TransferF_YplusC(Domain *D)
{
    int i,numberData,start,end,nx,nySub,ibegin,ibottom;
    int istart,iend,jstart,jend;
    int myrank, nTasks; 
    float *upF;
    MPI_Status status;         
   
    nx=D->nx;
    nySub=D->nySub;
    ibegin=1;
    ibottom=nx+3;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=4*(D->nx+2);
    upF=(float *)malloc(numberData*sizeof(float ));             

    //Transferring even ~ odd cores 
    start=0; 
    for(i=ibegin; i<ibottom; i++)
       upF[start+i]=field[i][nySub+1].PrC;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       upF[start+i]=field[i][nySub+1].PlC;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       upF[start+i]=field[i][nySub+1].SrC;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       upF[start+i]=field[i][nySub+1].SlC;
    
    if(myrank%2==1)
    {
       MPI_Recv(upF,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].PrC=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].PlC=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].SrC=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].SlC=upF[i+start];
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(upF,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=ibegin; i<ibottom; i++)
       upF[start+i]=field[i][nySub+1].PrC;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       upF[start+i]=field[i][nySub+1].PlC;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       upF[start+i]=field[i][nySub+1].SrC;
    start+=nx+2;
    for(i=ibegin; i<ibottom; i++)
       upF[start+i]=field[i][nySub+1].SlC;
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(upF,numberData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);  
       start=0;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].PrC=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].PlC=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].SrC=upF[i+start];
       start+=nx+2;
       for(i=ibegin; i<ibottom; i++)
          field[i][1].SlC=upF[i+start];
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(upF,numberData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(upF);
}

