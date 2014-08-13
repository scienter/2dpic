#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

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
    behindF[0]=D->field[D->nxSub].E1;
    behindF[1]=D->field[D->nxSub].Pr;
    behindF[2]=D->field[D->nxSub].Pl;
    behindF[3]=D->field[D->nxSub].Sr;
    behindF[4]=D->field[D->nxSub].Sl;
        
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
    behindF[0]=D->field[D->nxSub].E1;
    behindF[1]=D->field[D->nxSub].Pr;
    behindF[2]=D->field[D->nxSub].Pl;
    behindF[3]=D->field[D->nxSub].Sr;
    behindF[4]=D->field[D->nxSub].Sl;
        
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


void MPI_TransferF_Xminus(Domain *D)
{
    int i,rank,numberData;
    int myrank, nTasks; 
    float *frontF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=3;
    frontF=(float *)malloc(numberData*sizeof(float ));             

    //Transferring even ~ odd cores             
    frontF[0]=D->field[1].Pl;
    frontF[1]=D->field[1].Sl;
    frontF[2]=D->field[1].E1;
        
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(frontF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+1].Pl=frontF[0];
       D->field[D->nxSub+1].Sl=frontF[1];
       D->field[D->nxSub+1].E1=frontF[2];
    }
    else if(myrank%2==1)
       MPI_Send(frontF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    frontF[0]=D->field[1].Pl;
    frontF[1]=D->field[1].Sl;
    frontF[2]=D->field[1].E1;

        
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(frontF,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+1].Pl=frontF[0];
       D->field[D->nxSub+1].Sl=frontF[1];
       D->field[D->nxSub+1].E1=frontF[2];
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(frontF,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(frontF);
}

void MPI_TransferF_Xplus(Domain *D)
{
    int i,rank,numberData;
    int myrank, nTasks; 
    float *behindF;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=2;
    behindF=(float *)malloc(numberData*sizeof(float )); 
    
    //Arrange core orders.
    D->beforeCore=myrank-1;
    D->nextCore=myrank+1;           
    
    //Transferring even ~ odd cores             
    behindF[0]=D->field[D->nxSub+1].Pr;
    behindF[1]=D->field[D->nxSub+1].Sr;
        
    if(myrank%2==1)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[1].Pr=behindF[0];
       D->field[1].Sr=behindF[1];
    }
    else if(myrank%2==0 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    behindF[0]=D->field[D->nxSub+1].Pr;
    behindF[1]=D->field[D->nxSub+1].Sr;
        
    if(myrank%2==0 && myrank!=0)
    {
       MPI_Recv(behindF,numberData, MPI_FLOAT, D->beforeCore, D->beforeCore, MPI_COMM_WORLD,&status);  
       D->field[1].Pr=behindF[0];
       D->field[1].Sr=behindF[1];
    }
    else if(myrank%2==1 && myrank!=nTasks-1)
       MPI_Send(behindF,numberData, MPI_FLOAT, D->nextCore, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(behindF);
}

void MPI_TransferJ_Moving(Domain *D)
{
    int i,rank,numberData;
    int myrank, nTasks; 
    float *frontJ;
    MPI_Status status;         
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=3;
    frontJ=(float *)malloc(numberData*sizeof(float ));             

    //Transferring even ~ odd cores             
    frontJ[0]=D->field[1].J1;
    frontJ[1]=D->field[1].J2;
    frontJ[2]=D->field[1].J3;
        
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       MPI_Recv(frontJ,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+1].J1=frontJ[0];
       D->field[D->nxSub+1].J2=frontJ[1];
       D->field[D->nxSub+1].J3=frontJ[2];
    }
    else if(myrank%2==1)
       MPI_Send(frontJ,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    frontJ[0]=D->field[1].J1;
    frontJ[1]=D->field[1].J2;
    frontJ[2]=D->field[1].J3;
        
    if(myrank%2==1 && myrank!=nTasks-1)
    {
       MPI_Recv(frontJ,numberData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);  
       D->field[D->nxSub+1].J1=frontJ[0];
       D->field[D->nxSub+1].J2=frontJ[1];
       D->field[D->nxSub+1].J3=frontJ[2];
    }
    else if(myrank%2==0 && myrank!=0)
       MPI_Send(frontJ,numberData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    
    free(frontJ);
}
