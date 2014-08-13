#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>

FieldElement *field;

void MPI_TransferP_Moving(Domain *D)
{
    int i,numP=0,cnt,numberData,s;
    int myrank, nTasks;
    field=D->field;     
    double *frontP;
    ptclList *p,*tmp,*New;
    MPI_Status status;         
   
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=9;
    frontP=(double *)malloc(numberData*sizeof(double ));             
    for(i=0; i<9; i++)
       frontP[i]=0;

  for(s=0; s<D->nSpecies; s++)
  {
    //Even - odd
    if(myrank%2==0 && myrank!=0)
    {
       numP=0;
       p=field[1].head[s]->pt;
       while(p)   {
          p=p->next;
          numP++;
       }     
       MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==1 && myrank!=D->L-1) 
    {    
       MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    }
    MPI_Barrier(MPI_COMM_WORLD);    

    for(i=0; i<numP; i++)
    {
       if(myrank%2==0 && myrank!=0)
       {
          p=field[1].head[s]->pt;
          frontP[0]=p->x1;     
          frontP[1]=p->oldX1;  
          frontP[2]=p->p1;     
          frontP[3]=p->p2;     
          frontP[4]=p->p3;     
          frontP[5]=p->q;      
          frontP[6]=p->m;      
          frontP[7]=p->index;  
          frontP[8]=p->np2c; 
          tmp=p->next;
          p->next=NULL;
          free(p);
          field[1].head[s]->pt=tmp; 

          MPI_Send(frontP,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD); 
       }
    
       else if(myrank%2==1 && myrank!=D->L-1) 
       {    
          MPI_Recv(frontP,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
          New = (ptclList *)malloc(sizeof(ptclList)); 
          New->x1=frontP[0];
          New->oldX1=frontP[1]+D->nxSub;
          New->p1=frontP[2];
          New->p2=frontP[3];
          New->p3=frontP[4];
          New->q=frontP[5];
          New->m=frontP[6];
          New->index=frontP[7];
          New->np2c=frontP[8];

          New->next = field[D->nxSub+1].head[s]->pt;
          field[D->nxSub+1].head[s]->pt = New;      
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
       
    //Odd - evem
    if(myrank%2==1)
    {
       numP=0;
       p=field[1].head[s]->pt;
       while(p)   {
          p=p->next;
          numP++;
       }     
       MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==0 && myrank!=D->L-1) 
    {    
       MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    }
    MPI_Barrier(MPI_COMM_WORLD);    

    for(i=0; i<numP; i++)
    {
       if(myrank%2==1)
       {
          p=field[1].head[s]->pt;
          frontP[0]=p->x1;     
          frontP[1]=p->oldX1;  
          frontP[2]=p->p1;     
          frontP[3]=p->p2;     
          frontP[4]=p->p3;     
          frontP[5]=p->q;      
          frontP[6]=p->m;      
          frontP[7]=p->index;  
          frontP[8]=p->np2c;  
          tmp=p->next;
          p->next=NULL;
          free(p);
          field[1].head[s]->pt=tmp; 

          MPI_Send(frontP,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD); 
       }
    
       else if(myrank%2==0 && myrank!=D->L-1)  
       {    
          MPI_Recv(frontP,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
          New = (ptclList *)malloc(sizeof(ptclList)); 
          New->x1=frontP[0];
          New->oldX1=frontP[1]+D->nxSub;
          New->p1=frontP[2];
          New->p2=frontP[3];
          New->p3=frontP[4];
          New->q=frontP[5];
          New->m=frontP[6];
          New->index=frontP[7];
          New->np2c=frontP[8];

          New->next = field[D->nxSub+1].head[s]->pt;
          field[D->nxSub+1].head[s]->pt = New;      
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }	//End of for(s)    

  free(frontP);
}

void MPI_TransferP_Xplus(Domain *D)
{
    int i,s,numP=0,cnt,numberData;
    int myrank, nTasks;
    field=D->field;     
    double *behindP;
    ptclList *p,*tmp,*New;
    MPI_Status status;         

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=9;
    behindP=(double *)malloc(numberData*sizeof(double ));             
    for(i=0; i<9; i++)
       behindP[i]=0.0;

  for(s=0; s<D->nSpecies; s++)
  {
    //Even -> odd
    if(myrank%2==0 && myrank!=nTasks-1)
    {
       numP=0;
       p=field[D->nxSub+1].head[s]->pt;
       while(p)   {
          p=p->next;
          numP++;
       }     
       MPI_Send(&numP,1, MPI_INT, myrank+1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==1) 
    {    
       MPI_Recv(&numP,1, MPI_INT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
    }
    MPI_Barrier(MPI_COMM_WORLD);    
 
    for(i=0; i<numP; i++)
    {
       if(myrank%2==0 && myrank!=nTasks-1)
       {
          p=field[D->nxSub+1].head[s]->pt;
          behindP[0]=p->x1;     
          behindP[1]=p->oldX1-D->nxSub;  
          behindP[2]=p->p1;     
          behindP[3]=p->p2;     
          behindP[4]=p->p3;     
          behindP[5]=p->q;      
          behindP[6]=p->m;      
          behindP[7]=p->index;  
          behindP[8]=p->np2c;  
          tmp=p->next;
          p->next=NULL;
          free(p);
          field[D->nxSub+1].head[s]->pt=tmp; 

          MPI_Send(behindP,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD); 
       }
    
       else if(myrank%2==1) 
       {    
          MPI_Recv(behindP,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
          New = (ptclList *)malloc(sizeof(ptclList)); 
          New->x1=behindP[0];
          New->oldX1=behindP[1];
          New->p1=behindP[2];
          New->p2=behindP[3];
          New->p3=behindP[4];
          New->q=behindP[5];
          New->m=behindP[6];
          New->index=behindP[7];
          New->np2c=behindP[8];

          New->next = field[1].head[s]->pt;
          field[1].head[s]->pt = New;      
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //Odd -> evem
    if(myrank%2==1 && myrank!=D->L-1)
    {
       numP=0;
       p=field[D->nxSub+1].head[s]->pt;
       while(p)   {
          p=p->next;
          numP++;
       }     
       MPI_Send(&numP,1, MPI_INT, myrank+1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==0 && myrank!=0) 
    {    
       MPI_Recv(&numP,1, MPI_INT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
    }
    MPI_Barrier(MPI_COMM_WORLD);    

    for(i=0; i<numP; i++)
    {
       if(myrank%2==1 && myrank!=D->L-1)
       {
          p=field[D->nxSub+1].head[s]->pt;
          behindP[0]=p->x1;     
          behindP[1]=p->oldX1-D->nxSub;  
          behindP[2]=p->p1;     
          behindP[3]=p->p2;     
          behindP[4]=p->p3;     
          behindP[5]=p->q;      
          behindP[6]=p->m;      
          behindP[7]=p->index;  
          behindP[8]=p->np2c;  
          tmp=p->next;
          p->next=NULL;
          free(p);
          field[D->nxSub+1].head[s]->pt=tmp; 

          MPI_Send(behindP,numberData, MPI_DOUBLE, myrank+1, myrank, MPI_COMM_WORLD); 
       }
    
       else if(myrank%2==0 && myrank!=0) 
       {    
          MPI_Recv(behindP,numberData, MPI_DOUBLE, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
          New = (ptclList *)malloc(sizeof(ptclList)); 
          New->x1=behindP[0];
          New->oldX1=behindP[1];
          New->p1=behindP[2];
          New->p2=behindP[3];
          New->p3=behindP[4];
          New->q=behindP[5];
          New->m=behindP[6];
          New->index=behindP[7];
          New->np2c=behindP[8];

          New->next = field[1].head[s]->pt;
          field[1].head[s]->pt = New;      
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }	//End of for(s)

  free(behindP);
}



void MPI_TransferP_Xminus(Domain *D)
{
    int i,numP=0,cnt,numberData,s;
    int myrank, nTasks;
    field=D->field;     
    double *frontP;
    ptclList *p,*tmp,*New;
    MPI_Status status;         
   
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=9;
    frontP=(double *)malloc(numberData*sizeof(double ));             
    for(i=0; i<9; i++)
       frontP[i]=0;

  for(s=0; s<D->nSpecies; s++)
  {
    //Even - odd
    if(myrank%2==0 && myrank!=0)
    {
       numP=0;
       p=field[0].head[s]->pt;
       while(p)   {
          p=p->next;
          numP++;
       }     
       MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==1 && myrank!=D->L-1) 
    {    
       MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    }
    MPI_Barrier(MPI_COMM_WORLD);    

    for(i=0; i<numP; i++)
    {
       if(myrank%2==0 && myrank!=0)
       {
          p=field[0].head[s]->pt;
          frontP[0]=p->x1;     
          frontP[1]=p->oldX1;  
          frontP[2]=p->p1;     
          frontP[3]=p->p2;     
          frontP[4]=p->p3;     
          frontP[5]=p->q;      
          frontP[6]=p->m;      
          frontP[7]=p->index;  
          frontP[8]=p->np2c;  
          tmp=p->next;
          p->next=NULL;
          free(p);
          field[0].head[s]->pt=tmp; 

          MPI_Send(frontP,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD); 
       }
    
       else if(myrank%2==1 && myrank!=D->L-1) 
       {    
          MPI_Recv(frontP,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
          New = (ptclList *)malloc(sizeof(ptclList)); 
          New->x1=frontP[0];
          New->oldX1=frontP[1]+D->nxSub;
          New->p1=frontP[2];
          New->p2=frontP[3];
          New->p3=frontP[4];
          New->q=frontP[5];
          New->m=frontP[6];
          New->index=frontP[7];
          New->np2c=frontP[8];

          New->next = field[D->nxSub].head[s]->pt;
          field[D->nxSub].head[s]->pt = New;      
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
       
    //Odd - evem
    if(myrank%2==1)
    {
       numP=0;
       p=field[0].head[s]->pt;
       while(p)   {
          p=p->next;
          numP++;
       }     
       MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==0 && myrank!=D->L-1) 
    {    
       MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    }
    MPI_Barrier(MPI_COMM_WORLD);    

    for(i=0; i<numP; i++)
    {
       if(myrank%2==1)
       {
          p=field[0].head[s]->pt;
          frontP[0]=p->x1;     
          frontP[1]=p->oldX1;  
          frontP[2]=p->p1;     
          frontP[3]=p->p2;     
          frontP[4]=p->p3;     
          frontP[5]=p->q;      
          frontP[6]=p->m;      
          frontP[7]=p->index;  
          frontP[8]=p->np2c;  
          tmp=p->next;
          p->next=NULL;
          free(p);
          field[0].head[s]->pt=tmp; 

          MPI_Send(frontP,numberData, MPI_DOUBLE, myrank-1, myrank, MPI_COMM_WORLD); 
       }
    
       else if(myrank%2==0 && myrank!=D->L-1)  
       {    
          MPI_Recv(frontP,numberData, MPI_DOUBLE, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
          New = (ptclList *)malloc(sizeof(ptclList)); 
          New->x1=frontP[0];
          New->oldX1=frontP[1]+D->nxSub;
          New->p1=frontP[2];
          New->p2=frontP[3];
          New->p3=frontP[4];
          New->q=frontP[5];
          New->m=frontP[6];
          New->index=frontP[7];
          New->np2c=frontP[8];

          New->next = field[D->nxSub].head[s]->pt;
          field[D->nxSub].head[s]->pt = New;      
       }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }	//End of for(s)    

  free(frontP);
}

