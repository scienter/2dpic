#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>


void MPI_TransferP_Yplus(Domain *D)
{
   int i,j,n,s,numP,cnt,totalData,sendData=9,nxSub,nySub,nx;
   int istart,iend,jstart,jend;
   int myrank, nTasks;
   Particle **particle;
   particle=D->particle;     
   float *upP;
   ptclList *p,*tmp,*New;
   MPI_Status status;         

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

   nxSub=D->nxSub;
   nySub=D->nySub;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

    //Even -> odd
    if(myrank%2==0 && myrank!=nTasks-1)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart-1; i<=iend; i++)
        {
          p=particle[i][nySub+2].head[s]->pt;
          while(p)   {
            p=p->next;
            numP++;
          } 
        }
      MPI_Send(&numP,1, MPI_INT, myrank+1, myrank, MPI_COMM_WORLD);    
    }    
    else if(myrank%2==1) 
      MPI_Recv(&numP,1, MPI_INT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    upP=(float *)malloc(totalData*sizeof(float ));             
//     for(i=0; i<totalData; i++)
//       upP[i]=0.0;
 

    if(myrank%2==0 && myrank!=nTasks-1)
    {    
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart-1; i<=iend; i++)
        {
          p=particle[i][nySub+2].head[s]->pt;
          while(p)   
          {
            upP[n*sendData+0]=p->x+i;     
            upP[n*sendData+1]=p->oldX;  
            upP[n*sendData+2]=p->y;  
            upP[n*sendData+3]=p->oldY-nySub;  
            upP[n*sendData+4]=p->p1;  
            upP[n*sendData+5]=p->p2;  
            upP[n*sendData+6]=p->p3;  
            upP[n*sendData+7]=p->index;  
            upP[n*sendData+8]=(float)s;  
            p=p->next;
            n++;
          }
        }
      MPI_Send(upP,totalData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD); 
    }

    else if(myrank%2==1) 
    {    
      MPI_Recv(upP,totalData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(upP[n*sendData+0]);
        s=(int)(upP[n*sendData+8]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][jstart].head[s]->pt;
        particle[i][jstart].head[s]->pt = New;             
        New->x=upP[n*sendData+0]-i;     
        New->oldX=upP[n*sendData+1];     
        New->y=upP[n*sendData+2];     
        New->oldY=upP[n*sendData+3];     
        New->p1=upP[n*sendData+4];
        New->p2=upP[n*sendData+5];
        New->p3=upP[n*sendData+6];
        New->index=upP[n*sendData+7];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(upP);

    //Odd -> evem
    if(myrank%2==1 && myrank!=nTasks-1)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart-1; i<=iend; i++)
        {
          p=particle[i][nySub+2].head[s]->pt;
          while(p)   {
            p=p->next;
            numP++;
          }
        } 
      MPI_Send(&numP,1, MPI_INT, myrank+1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==0 && myrank!=0) 
      MPI_Recv(&numP,1, MPI_INT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    upP=(float *)malloc(totalData*sizeof(float ));             
//     for(i=0; i<totalData; i++)
//       upP[i]=0.0;


    if(myrank%2==1 && myrank!=nTasks-1)
    {
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart-1; i<=iend; i++)
        {
          p=particle[i][nySub+2].head[s]->pt;
          while(p)   
          {
            upP[n*sendData+0]=p->x+i;     
            upP[n*sendData+1]=p->oldX;  
            upP[n*sendData+2]=p->y;  
            upP[n*sendData+3]=p->oldY-nySub;  
            upP[n*sendData+4]=p->p1;  
            upP[n*sendData+5]=p->p2;  
            upP[n*sendData+6]=p->p3;  
            upP[n*sendData+7]=p->index;   
            upP[n*sendData+8]=(float)s;   
            p=p->next;
            n++;
          }
        }
      MPI_Send(upP,totalData, MPI_FLOAT, myrank+1, myrank, MPI_COMM_WORLD); 
    }
 
    else if(myrank%2==0 && myrank!=0) 
    {    
      MPI_Recv(upP,totalData, MPI_FLOAT, myrank-1, myrank-1, MPI_COMM_WORLD,&status);      
      for(n=0; n<numP; n++)
      {
        i=(int)(upP[n*sendData+0]);
        s=(int)(upP[n*sendData+8]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][jstart].head[s]->pt;
        particle[i][jstart].head[s]->pt = New;             
        New->x=upP[n*sendData+0]-i;     
        New->oldX=upP[n*sendData+1];     
        New->y=upP[n*sendData+2];     
        New->oldY=upP[n*sendData+3];     
        New->p1=upP[n*sendData+4];
        New->p2=upP[n*sendData+5];
        New->p3=upP[n*sendData+6];
        New->index=upP[n*sendData+7];
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(upP);

}


void MPI_TransferP_Yminus(Domain *D)
{
   int i,j,n,s,numP,cnt,totalData,sendData=9,nxSub,nySub;
   int istart,iend,jstart,jend;
   int myrank, nTasks;
   Particle **particle;
   particle=D->particle;     
   float *btP;
   ptclList *p,*tmp,*New;
   MPI_Status status;         

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);   
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   nxSub=D->nxSub;
   nySub=D->nySub;

    //Even -> odd
    if(myrank%2==0 && myrank!=0)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart-1; i<=iend; i++)
        {
          p=particle[i][istart-1].head[s]->pt;
          while(p)   {
            p=p->next;
            numP++;
          } 
        }
      MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }    
    else if(myrank%2==1 && myrank!=nTasks-1) 
      MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);        
    MPI_Barrier(MPI_COMM_WORLD);    

    totalData=numP*sendData;
    btP=(float *)malloc(totalData*sizeof(float ));             
//     for(i=0; i<totalData; i++)
//       upP[i]=0.0;
 
    if(myrank%2==0 && myrank!=0)
    {
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart-1; i<=iend; i++)
        {
          p=particle[i][istart-1].head[s]->pt;
          while(p)   
          {
            btP[n*sendData+0]=p->x+i;     
            btP[n*sendData+1]=p->oldX;  
            btP[n*sendData+2]=p->y;  
            btP[n*sendData+3]=p->oldY;  
            btP[n*sendData+4]=p->p1;  
            btP[n*sendData+5]=p->p2;  
            btP[n*sendData+6]=p->p3;  
            btP[n*sendData+7]=p->index;  
            btP[n*sendData+8]=(float)s;  
            p=p->next;
            n++;
          }
        }
      MPI_Send(btP,totalData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD); 
    }
    else if(myrank%2==1 && myrank!=nTasks-1) 
    {    
      MPI_Recv(btP,totalData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
      n=0;
      while(numP>0)
      {
        i=(int)(btP[n*sendData+0]);
        s=(int)(btP[n*sendData+8]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][D->nySub+1].head[s]->pt;
        particle[i][D->nySub+1].head[s]->pt = New;             
        New->x=btP[n*sendData+0]-i;     
        New->oldX=btP[n*sendData+1];     
        New->y=btP[n*sendData+2];     
        New->oldY=btP[n*sendData+3]+nySub;     
        New->p1=btP[n*sendData+4];
        New->p2=btP[n*sendData+5];
        New->p3=btP[n*sendData+6];
        New->index=btP[n*sendData+7];
        n++;
        numP--;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(btP);

    //Odd -> evem
    if(myrank%2==1)
    {
      numP=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart-1; i<=iend; i++)
        {
          p=particle[i][1].head[s]->pt;
          while(p)   {
            p=p->next;
            numP++;
          }
        } 
      MPI_Send(&numP,1, MPI_INT, myrank-1, myrank, MPI_COMM_WORLD);    
    }
    else if(myrank%2==0 && myrank!=nTasks-1) 
      MPI_Recv(&numP,1, MPI_INT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
    MPI_Barrier(MPI_COMM_WORLD);    


    totalData=numP*sendData;
    btP=(float *)malloc(totalData*sizeof(float ));             
//    for(i=0; i<totalData; i++)
//      btP[i]=0.0;

    if(myrank%2==1)
    {
      n=0;
      for(s=0; s<D->nSpecies; s++)
        for(i=istart-1; i<=iend; i++)
        {
          p=particle[i][1].head[s]->pt;
          while(p)   
          {
            btP[n*sendData+0]=p->x+i;     
            btP[n*sendData+1]=p->oldX;  
            btP[n*sendData+2]=p->y;  
            btP[n*sendData+3]=p->oldY;  
            btP[n*sendData+4]=p->p1;  
            btP[n*sendData+5]=p->p2;  
            btP[n*sendData+6]=p->p3;  
            btP[n*sendData+7]=p->index;  
            btP[n*sendData+8]=(float)s;  
            p=p->next;
            n++;
          }
        }
      MPI_Send(btP,totalData, MPI_FLOAT, myrank-1, myrank, MPI_COMM_WORLD); 
    }

    else if(myrank%2==0 && myrank!=nTasks-1) 
    {    
      MPI_Recv(btP,totalData, MPI_FLOAT, myrank+1, myrank+1, MPI_COMM_WORLD,&status);      
      n=0;
      while(numP>0)
      {

        i=(int)(btP[n*sendData+0]);
        s=(int)(btP[n*sendData+8]);
        New = (ptclList *)malloc(sizeof(ptclList)); 
        New->next = particle[i][D->nySub+1].head[s]->pt;
        particle[i][D->nySub+1].head[s]->pt = New;             
        New->x=btP[n*sendData+0]-i;     
        New->oldX=btP[n*sendData+1];     
        New->y=btP[n*sendData+2];     
        New->oldY=btP[n*sendData+3]+nySub;     
        New->p1=btP[n*sendData+4];
        New->p2=btP[n*sendData+5];
        New->p3=btP[n*sendData+6];
        New->index=btP[n*sendData+7];
        numP--;
        n++;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(btP);

}

