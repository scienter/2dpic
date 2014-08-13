#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"

void periodY(Domain *D)
{
    int i,j,start,end,numberData;  
    float dx,dy,dt,preB1C;
    FieldElement **field;
    field=D->field;
    int myrank, nTasks; 
    MPI_Status status;         
    float *btF, *upF;
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=4*(D->nx+2);
    btF=(float *)malloc(numberData*sizeof(float ));             

	 //rank(nTasks-1) -> rank(0)
	 start=0; 
	 for(i=0; i<=D->nx+1; i++)
	    btF[start+i]=D->field[i][D->nySub].Pr;
	 start+=D->nx;
	 for(i=0; i<=D->nx+1; i++)
		 btF[start+i]=D->field[i][D->nySub].Pl;
	 start+=D->nx;
	 for(i=0; i<=D->nx+1; i++)
		 btF[start+i]=D->field[i][D->nySub].Sr;
	 start+=D->nx;
	 for(i=0; i<=D->nx+1; i++)
	    btF[start+i]=D->field[i][D->nySub].Sl;

    if(myrank==nTasks-1)
       MPI_Send(btF,numberData, MPI_FLOAT, 0, myrank, MPI_COMM_WORLD);
    else if(myrank==0)         
	 {
	    MPI_Recv(btF,numberData, MPI_FLOAT, nTasks-1, nTasks-1, MPI_COMM_WORLD,&status);  
		 start=0;
		 end=start+D->nx;
		 for(i=0; i<=D->nx+1; i++)
		    D->field[i][0].Pr=btF[i+start];
		 start=end;
		 end=start+D->nx;
		 for(i=0; i<=D->nx+1; i++)
		    D->field[i][0].Pl=btF[i+start];
		 start=end;
		 end=start+D->nx;
		 for(i=0; i<=D->nx+1; i++)
		    D->field[i][0].Sr=btF[i+start];
		 start=end;
		 end=start+D->nx;
		 for(i=0; i<=D->nx+1; i++)
		    D->field[i][0].Sl=btF[i+start];
	 }

    numberData=2*(D->nx+2);
    upF=(float *)malloc(numberData*sizeof(float ));             

    //rank(0) -> rank(nTasks-1) 
    start=0; 
    for(i=0; i<=D->nx+1; i++)
       upF[start+i]=D->field[i][1].B1;
    start+=D->nx;
    for(i=0; i<=D->nx+1; i++)
       upF[start+i]=D->field[i][1].E1;
        
    if(myrank==0)
       MPI_Send(upF,numberData, MPI_FLOAT, nTasks-1, myrank, MPI_COMM_WORLD); 
    else if(myrank==nTasks-1)         
    {
       MPI_Recv(upF,numberData, MPI_FLOAT, 0, 0, MPI_COMM_WORLD,&status);  
       start=0;
       end=start+D->nx;
       for(i=0; i<=D->nx+1; i++)
          D->field[i][D->nySub+1].B1=upF[i+start];
       start=end;
       end=start+D->nx;
       for(i=0; i<=D->nx+1; i++)
          D->field[i][D->nySub+1].E1=upF[i+start];
    }

    free(btF);
    free(upF);
}

void periodYC(Domain *D)
{
    int i,j,start,end,numberData;  
    float dx,dy,dt,preB1C;
    FieldElement **field;
    field=D->field;
    int myrank, nTasks; 
    MPI_Status status;         
    float *btF, *upF;
   
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    numberData=4*(D->nx+2);
    btF=(float *)malloc(numberData*sizeof(float ));             

	 //rank(nTasks-1) -> rank(0)
	 start=0; 
	 for(i=0; i<=D->nx+1; i++)
	    btF[start+i]=D->field[i][D->nySub].PrC;
	 start+=D->nx;
	 for(i=0; i<=D->nx+1; i++)
		 btF[start+i]=D->field[i][D->nySub].PlC;
	 start+=D->nx;
	 for(i=0; i<=D->nx+1; i++)
		 btF[start+i]=D->field[i][D->nySub].SrC;
	 start+=D->nx;
	 for(i=0; i<=D->nx+1; i++)
	    btF[start+i]=D->field[i][D->nySub].SlC;

    if(myrank==nTasks-1)
       MPI_Send(btF,numberData, MPI_FLOAT, 0, myrank, MPI_COMM_WORLD);
    else if(myrank==0)         
	 {
	    MPI_Recv(btF,numberData, MPI_FLOAT, nTasks-1, nTasks-1, MPI_COMM_WORLD,&status);  
		 start=0;
		 end=start+D->nx;
		 for(i=0; i<=D->nx+1; i++)
		    D->field[i][0].PrC=btF[i+start];
		 start=end;
		 end=start+D->nx;
		 for(i=0; i<=D->nx+1; i++)
		    D->field[i][0].PlC=btF[i+start];
		 start=end;
		 end=start+D->nx;
		 for(i=0; i<=D->nx+1; i++)
		    D->field[i][0].SrC=btF[i+start];
		 start=end;
		 end=start+D->nx;
		 for(i=0; i<=D->nx+1; i++)
		    D->field[i][0].SlC=btF[i+start];
	 }

    numberData=2*(D->nx+2);
    upF=(float *)malloc(numberData*sizeof(float ));             

    //rank(0) -> rank(nTasks-1) 
    start=0; 
    for(i=0; i<=D->nx+1; i++)
       upF[start+i]=D->field[i][1].B1C;
    start+=D->nx;
    for(i=0; i<=D->nx+1; i++)
       upF[start+i]=D->field[i][1].E1C;
        
    if(myrank==0)
       MPI_Send(upF,numberData, MPI_FLOAT, nTasks-1, myrank, MPI_COMM_WORLD); 
    else if(myrank==nTasks-1)         
    {
       MPI_Recv(upF,numberData, MPI_FLOAT, 0, 0, MPI_COMM_WORLD,&status);  
       start=0;
       end=start+D->nx;
       for(i=0; i<=D->nx+1; i++)
          D->field[i][D->nySub+1].B1C=upF[i+start];
       start=end;
       end=start+D->nx;
       for(i=0; i<=D->nx+1; i++)
          D->field[i][D->nySub+1].E1C=upF[i+start];
    }

    free(btF);
    free(upF);
}

void periodY1core(Domain *D)
{
    int i,j,start,end,numberData;  
    float dx,dy,dt,preB1C;
    FieldElement **field;
    field=D->field;
   
    for(i=0; i<=D->nxSub+1; i++)
    {
       field[i][0].Pr=field[i][D->nySub].Pr;
       field[i][0].Pl=field[i][D->nySub].Pl;
       field[i][0].Sr=field[i][D->nySub].Sr;
       field[i][0].Sl=field[i][D->nySub].Sl;
       field[i][D->nySub+1].B1=field[i][1].B1;
       field[i][D->nySub+1].E1=field[i][1].E1;
    }
}

void periodY1coreC(Domain *D)
{
    int i,j,start,end,numberData;  
    float dx,dy,dt,preB1C;
    FieldElement **field;
    field=D->field;
   
    for(i=0; i<=D->nxSub+1; i++)
    {
       field[i][D->nySub+1].PrC=field[i][1].PrC;
       field[i][D->nySub+1].PlC=field[i][1].PlC;
       field[i][D->nySub+1].SrC=field[i][1].SrC;
       field[i][D->nySub+1].SlC=field[i][1].SlC;
       field[i][D->nySub+1].E1C=field[i][1].E1C;
       field[i][D->nySub+1].B1C=field[i][1].B1C;
    }
}
