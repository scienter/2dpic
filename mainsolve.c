#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "mesh.h"
#include "laser.h"
#include <math.h>
#include "mpi.h"

float source(float t,float f,float rU);
void boundary(Domain *D);
void yeeSolveHfield(Domain *D);
void yeeSolveEfield(Domain *D);
void laserLoad(Domain *D, Laser *L, float t);
void pointLaserLoad(Domain *D, Laser *L, float t);

int main(int argc, char *argv[])
{
    int myrank, nTasks;
    int i,j,k=0,iteration=0;
    int maxStep=200,saveStep=20;
    float dx,dy,ds,dt,t,rU,divisionLamda;
    FILE *out;
    char name[100];
    Domain D;
    Laser L;
    
    //MPI setting
    D.L=2;
    D.M=1;
    D.N=1;
    
    //setting parameters        
    D.nx=201+2;
    D.ny=201+2;   
    D.nz=1;
    L.divisionLambda=10; 
    L.lambda=1e-6;  
    L.spotsize=2e-6; 
     
    dx=L.lambda/L.divisionLambda;
    ds=dy=dx;
    dt=0.99*1.0/sqrt(1.0/dx/dx+1.0/dy/dy)/velocityC;    
    
    D.dt=dt;
    D.ds=ds;
    D.dimension=2;
    L.polarity=3;
    L.load=4;

    //---------------- MPI initiation --------------------------------------
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    //----------------------------------------------------------------------    
  
    //create mesh;   
    boundary(&D);

    //rooping time
    while(iteration<=maxStep)
    {                           
       // calculation H field
//       yeeSolveHfield(&D);       
//       yeeSolveEfield(&D);
       
       // laser load
//       laserLoad(&D,&L,t);
//       pointLaserLoad(&D,&L,t);       
    
       //save File      
/*       if(iteration%saveStep==0)   
       {
          sprintf(name,"field%d",iteration);
          out = fopen(name,"w");
          for(i=0; i<D.nx; i++)
          {
             for(j=0; j<D.ny; j++)
             {
                fprintf(out,"%d %d %g %g %g %g %g %g\n",i,j,
                        D.field[i][j][k].E1,D.field[i][j][k].E2,D.field[i][j][k].E3,
                        D.field[i][j][k].H1,D.field[i][j][k].H2,D.field[i][j][k].H3);
             }
             fprintf(out,"\n");
          }             
          fclose(out);
          printf("field%d is made.\n",iteration);
       }
*/       
       t=dt*iteration;        
       iteration+=1;      
    }     //end of time roop                                    
   MPI_Finalize();
   
   free(D.field);



}


