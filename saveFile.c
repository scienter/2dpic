#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void saveProbe(Domain *D,int iteration)
{
    int i,j,n;
    char name[100];
    float t,Ex,Ey,Ez,Bx,By,Bz,Pr,Pl,Sr,Sl,x,y;
    float omega,frequency,dt;
    FILE *out;
    int myrank, nprocs;    
    Probe **probe;
    probe=D->probe;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    omega=2*pi*velocityC/D->lambda;
    frequency=omega/2.0/pi;   
    dt=1.0/frequency*D->dt;

    for(n=0; n<D->probeNum; n++)
    {
      if(D->probeX[n]>=D->minXSub && D->probeX[n]<D->maxXSub && 
         D->probeY[n]>=D->minYSub && D->probeY[n]<D->maxYSub)
      {
        x=D->probeX[n]*D->dx*D->lambda;
        y=D->probeY[n]*D->dy*D->lambda;
        sprintf(name,"probeRaman%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          Pr=probe[n][i].Pr;
          Pl=probe[n][i].Pl;
          Sr=probe[n][i].Sr;
          Sl=probe[n][i].Sl;
          fprintf(out,"%g %g %g %g %g %g %g\n",t,Pr,Pl,Sr,Sl,x,y);
        }
        fclose(out);

        sprintf(name,"probe%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          Ex=probe[n][i].E1;
          Bx=probe[n][i].B1;
          Ey=probe[n][i].Pr+probe[n][i].Pl;
          Ez=probe[n][i].Sr+probe[n][i].Sl;
          By=probe[n][i].Sl-probe[n][i].Sr;
          Bz=probe[n][i].Pr-probe[n][i].Pl;
          fprintf(out,"%g %g %g %g %g %g %g %g %g\n",t,Ex,Ey,Ez,Bx,By,Bz,x,y);
        }             
        fclose(out);
      }
    }
}

void saveRho2D(Domain *D,int iteration)
{
    int i,j,istart,iend,jstart,jend,s;
    float x,y,a00,a11,a01,a10;
    char name[100];
    Particle **particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    float rho0[D->nSpecies];
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       rho0[s]=LL->charge*LL->density/((float)LL->numberInCell);
       LL=LL->next;
       s++;
    }

    //initializing density
    for(i=istart; i<=iend; i++)
      for(j=jstart; j<=jend; j++)
        particle[i][j].rho=0.0;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
//        for(s=0; s<D->nSpecies; s++)
        for(s=0; s<1; s++)
        {
          p=particle[i][j].head[s]->pt;
          while(p)
          {
            x=p->x;
            y=p->y;
            a00=(1-x)*(1-y);
            a11=x*y;
            a10=(1-x)*y;
            a01=x*(1-y);
            particle[i][j].rho+=a00*rho0[s];
            particle[i+1][j+1].rho+=a11*rho0[s];
            particle[i+1][j].rho+=a10*rho0[s];
            particle[i][j+1].rho+=a01*rho0[s];
            p=p->next;
          }
        }

    sprintf(name,"rho%d_%d",iteration,myrank);
    out = fopen(name,"w");    
    for(i=istart; i<iend; i++)
    {
      for(j=jstart; j<jend; j++)
      {
        x=(i-istart+D->minXSub)*D->dx*D->lambda;
        y=(j-jstart+D->minYSub)*D->dy*D->lambda;
          fprintf(out,"%g %g %g\n",x,y,particle[i][j].rho);    
      }           
      fprintf(out,"\n");    
    }
    fclose(out);

}


void boostSaveField(Domain *D,int labSaveStep)
{
    int i,j,istart,iend,jstart,jend,show;
    char name[100];
    float x,y,e1,pr,pl,b1,sr,sl;
    float factor,dx;
    FILE *out;
    int myrank, nprocs;    
    Boost **boost;
    boost=D->boost;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    sprintf(name,"bField%d_%d",labSaveStep,myrank);
    out = fopen(name,"w");

    for(i=istart; i<iend; i++)
    {
      for(j=jstart; j<jend; j++)
      {
        x=boost[i][j].x;
        y=boost[i][j].y;
        e1=boost[i][j].E1;
        pr=boost[i][j].Pr;
        pl=boost[i][j].Pl;
        b1=boost[i][j].B1;
        sr=boost[i][j].Sr;
        sl=boost[i][j].Sl;
        if(x>0) {
          show=1;
          fprintf(out,"%g %g %g %g %g %g %g %g\n",x,y,e1,pr,pl,b1,sr,sl);
        }
        else show=0;
      }  
      if(show==1)           
        fprintf(out,"\n");
    }
    fclose(out);
    
    if(myrank==0)
      printf("bField%d is saved.\n",labSaveStep);
}

void saveField2D(Domain *D,int iteration)
{
    int i,j,istart,iend,jstart,jend;
    char name[100];
    float x,y,Ex,Ey,Ez,Bx,By,Bz,factor;
    FILE *out;
    int myrank, nprocs;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    if(D->fieldType==1)
    {
      FieldDSX **field;
      field=D->fieldDSX;
      sprintf(name,"field%d_%d",iteration,myrank);
      out = fopen(name,"w");
      factor=D->gamma*(1+D->beta);
    
      for(i=istart; i<iend; i++)
      {
        for(j=jstart; j<jend; j++)
        {
          x=(i-2+D->minXSub+0.5)*D->dx*D->lambda;
          y=(j-2+D->minYSub+0.5)*D->dy*D->lambda;
          Ex=field[i][j].E1;    
          Ey=field[i][j].Pr+field[i][j].Pl;
          Ez=field[i][j].Sr+field[i][j].Sl;
          Bx=field[i][j].B1;    
          By=field[i][j].Sl-field[i][j].Sr;
          Bz=field[i][j].Pr-field[i][j].Pl;      
          fprintf(out,"%g %g %g %g %g %g %g %g\n",x,y,Ex,Ey,Ez,Bx,By,Bz);
        }
        fprintf(out,"\n");                 
      }
      fclose(out);
    }
/*
    else if(D->fieldType==12)  {
       FieldDSXY **field;
       field=D->fieldDSXY;
    }
*/

}

void saveRaman2D(Domain *D,int iteration)
{
    int i,j,istart,iend,jstart,jend;
    char name[100];
    float x,y,Pr,Pl,Pu,Pd,Sr,Sl,Su,Sd;
    FILE *out;
    int myrank, nprocs;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    
    if(D->fieldType==1)
    {
      FieldDSX **field;
      field=D->fieldDSX;

      sprintf(name,"raman%d_%d",iteration,myrank);
      out = fopen(name,"w");

      for(i=istart; i<iend; i++)
      {
        for(j=jstart; j<jend; j++)
        {
          x=(i-2+D->minXSub)*D->dx*D->lambda;
          y=(j-2+D->minYSub)*D->dy*D->lambda;
          Pr=field[i][j].Pr;
          Pl=field[i][j].Pl;
          Sr=field[i][j].Sr;
          Sl=field[i][j].Sl;
          fprintf(out,"%g %g %g %g %g %g\n",x,y,Pr,Pl,Sr,Sl);
        }
        fprintf(out,"\n");
      }             
      fclose(out);
    }
}

void saveParticle(Domain *D,int iteration)
{
    int i,j,istart,iend,jstart,jend,s;
    char name[100];
    float x,y,p1,p2,p3,index,gamma,mc;
    float Pr,Pl,E1,Sr,Sl,B1;
    Particle **particle;
    particle=D->particle;
    ptclList *p;
    FILE *out;
    int myrank, nprocs;    

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    for(s=0; s<D->nSpecies; s++)
    {
      sprintf(name,"%dParticle%d_%d",s,iteration,myrank);
      out = fopen(name,"w");    
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
        {
          p=particle[i][j].head[s]->pt;
          while(p)
          {
             x=((i-istart+D->minXSub)+p->x)*D->dx*D->lambda; 
             y=((j-jstart+D->minYSub)+p->y)*D->dy*D->lambda; 
             p1=p->p1;
             p2=p->p2;    
             p3=p->p3;
             gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
             index=p->index;
             fprintf(out,"%g %g %g %g %g %g %g\n",x,y,p1,p2,p3,gamma,index);               
             p=p->next;
          }
        }		//End of for(i,j)
      fclose(out);
    }				//End of for(s)
}

