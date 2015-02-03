#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include "plasma.h"


void boostShot(Domain *D,int iteration)
{
    int i,j,s,labStep,istart,iend,jstart,jend,time;
    char name[100];
    FILE *out;
    float x,y,e1,e2,e3,b1,b2,b3;
    float xx,yy,E1,E2,E3,B1,B2,B3;
    float tmp,factor,factor1;
    int myrank, nprocs, cnt;    
    float rho,density,J1;	//save density
    Particle **particle;
    particle=D->particle;
    ptclList *p;
    float p1,p2,p3,pp1,gamma,index;

    LoadList *LL;

    LL=D->loadList;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    factor=D->gamma*(1+D->beta);
    factor1=1+D->beta;

    if(D->fieldType==1)
    {   
      FieldDSX **field;
      field=D->fieldDSX;
//      boost=D->boost;

      for(time=D->saveStart; time<=D->maxStep; time+=D->saveStep)
      {
        i=(int)(time/D->gamma/D->beta/factor-iteration/D->beta-D->minXSub+istart);
        if(i>=istart && i<iend)
        {
          sprintf(name,"bField%d_%d_%d",time,i-istart,myrank);
          remove(name);  
          out = fopen(name,"w");
          for(j=jstart; j<jend; j++) 
          {       
            //boost frame data
            x=(i-istart+D->minXSub)*D->dx*D->lambda;
            y=(j-jstart+D->minYSub)*D->dy*D->lambda;
            e1=field[i][j].E1;
            e2=field[i][j].Pr+field[i][j].Pl;
            e3=field[i][j].Sr+field[i][j].Sl;
            b1=field[i][j].B1;
            b2=field[i][j].Sl-field[i][j].Sr;
            b3=field[i][j].Pr-field[i][j].Pl;

            //lab frame data
            xx=D->gamma*(x+D->beta*iteration*D->dt*D->lambda);
            yy=y;
            E1=e1/factor;
            E2=(e2+D->beta*b3)/factor1;
            E3=(e3-D->beta*b2)/factor1;
            B1=b1/factor;
            B2=(b2-D->beta*b3)/factor1;
            B3=(b3+D->beta*e2)/factor1;
  
            fprintf(out,"%g %g %g %g %g %g %g %g\n",xx,yy,E1,E2,E3,B1,B2,B3);
          }	//End of for(j)
          fclose(out);

        }	//End of if(i)
      }		//End of for(time)
    }			//End of fieldType=1

    //save particles
    for(time=D->saveStart; time<=D->maxStep; time+=D->saveStep)
    {
      i=(int)(time/D->gamma/D->beta/factor-iteration/D->beta-D->minXSub+istart);
      sprintf(name,"bParticle%d_%d",time,myrank);        

      if(i>=istart && i<iend)
      {
        out = fopen(name,"a");
        for(j=jstart; j<jend; j++) 
        {       
          p=particle[i][j].head[0]->pt;
          while(p)
          {
            x=(i-istart+D->minXSub+p->x)*D->dx*D->lambda;
            y=(j-jstart+D->minYSub+p->y)*D->dy*D->lambda;
            p1=p->p1;
            p2=p->p2;
            p3=p->p3;
            gamma=sqrt(1.0+p1*p1+p2*p2+p3*p3);
            index=p->index;
 
            xx=D->gamma*(x+D->beta*iteration*D->dt*D->lambda);
            yy=y;
            pp1=D->gamma*(p1+D->beta*gamma);
            gamma=sqrt(1.0+pp1*pp1+p2*p2+p3*p3);
  
            fprintf(out,"%g %g %g %g %g %g %g\n",xx,yy,pp1,p2,p3,gamma,index);
            p=p->next;
          }
        }	//End of for(j)
        fclose(out);
      }		//End of if(i)
    }			//End of for(time)      

    //save density
    float rho0[D->nSpecies];
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       rho0[s]=LL->charge*LL->density/((float)LL->numberInCell);
       LL=LL->next;
       s++;
    }

    for(time=D->saveStart; time<=D->maxStep; time+=D->saveStep)
    {
      i=(int)(time/D->gamma/D->beta/factor-iteration/D->beta-D->minXSub+istart);

      if(i>=istart && i<iend)
      {
        sprintf(name,"bDensity%d_%d_%d",time,i,myrank);
        remove(name);
        out = fopen(name,"w");
        for(j=jstart; j<jend; j++) 
        {
          density=0;
          for(s=0; s<1; s++) 
          {       
            p=particle[i][j].head[s]->pt;
            cnt=0;
            while(p)
            {
              cnt++;
              p=p->next;
            }
            density+=cnt*rho0[s];
          }
          x=(i-istart+D->minXSub)*D->dx*D->lambda;
          y=(j-jstart+D->minYSub)*D->dy*D->lambda;
          xx=D->gamma*(x+D->beta*iteration*D->dt*D->lambda);
          yy=y;
          rho=D->gamma*(density+D->beta*D->fieldDSX[i][j].J1*LL->criticalDensity);
  
          fprintf(out,"%g %g %g\n",xx,yy,rho);
        }	//End of for(j)
        fclose(out);
      }		//End of if(i)
    }			//End of for(time)      

}

