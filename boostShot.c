#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include "plasma.h"


void boostShot(Domain *D,int iteration)
{
    int i,j,s,labStep,istart,iend,jstart,jend,time,shiftIndex;
    char name[100];
    FILE *out;
    float x,y,e1,e2,e3,b1,b2,b3;
    float xx,yy,E1,E2,E3,B1,B2,B3;
    float tmp,factor,factor1,weight,current;
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
        tmp=time/D->gamma/D->beta/factor-iteration/D->beta;
        weight=tmp-(int)tmp;
        if(i>=istart && i<iend)
        {
          sprintf(name,"bField%d_%d_%d",time,i-istart,myrank);
          remove(name);  
          out = fopen(name,"w");
          x=tmp*D->dx*D->lambda;
//          x=(i-istart+D->minXSub)*D->dx*D->lambda;
          xx=D->gamma*(x+D->beta*iteration*D->dt*D->lambda);
          for(j=jstart; j<jend; j++) 
          {       
            //boost frame data
            y=(j-jstart+D->minYSub)*D->dy*D->lambda;
            e1=(1-weight)*field[i][j].E1+weight*field[i+1][j].E1;
            e2=(1-weight)*field[i][j].Pr+field[i][j].Pl
              +(weight)*(field[i+1][j].Pr+field[i+1][j].Pl);
            e3=(1-weight)*(field[i][j].Sr+field[i][j].Sl)
              +(weight)*(field[i+1][j].Sr+field[i+1][j].Sl);
            b1=(1-weight)*field[i][j].B1+weight*field[i+1][j].B1;
            b2=(1-weight)*(field[i][j].Sl-field[i][j].Sr)
              +(weight)*(field[i+1][j].Sl-field[i+1][j].Sr);
            b3=(1-weight)*(field[i][j].Pr-field[i][j].Pl)
              +(weight)*(field[i+1][j].Pr-field[i+1][j].Pl);

            //lab frame data
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
      tmp=time/D->gamma/D->beta/factor-iteration/D->beta;
      weight=tmp-(int)tmp;
      if(i>=istart && i<iend)
      {
        shiftIndex=(int)(weight+0.5)-1;
        sprintf(name,"bDensity%d_%d_%d",time,i-istart,myrank);
        remove(name);
        out = fopen(name,"w");
        x=tmp*D->dx*D->lambda;
        xx=D->gamma*(x+D->beta*iteration*D->dt*D->lambda);
        for(j=jstart; j<jend; j++) 
        {
          density=0;
          for(s=0; s<1; s++) 
          {       
            cnt=0;
            p=particle[i+shiftIndex][j].head[s]->pt;
            while(p)
            {
              if(p->x>0.5+weight-(int)(weight+0.5))
                cnt++;
              p=p->next;
            }
            p=particle[i+shiftIndex+1][j].head[s]->pt;
            while(p)
            {
              if(p->x<0.5+weight-(int)(weight+0.5))
                cnt++;
              p=p->next;
            }
            density+=cnt*rho0[s];
          }

          y=(j-jstart+D->minYSub)*D->dy*D->lambda;
          yy=y;
          weight=weight+0.5-(int)(weight+0.5);
          current=(1-weight)*D->fieldDSX[i+shiftIndex-1][j].J1+weight*D->fieldDSX[i+shiftIndex][j].J1;
          rho=D->gamma*(density+D->beta*current*LL->criticalDensity);
  
          fprintf(out,"%g %g %g\n",xx,yy,rho);
        }	//End of for(j)
        fclose(out);
      }		//End of if(i)
    }			//End of for(time)      

}

