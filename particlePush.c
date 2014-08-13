#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "plasma.h"
#include <math.h>

void particlePush2D(Domain *D)
{
    int i,j,istart,iend,jstart,jend,l,m,s,shift,cnt;
    float x,shiftX,shiftY,dt,dx,dy,gamma,sqrT,coef,dxOverdy;
    float E1,Pr,Pl,B1,Sr,Sl;
    float pMinus[3],T[3],S[3],operate[3][3],pPlus[3];
    Particle **particle;
    particle=D->particle;
    LoadList *LL;
    ptclList *p, *New, *tmp, *prev;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    for(i=0; i<3; i++)   {
       pMinus[i]=0.0;
       T[i]=0.0;
       S[i]=0.0;
       pPlus[i]=0.0;
    }

    dt=D->dt;
    dx=D->dx;
    dy=D->dy;
    dxOverdy=D->dx/D->dy;

    float mass[D->nSpecies];
    int charge[D->nSpecies];
    LL=D->loadList;
    s=0;

    while(LL->next)
    {
       mass[s]=LL->mass;
       charge[s]=LL->charge;
       s++;
       LL=LL->next;
    }

    shiftX=shiftY=0;

    for(i=istart; i<=iend; i++)
      for(j=jstart; j<jend; j++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          coef=pi*charge[s]/mass[s]*dt;
          p=particle[i][j].head[s]->pt;     
          while(p)
          {    
             //Calculate vector P- 
             pMinus[0]=p->p1+coef*(p->E1);
             pMinus[1]=p->p2+coef*(p->Pr+p->Pl);    
             pMinus[2]=p->p3+coef*(p->Sr+p->Sl);
 
             //Calculate vector T 
             gamma=sqrt(1.0+pMinus[0]*pMinus[0]+pMinus[1]*pMinus[1]+pMinus[2]*pMinus[2]);           
             T[0]=coef/gamma*(p->B1);   
             T[1]=coef/gamma*(p->Sl-p->Sr);
             T[2]=coef/gamma*(p->Pr-p->Pl);

             //Calculate vector S
             sqrT=1.0+T[0]*T[0]+T[1]*T[1]+T[2]*T[2];
             for(l=0; l<3; l++)  
                S[l]=2.0*T[l]/sqrT;
  
             //Calculate operator A from P+=A.P-
             operate[0][0]=1.0-S[2]*T[2]-S[1]*T[1];      
             operate[0][1]=S[1]*T[0]+S[2];    
             operate[0][2]=S[2]*T[0]-S[1];     
             operate[1][0]=S[0]*T[1]-S[2];        
             operate[1][1]=1.0-S[0]*T[0]-S[2]*T[2];          
             operate[1][2]=S[2]*T[1]+S[0];         
             operate[2][0]=S[0]*T[2]+S[1];    
             operate[2][1]=S[1]*T[2]-S[0];    
             operate[2][2]=1-S[0]*T[0]-S[1]*T[1]; 
             //Calculate vector P+
             for(l=0; l<3; l++)  {
                pPlus[l]=0.0;
                for(m=0; m<3; m++)   
                   pPlus[l]+=operate[l][m]*pMinus[m];
                }
             //Updated momentum              
             p->p1=pPlus[0]+coef*(p->E1); 
             p->p2=pPlus[1]+coef*(p->Pr+p->Pl);    
             p->p3=pPlus[2]+coef*(p->Sr+p->Sl); 
    
             //Translation
             gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
             shiftX=p->p1/gamma;    //dt is ignored because of dx=dt=1 in cell.
             shiftY=p->p2/gamma*dxOverdy;    //dt is ignored because of dx=dt=1 in cell.
             if(shiftX>=1 || shiftY>=1*dxOverdy)  {
                printf("particle's movement exceeds C velocity\n");
                printf("i=%d, j=%d,shiftX=%g, shiftY=%g\n",i,j,shiftX,shiftY);
                exit(0);
             } 
             p->oldX=i+p->x;
             p->x+=shiftX;
             p->oldY=j+p->y;
             p->y+=shiftY;
             p=p->next;
          }		//End of while(p)
        }		//End of for(s)
      }      	//End of for(i,j)

}

             
