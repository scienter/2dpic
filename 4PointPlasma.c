#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>

void loadMovingPlasma2D(Domain *D)
{
   int i,j,s,l,intNum,cnt,np,nc,leftIndex,rightIndex;
   float wp,pDt,v1,v2,v3,gamma;

   float ne,randTest=0,positionX,positionY;
   double maxwellianVelocity();
   double randomValue();

   LoadList *LL;
   FieldElement **field;
   field=D->field;
   ptclList *New,*p;   

   i=D->nxSub+1;

   //position define      
     for(j=1; j<=D->nySub; j++)
     {
       LL=D->loadList;
       s=0;
       while(LL->next)
       {
         for(l=0; l<LL->lnodes-1; l++)
         {
           if(i+D->minXSub>=LL->lpoint[l] && i+D->minXSub<LL->lpoint[l+1])
           {
             ne=((LL->ln[l+1]-LL->ln[l])/(LL->lpoint[l+1]-LL->lpoint[l])
                *(i+D->minXSub-LL->lpoint[l])+LL->ln[l]);
             ne*=LL->numberInCell;	//it is the float number of superparticles.
             intNum=(int)ne;
             randTest=ne-intNum;
             cnt=0;
             while(cnt<intNum)
             {               
               positionX=randomValue();
               positionY=randomValue();

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s]->pt;
                  field[i][j].head[s]->pt = New;

                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            
               } 
               else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s]->pt;
                  field[i][j].head[s]->pt = New;
 
                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            

                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s+1]->pt;
                  field[i][j].head[s+1]->pt = New;

                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            
               } 
               
               cnt++;
             }		//end of while(cnt)

             if(randTest>randomValue())
             {
               positionX=randomValue();
               positionY=randomValue();

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s]->pt;
                  field[i][j].head[s]->pt = New;

                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            
               } 
               else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s]->pt;
                  field[i][j].head[s]->pt = New;
 
                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            

                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s+1]->pt;
                  field[i][j].head[s+1]->pt = New;

                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            
               } 
             }		//end of if(randTest)

           }
         } 			//end of for(lnodes)  

         LL=LL->next;
         s++;
       }				//End of while(LL)   
     } 				//End of for(j)


   //other define define      
   LL=D->loadList;
   s=0;
   while(LL->next)
   {
       for(j=1; j<=D->nySub; j++)
       {
         p=field[i][j].head[s]->pt;
         while(p)
         {
           p->E1=p->B1=p->Pr=p->Pl=p->Sr=p->Sl=0.0;
           v1=0; //maxwellianVelocity(LL->temperature)/velocityC;
           v2=0; //maxwellianVelocity(LL->temperature)/velocityC;
           v3=0; //maxwellianVelocity(LL->temperature)/velocityC;
           gamma=1.0/sqrt(1.0-(v1*v1+v2*v2+v3*v3));
           p->p1=-D->gamma*D->beta+gamma*v1;
           p->p2=gamma*v2;
           p->p3=gamma*v3;

           p=p->next;
         } 
       }				//end of for(j)
       LL=LL->next;
       s++;
   }					//End of while(LL)

}


void loadPlasma2D(Domain *D)
{
   int i,j,s,l,intNum,cnt,np,nc,leftIndex,rightIndex;
   float wp,pDt,v1,v2,v3,gamma;
   float ne,randTest,positionX,positionY;
   double maxwellianVelocity();
   double randomValue();

   LoadList *LL;
   FieldElement **field;
   field=D->field;
   ptclList *New,*p;   

   for(i=0; i<D->nxSub+2; i++)
     for(j=0; j<D->nySub+2; j++)
     {
        field[i][j].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
        for(s=0; s<D->nSpecies; s++)
        {
           field[i][j].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
           field[i][j].head[s]->pt = NULL;
        }
     }
   
   //position define      
   for(i=1; i<=D->nxSub; i++)
     for(j=1; j<=D->nySub; j++)
     {
       LL=D->loadList;
       s=0;
       while(LL->next)
       {
         for(l=0; l<LL->lnodes-1; l++)
         {
           if(i+D->minXSub>=LL->lpoint[l] && i+D->minXSub<LL->lpoint[l+1])
           {
             ne=((LL->ln[l+1]-LL->ln[l])/(LL->lpoint[l+1]-LL->lpoint[l])
                *(i+D->minXSub-LL->lpoint[l])+LL->ln[l]);
             ne*=LL->numberInCell;	//it is the float number of superparticles.
             intNum=(int)ne;
             randTest=ne-intNum;
             
             cnt=0;
             while(cnt<intNum)
             {               
               positionX=randomValue();
               positionY=randomValue();

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s]->pt;
                  field[i][j].head[s]->pt = New;

                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            
               } 
               else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s]->pt;
                  field[i][j].head[s]->pt = New;
 
                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            

                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s+1]->pt;
                  field[i][j].head[s+1]->pt = New;

                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            
               } 
               
               cnt++;
             }		//end of while(cnt)

             if(randTest>randomValue())
             {
               positionX=randomValue();
               positionY=randomValue();

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s]->pt;
                  field[i][j].head[s]->pt = New;

                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            
               } 
               else if(LL->withNextSpcs==1 && LL->withPrevSpcs==0)
               {
                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s]->pt;
                  field[i][j].head[s]->pt = New;
 
                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            

                  New = (ptclList *)malloc(sizeof(ptclList)); 
                  New->next = field[i][j].head[s+1]->pt;
                  field[i][j].head[s+1]->pt = New;

                  New->x = positionX;
                  New->oldX=i+positionX;
                  New->y = positionY;
                  New->oldY=j+positionY;
                  LL->index++;
                  New->index=LL->index;            
               } 
             }		//end of if(randTest)
           }
         } 			//end of for(lnodes)  

         LL=LL->next;
         s++;
       }				//End of while(LL)   
     }				//End of for(i,j)
    
   //other define define      
   LL=D->loadList;
   s=0;
   while(LL->next)
   {
     LL->criticalDensity=eps0*eMass*D->omega*D->omega/eCharge/eCharge;
     LL->superP=LL->density*D->lambda*D->dx*D->lambda*D->dy/LL->numberInCell;

     for(i=1; i<=D->nxSub; i++)
       for(j=1; j<=D->nySub; j++)
       {
         p=field[i][j].head[s]->pt;
         while(p)
         {
           p->E1=p->B1=p->Pr=p->Pl=p->Sr=p->Sl=0.0;
           v1=0; //maxwellianVelocity(LL->temperature)/velocityC;
           v2=0; //maxwellianVelocity(LL->temperature)/velocityC;
           v3=0; //maxwellianVelocity(LL->temperature)/velocityC;
           gamma=1.0/sqrt(1.0-(v1*v1+v2*v2+v3*v3));
           p->p1=-D->gamma*D->beta+gamma*v1;
           p->p2=gamma*v2;
           p->p3=gamma*v3;

           p=p->next;
         } 
       }				//end of for(i,j)
       LL=LL->next;
       s++;
   }					//End of while(LL)
}

double maxwellianVelocity(double temperature)
{
   float vth,r,prob,v,random;
   int intRand,randRange=1e5;

   vth=sqrt(2.0*eCharge*temperature/eMass);
   
   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((float)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((float)intRand)/randRange;
      v = 6.0*(random-0.5);
      prob=exp(-v*v);
   }
   return vth*v;
}

double randomValue()
{
   double r;
   int intRand, randRange=100;

   intRand = rand() % randRange;
   r = ((float)intRand)/randRange;

   return r;
}
