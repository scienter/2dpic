#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>

void loadMovingPlasma2DBoost(Domain *D)
{
   int i,j,istart,iend,jstart,jend,s,l,intNum,cnt,np,nc,leftIndex,rightIndex;
   int posX,posY,iter,t;
   float tmp,dx,dy,cx,cy;
   float wp,pDt,v1,v2,v3,gamma,beta[2];

   float ne,randTest=0,positionX,positionY;
   double maxwellianVelocity();
   double randomValue();

   LoadList *LL;
   ptclList *New,*p;   
   Particle **particle;
   particle=D->particle;

   dx=D->dx;
   dy=D->dy;

   beta[0]=D->beta;
   beta[1]=1.0;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   //position define      
   iter=0;
 for(i=iend-1; i<=iend; i++)
 {
   for(j=jstart; j<jend; j++)
   {
     LL=D->loadList;
     s=0;
     while(LL->next)
     {
       cx=LL->cx;
       cy=LL->cy;
       for(l=0; l<LL->lnodes-1; l++)
         for(t=0; t<LL->tnodes-1; t++)
         {
           
           if(LL->type==Circle) {
             tmp=sqrt((cx-(i-istart+D->minXSub)*dx)*(cx-(i-istart+D->minXSub)*dx)+
                    (cy-(j-jstart+D->minYSub)*dy)*(cy-(j-jstart+D->minYSub)*dy));
             posX=(int)((cx-tmp)/dx+istart)-istart;
             posY=j+D->minYSub-jstart;
           }
           else if(LL->type==Point4) {
             posX=i+D->minXSub-istart;
             posY=j+D->minYSub-jstart;
           }
 
           if((posX>=LL->lpoint[l] && posX<LL->lpoint[l+1]) && 
              (posY>=LL->tpoint[t] && posY<LL->tpoint[t+1]))
           {
             ne=((LL->ln[l+1]-LL->ln[l])/(LL->lpoint[l+1]-LL->lpoint[l])
                *(posX-LL->lpoint[l])+LL->ln[l]);
             ne*=((LL->tn[t+1]-LL->tn[t])/(LL->tpoint[t+1]-LL->tpoint[t])
                *(posY-LL->tpoint[t])+LL->tn[t]);
             ne*=LL->numberInCell*beta[iter];	//it is the float number of superparticles.
             intNum=(int)ne;
             randTest=ne-intNum;
             cnt=0;
             while(cnt<intNum)
             {               
               positionX=randomValue(beta[iter]);
               positionY=randomValue(1.0);
  
               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
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
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;

                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            

                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s+1]->pt;
                 particle[i][j].head[s+1]->pt = New;

                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
               } 
             
               cnt++;
             }		//end of while(cnt)

             if(randTest>randomValue(1.0))
             {
               positionX=randomValue(beta[iter]);
               positionY=randomValue(1.0);

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
 
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
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
 
                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
 
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s+1]->pt;
                 particle[i][j].head[s+1]->pt = New;
 
                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
               } 
             }		//end of if(randTest)
           }		//end of if(l,t nodes)
         }		//end of for(lnodes,tnodes)  

       LL=LL->next;
       s++;
     }				//End of while(LL)   
   } 				//End of for(j)
   iter=iter+1;
 }				//End of for(i)
   //other define define      
   LL=D->loadList;
   s=0;
   while(LL->next)
   {
     for(i=iend-1; i<=iend; i++)
       for(j=jstart; j<jend; j++)
       {
         p=particle[i][j].head[s]->pt;
         while(p)
         {
           p->E1=p->B1=p->Pr=p->Pl=p->Sr=p->Sl=0.0;
           v1=maxwellianVelocity(LL->temperature)/velocityC;
           v2=maxwellianVelocity(LL->temperature)/velocityC;
           v3=maxwellianVelocity(LL->temperature)/velocityC;
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

void loadMovingPlasma2D(Domain *D)
{
   int i,j,istart,iend,jstart,jend,s,l,intNum,cnt,np,nc,leftIndex,rightIndex;
   int posX,posY,iter,t;
   float tmp,dx,dy,cx,cy;
   float wp,pDt,v1,v2,v3,gamma;

   float ne,randTest=0,positionX,positionY;
   double maxwellianVelocity();
   double randomValue();

   LoadList *LL;
   ptclList *New,*p;   
   Particle **particle;
   particle=D->particle;

   dx=D->dx;
   dy=D->dy;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   //position define      
   i=iend-1;
   for(j=jstart; j<jend; j++)
   {
     LL=D->loadList;
     s=0;
     while(LL->next)
     {
       cx=LL->cx;
       cy=LL->cy;
       for(l=0; l<LL->lnodes-1; l++)
         for(t=0; t<LL->tnodes-1; t++)
         {
           if(LL->type==Circle) {
             tmp=sqrt((cx-(i-istart+D->minXSub)*dx)*(cx-(i-istart+D->minXSub)*dx)+
                    (cy-(j-jstart+D->minYSub)*dy)*(cy-(j-jstart+D->minYSub)*dy));
             posX=(int)((cx-tmp)/dx+istart)-istart;
             posY=j+D->minYSub-jstart;
           }
           else if(LL->type==Point4)  {
             posX=i+D->minXSub-istart;
             posY=j+D->minYSub-jstart;
           }
 
           if(posX>=LL->lpoint[l] && posX<LL->lpoint[l+1] &&
              posY>=LL->tpoint[t] && posY<LL->tpoint[t+1])
           {
             ne=((LL->ln[l+1]-LL->ln[l])/(LL->lpoint[l+1]-LL->lpoint[l])
                *(posX-LL->lpoint[l])+LL->ln[l]);
             ne*=((LL->tn[t+1]-LL->tn[t])/(LL->tpoint[t+1]-LL->tpoint[t])
                *(posY-LL->tpoint[t])+LL->tn[t]);
             ne*=LL->numberInCell;	//it is the float number of superparticles.
             intNum=(int)ne;
             randTest=ne-intNum;
             cnt=0;
             while(cnt<intNum)
             {               
               positionX=randomValue(1.0);
               positionY=randomValue(1.0);

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
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
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;

                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
 
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s+1]->pt;
                 particle[i][j].head[s+1]->pt = New;
 
                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
               } 
             
               cnt++;
             }		//end of while(cnt)

             if(randTest>randomValue(1.0))
             {
               positionX=randomValue(1.0);
               positionY=randomValue(1.0);

               if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
               {
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
 
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
                 New->next = particle[i][j].head[s]->pt;
                 particle[i][j].head[s]->pt = New;
 
                 New->x = positionX;
                 New->oldX=i+positionX;
                 New->y = positionY;
                 New->oldY=j+positionY;
                 LL->index++;
                 New->index=LL->index;            
 
                 New = (ptclList *)malloc(sizeof(ptclList)); 
                 New->next = particle[i][j].head[s+1]->pt;
                 particle[i][j].head[s+1]->pt = New;
 
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
       for(j=jstart; j<jend; j++)
       {
         p=particle[i][j].head[s]->pt;
         while(p)
         {
           p->E1=p->B1=p->Pr=p->Pl=p->Sr=p->Sl=0.0;
           v1=maxwellianVelocity(LL->temperature)/velocityC;
           v2=maxwellianVelocity(LL->temperature)/velocityC;
           v3=maxwellianVelocity(LL->temperature)/velocityC;
           gamma=sqrt(1.0+(v1*v1+v2*v2+v3*v3));
           p->p1=-D->gamma*D->beta+v1;
           p->p2=v2;
           p->p3=v3;

           p=p->next;
         } 
       }				//end of for(j)
       LL=LL->next;
       s++;
   }					//End of while(LL)

}

void loadPlasma2D(Domain *D)
{
   int i,j,istart,iend,jstart,jend;
   int s,l,intNum,cnt,np,nc,leftIndex,rightIndex;
   int posX,posY,t;
   float cx,cy,tmp,dx,dy;
   float wp,pDt,v1,v2,v3,gamma;
   float ne,randTest,positionX,positionY;
   double maxwellianVelocity();
   double randomValue();
   Particle **particle;
   particle=D->particle;

   LoadList *LL;
//   FieldDSX **field;
//   field=D->fieldDSX;
   ptclList *New,*p;   

   dx=D->dx;
   dy=D->dy;

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   //position define      
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       LL=D->loadList;
       s=0;
       while(LL->next)
       {
         cx=LL->cx;
         cy=LL->cy;
         for(l=0; l<LL->lnodes-1; l++)
           for(t=0; t<LL->tnodes-1; t++)
           {
             if(LL->type==Circle) {
               tmp=sqrt((cx-(i-istart+D->minXSub)*dx)*(cx-(i-istart+D->minXSub)*dx)+
                      (cy-(j-jstart+D->minYSub)*dy)*(cy-(j-jstart+D->minYSub)*dy));
               posX=(int)((cx-tmp)/dx+istart)-istart;
               posY=j+D->minYSub-jstart;
             }
             else if(LL->type==Point4) {
               posX=i+D->minXSub-istart;
               posY=j+D->minYSub-jstart;
             }
 
             if(posX>=LL->lpoint[l] && posX<LL->lpoint[l+1] &&
                posY>=LL->tpoint[t] && posY<LL->tpoint[t+1])
             {
               ne=((LL->ln[l+1]-LL->ln[l])/(LL->lpoint[l+1]-LL->lpoint[l])*(posX-LL->lpoint[l])+LL->ln[l]);
               ne*=((LL->tn[t+1]-LL->tn[t])/(LL->tpoint[t+1]-LL->tpoint[t])*(posY-LL->tpoint[t])+LL->tn[t]);
               ne*=LL->numberInCell;	//it is the float number of superparticles.
               intNum=(int)ne;
               randTest=ne-intNum;
             
               cnt=0;
               while(cnt<intNum)
               {               
                 positionX=randomValue(1.0);
                 positionY=randomValue(1.0);

                 if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
                 {
                   New = (ptclList *)malloc(sizeof(ptclList)); 
                   New->next = particle[i][j].head[s]->pt;
                   particle[i][j].head[s]->pt = New;
 
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
                   New->next = particle[i][j].head[s]->pt;
                   particle[i][j].head[s]->pt = New;
  
                   New->x = positionX;
                   New->oldX=i+positionX;
                   New->y = positionY;
                   New->oldY=j+positionY;
                   LL->index++;
                   New->index=LL->index;            
 
                   New = (ptclList *)malloc(sizeof(ptclList)); 
                   New->next = particle[i][j].head[s+1]->pt;
                   particle[i][j].head[s+1]->pt = New;
 
                   New->x = positionX;
                   New->oldX=i+positionX;
                   New->y = positionY;
                   New->oldY=j+positionY;
                   LL->index++;
                   New->index=LL->index;            
                 } 
               
                 cnt++;
               }		//end of while(cnt)

               if(randTest>randomValue(1.0))
               {
                 positionX=randomValue(1.0);
                 positionY=randomValue(1.0);

                 if(LL->withNextSpcs==0 && LL->withPrevSpcs==0)
                 {
                   New = (ptclList *)malloc(sizeof(ptclList)); 
                   New->next = particle[i][j].head[s]->pt;
                   particle[i][j].head[s]->pt = New;
 
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
                   New->next = particle[i][j].head[s]->pt;
                   particle[i][j].head[s]->pt = New;
  
                   New->x = positionX;
                   New->oldX=i+positionX;
                   New->y = positionY;
                   New->oldY=j+positionY;
                   LL->index++;
                   New->index=LL->index;            
 
                   New = (ptclList *)malloc(sizeof(ptclList)); 
                   New->next = particle[i][j].head[s+1]->pt;
                   particle[i][j].head[s+1]->pt = New;
 
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

     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++)
       {
         p=particle[i][j].head[s]->pt;
         while(p)
         {
           p->E1=p->B1=p->Pr=p->Pl=p->Sr=p->Sl=0.0;
           v1=maxwellianVelocity(LL->temperature)/velocityC;
           v2=maxwellianVelocity(LL->temperature)/velocityC;
           v3=maxwellianVelocity(LL->temperature)/velocityC;
//           gamma=1.0/sqrt(1.0-(v1*v1+v2*v2+v3*v3));
           p->p1=-D->gamma*D->beta+v1;
           p->p2=v2;
           p->p3=v3;

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

double randomValue(double beta)
{
   double r;
   int intRand, randRange=100, rangeDev;

   rangeDev=(int)(randRange*(1.0-beta));
   intRand = rand() % (randRange-rangeDev);
   r = ((float)intRand)/randRange+(1.0-beta);

   return r;
}
