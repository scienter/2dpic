#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void interpolation2D_3rd(Domain *D,External *Ext)
{
   int i,j,i1,j1,istart,iend,jstart,jend,s,cnt;
   float E1,Pr,Pl,B1,Sr,Sl,extE1,extB1,extPr,extPl,extSr,extSl,x,y;
   float Wx1,Wx2,Wx3,Wx4,Wy1,Wy2,Wy3,Wy4;
   float a,x1,x2,x3,x4,y1,y2,y3,y4;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   Particle **particle;
   particle=D->particle;

   extE1=Ext->E1;
   extB1=Ext->B1;
   extPr=Ext->Pr;
   extPl=Ext->Pl;
   extSr=Ext->Sr;
   extSl=Ext->Sl;
   if(D->fieldType==1)
   {   
     FieldDSX **field;
     field=D->fieldDSX;

     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++)
       {
         for(s=0; s<D->nSpecies; s++)
         {
           cnt=0;
           p=particle[i][j].head[s]->pt;
           while(p)
           {
             x=p->x+i;          y=p->y+j;
             i1=(int)(x+0.5);   j1=(int)(y+0.5);
             x=x-i1;            y=y-j1;
             x1=1.5+x;
             x2=0.5+x;
             x3=0.5-x;
             x4=1.5-x;
             y1=1.5+y;
             y2=0.5+y;
             y3=0.5-y;
             y4=1.5-y;
             Wx1=(2-x1)*(2-x1)*(2-x1)/6.0;
             Wx2=(4-6*x2*x2+3*x2*x2*x2)/6.0;
             Wx3=(4-6*x3*x3+3*x3*x3*x3)/6.0;
             Wx4=(2-x4)*(2-x4)*(2-x4)/6.0;
             Wy1=(2-y1)*(2-y1)*(2-y1)/6.0;
             Wy2=(4-6*y2*y2+3*y2*y2*y2)/6.0;
             Wy3=(4-6*y3*y3+3*y3*y3*y3)/6.0;
             Wy4=(2-y4)*(2-y4)*(2-y4)/6.0;
/*
             Wx1=a*(x1*x1*x1-5*x1*x1+8*x1-4);
             Wx2=(a+2)*x2*x2*x2-(a+3)*x2*x2+1;
             Wx3=(a+2)*x3*x3*x3-(a+3)*x3*x3+1;
             Wx4=a*(x4*x4*x4-5*x4*x4+8*x4-4);
             Wy1=a*(y1*y1*y1-5*y1*y1+8*y1-4);
             Wy2=(a+2)*y2*y2*y2-(a+3)*y2*y2+1;
             Wy3=(a+2)*y3*y3*y3-(a+3)*y3*y3+1;
             Wy4=a*(y4*y4*y4-5*y4*y4+8*y4-4);
*/
             Pr=Wx1*Wy1*field[i1-2][j1-2].Pr+Wx2*Wy1*field[i1-1][j1-2].Pr+Wx3*Wy1*field[i1+0][j1-2].Pr+Wx4*Wy1*field[i1+1][j1-2].Pr
               +Wx1*Wy2*field[i1-2][j1-1].Pr+Wx2*Wy2*field[i1-1][j1-1].Pr+Wx3*Wy2*field[i1+0][j1-1].Pr+Wx4*Wy2*field[i1+1][j1-1].Pr
               +Wx1*Wy3*field[i1-2][j1+0].Pr+Wx2*Wy3*field[i1-1][j1+0].Pr+Wx3*Wy3*field[i1+0][j1+0].Pr+Wx4*Wy3*field[i1+1][j1+0].Pr
               +Wx1*Wy4*field[i1-2][j1+1].Pr+Wx2*Wy4*field[i1-1][j1+1].Pr+Wx3*Wy4*field[i1+0][j1+1].Pr+Wx3*Wy4*field[i1+1][j1+1].Pr;

             Pl=Wx1*Wy1*field[i1-2][j1-2].Pl+Wx2*Wy1*field[i1-1][j1-2].Pl+Wx3*Wy1*field[i1+0][j1-2].Pl+Wx4*Wy1*field[i1+1][j1-2].Pl
               +Wx1*Wy2*field[i1-2][j1-1].Pl+Wx2*Wy2*field[i1-1][j1-1].Pl+Wx3*Wy2*field[i1+0][j1-1].Pl+Wx4*Wy2*field[i1+1][j1-1].Pl
               +Wx1*Wy3*field[i1-2][j1+0].Pl+Wx2*Wy3*field[i1-1][j1+0].Pl+Wx3*Wy3*field[i1+0][j1+0].Pl+Wx4*Wy3*field[i1+1][j1+0].Pl
               +Wx1*Wy4*field[i1-2][j1+1].Pl+Wx2*Wy4*field[i1-1][j1+1].Pl+Wx3*Wy4*field[i1+0][j1+1].Pl+Wx3*Wy4*field[i1+1][j1+1].Pl;

             Sr=Wx1*Wy1*field[i1-2][j1-2].Sr+Wx2*Wy1*field[i1-1][j1-2].Sr+Wx3*Wy1*field[i1+0][j1-2].Sr+Wx4*Wy1*field[i1+1][j1-2].Sr
               +Wx1*Wy2*field[i1-2][j1-1].Sr+Wx2*Wy2*field[i1-1][j1-1].Sr+Wx3*Wy2*field[i1+0][j1-1].Sr+Wx4*Wy2*field[i1+1][j1-1].Sr
               +Wx1*Wy3*field[i1-2][j1+0].Sr+Wx2*Wy3*field[i1-1][j1+0].Sr+Wx3*Wy3*field[i1+0][j1+0].Sr+Wx4*Wy3*field[i1+1][j1+0].Sr
               +Wx1*Wy4*field[i1-2][j1+1].Sr+Wx2*Wy4*field[i1-1][j1+1].Sr+Wx3*Wy4*field[i1+0][j1+1].Sr+Wx3*Wy4*field[i1+1][j1+1].Sr;

             Sl=Wx1*Wy1*field[i1-2][j1-2].Sl+Wx2*Wy1*field[i1-1][j1-2].Sl+Wx3*Wy1*field[i1+0][j1-2].Sl+Wx4*Wy1*field[i1+1][j1-2].Sl
               +Wx1*Wy2*field[i1-2][j1-1].Sl+Wx2*Wy2*field[i1-1][j1-1].Sl+Wx3*Wy2*field[i1+0][j1-1].Sl+Wx4*Wy2*field[i1+1][j1-1].Sl
               +Wx1*Wy3*field[i1-2][j1+0].Sl+Wx2*Wy3*field[i1-1][j1+0].Sl+Wx3*Wy3*field[i1+0][j1+0].Sl+Wx4*Wy3*field[i1+1][j1+0].Sl
               +Wx1*Wy4*field[i1-2][j1+1].Sl+Wx2*Wy4*field[i1-1][j1+1].Sl+Wx3*Wy4*field[i1+0][j1+1].Sl+Wx3*Wy4*field[i1+1][j1+1].Sl;

             x=p->x;           y=p->y; 
             x1=1+x;
             x2=x;
             x3=1-x;
             x4=2-x;
             y1=1+y;
             y2=y;
             y3=1-y;
             y4=2-y;
             Wx1=(2-x1)*(2-x1)*(2-x1)/6.0;
             Wx2=(4-6*x2*x2+3*x2*x2*x2)/6.0;
             Wx3=(4-6*x3*x3+3*x3*x3*x3)/6.0;
             Wx4=(2-x4)*(2-x4)*(2-x4)/6.0;
             Wy1=(2-y1)*(2-y1)*(2-y1)/6.0;
             Wy2=(4-6*y2*y2+3*y2*y2*y2)/6.0;
             Wy3=(4-6*y3*y3+3*y3*y3*y3)/6.0;
             Wy4=(2-y4)*(2-y4)*(2-y4)/6.0;
/*
             Wx1=a*(x1*x1*x1-5*x1*x1+8*x1-4);
             Wx2=(a+2)*x2*x2*x2-(a+3)*x2*x2+1;
             Wx3=(a+2)*x3*x3*x3-(a+3)*x3*x3+1;
             Wx4=a*(x4*x4*x4-5*x4*x4+8*x4-4);
             Wy1=a*(y1*y1*y1-5*y1*y1+8*y1-4);
             Wy2=(a+2)*y2*y2*y2-(a+3)*y2*y2+1;
             Wy3=(a+2)*y3*y3*y3-(a+3)*y3*y3+1;
             Wy4=a*(y4*y4*y4-5*y4*y4+8*y4-4);
*/            
             E1=Wx1*Wy1*field[i-1][j-1].E1+Wx2*Wy1*field[i+0][j-1].E1+Wx3*Wy1*field[i+1][j-1].E1+Wx4*Wy1*field[i+2][j-1].E1
               +Wx1*Wy2*field[i-1][j+0].E1+Wx2*Wy2*field[i+0][j+0].E1+Wx3*Wy2*field[i+1][j+0].E1+Wx4*Wy2*field[i+2][j+0].E1
               +Wx1*Wy3*field[i-1][j+1].E1+Wx2*Wy3*field[i+0][j+1].E1+Wx3*Wy3*field[i+1][j+1].E1+Wx4*Wy3*field[i+2][j+1].E1
               +Wx1*Wy4*field[i-1][j+2].E1+Wx2*Wy4*field[i+0][j+2].E1+Wx3*Wy4*field[i+1][j+2].E1+Wx4*Wy4*field[i+2][j+2].E1;

             B1=Wx1*Wy1*field[i-1][j-1].B1+Wx2*Wy1*field[i+0][j-1].B1+Wx3*Wy1*field[i+1][j-1].B1+Wx4*Wy1*field[i+2][j-1].B1
               +Wx1*Wy2*field[i-1][j+0].B1+Wx2*Wy2*field[i+0][j+0].B1+Wx3*Wy2*field[i+1][j+0].B1+Wx4*Wy2*field[i+2][j+0].B1
               +Wx1*Wy3*field[i-1][j+1].B1+Wx2*Wy3*field[i+0][j+1].B1+Wx3*Wy3*field[i+1][j+1].B1+Wx4*Wy3*field[i+2][j+1].B1
               +Wx1*Wy4*field[i-1][j+2].B1+Wx2*Wy4*field[i+0][j+2].B1+Wx3*Wy4*field[i+1][j+2].B1+Wx4*Wy4*field[i+2][j+2].B1;

             p->E1=E1+extE1; p->Pr=Pr+extPr; p->Pl=Pl+extPl;
             p->B1=B1+extB1; p->Sr=Sr+extSr; p->Sl=Sl+extSl;

             p=p->next;
             cnt++;
           }
         }		//for(s)        
       }		   //for(i,j)
   }           //End of fieldType=1
}

void interpolation2D_2nd(Domain *D,External *Ext)  //bicubic
{
   int i,j,ii,jj,istart,iend,jstart,jend,s,cnt;
   float x,y,xx,yy,extE1,extB1,extPr,extPl,extSr,extSl,alpha,beta;
   float Wx1,Wx2,Wx3,Wy1,Wy2,Wy3;
   ptclList *p;
   int myrank, nprocs;    
   double bicubic();

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   Particle **particle;
   particle=D->particle;

   extE1=Ext->E1;
   extB1=Ext->B1;
   extPr=Ext->Pr;
   extPl=Ext->Pl;
   extSr=Ext->Sr;
   extSl=Ext->Sl;

   int wt[16][16] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                    -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
                     2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
                     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                     0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
                     0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
                    -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
                     9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
                    -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
                     2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
                    -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
                     4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1};
 
   if(D->fieldType==1)
   {   
     FieldDSX **field;
     field=D->fieldDSX;
     double Pr[16],Pl[16],Sr[16],Sl[16],E1[16],B1[16];

     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++)
       {
         for(s=0; s<D->nSpecies; s++)
         {
           cnt=0;
           p=particle[i][j].head[s]->pt;
           while(p)
           {
             x=p->x;  y=p->y;
             xx=x+0.5-((int)(x+0.5));
             yy=y+0.5-((int)(y+0.5));
             ii=((int)(i+x+0.5))-1;
             jj=((int)(j+y+0.5))-1;

             Pr[0]=field[ii][jj].Pr;
             Pr[1]=field[ii][jj+1].Pr;
             Pr[2]=field[ii+1][jj+1].Pr;
             Pr[3]=field[ii+1][jj].Pr;
             Pr[4]=(field[ii+1][jj].Pr-field[ii-1][jj].Pr)*0.5;
             Pr[5]=(field[ii+1][jj+1].Pr-field[ii-1][jj+1].Pr)*0.5;
             Pr[6]=(field[ii+2][jj+1].Pr-field[ii][jj+1].Pr)*0.5;
             Pr[7]=(field[ii+2][jj].Pr-field[ii][jj].Pr)*0.5;
             Pr[8]=(field[ii][jj+1].Pr-field[ii][jj-1].Pr)*0.5;
             Pr[9]=(field[ii][jj+2].Pr-field[ii][jj].Pr)*0.5;
             Pr[10]=(field[ii+1][jj+2].Pr-field[ii+1][jj].Pr)*0.5;
             Pr[11]=(field[ii+1][jj+1].Pr-field[ii+1][jj-1].Pr)*0.5;
             Pr[12]=(field[ii+1][jj+1].Pr+field[ii-1][jj-1].Pr-field[ii-1][jj+1].Pr-field[ii+1][jj-1].Pr)*0.25;
             Pr[13]=(field[ii+1][jj+2].Pr+field[ii-1][jj].Pr-field[ii-1][jj+2].Pr-field[ii+1][jj].Pr)*0.25;
             Pr[14]=(field[ii+2][jj+2].Pr+field[ii][jj].Pr-field[ii][jj+2].Pr-field[ii+2][jj].Pr)*0.25;
             Pr[15]=(field[ii+2][jj+1].Pr+field[ii][jj-1].Pr-field[ii][jj+1].Pr-field[ii+2][jj-1].Pr)*0.25;

             Pl[0]=field[ii][jj].Pl;
             Pl[1]=field[ii][jj+1].Pl;
             Pl[2]=field[ii+1][jj+1].Pl;
             Pl[3]=field[ii+1][jj].Pl;
             Pl[4]=field[ii+1][jj].Pl-field[ii-1][jj].Pl;
             Pl[5]=field[ii+1][jj+1].Pl-field[ii-1][jj+1].Pl;
             Pl[6]=field[ii+2][jj+1].Pl-field[ii][jj+1].Pl;
             Pl[7]=field[ii+2][jj].Pl-field[ii][jj].Pl;
             Pl[8]=field[ii][jj+1].Pl-field[ii][jj-1].Pl;
             Pl[9]=field[ii][jj+2].Pl-field[ii][jj].Pl;
             Pl[10]=field[ii+1][jj+2].Pl-field[ii+1][jj].Pl;
             Pl[11]=field[ii+1][jj+1].Pl-field[ii+1][jj-1].Pl;
             Pl[12]=field[ii+1][jj+1].Pl+field[ii-1][jj-1].Pl-field[ii-1][jj+1].Pl-field[ii+1][jj-1].Pl;
             Pl[13]=field[ii+1][jj+2].Pl+field[ii-1][jj].Pl-field[ii-1][jj+2].Pl-field[ii+1][jj].Pl;
             Pl[14]=field[ii+2][jj+2].Pl+field[ii][jj].Pl-field[ii][jj+2].Pl-field[ii+2][jj].Pl;
             Pl[15]=field[ii+2][jj+1].Pl+field[ii][jj-1].Pl-field[ii][jj+1].Pl-field[ii+2][jj-1].Pl;

             Sr[0]=field[ii][jj].Sr;
             Sr[1]=field[ii][jj+1].Sr;
             Sr[2]=field[ii+1][jj+1].Sr;
             Sr[3]=field[ii+1][jj].Sr;
             Sr[4]=field[ii+1][jj].Sr-field[ii-1][jj].Sr;
             Sr[5]=field[ii+1][jj+1].Sr-field[ii-1][jj+1].Sr;
             Sr[6]=field[ii+2][jj+1].Sr-field[ii][jj+1].Sr;
             Sr[7]=field[ii+2][jj].Sr-field[ii][jj].Sr;
             Sr[8]=field[ii][jj+1].Sr-field[ii][jj-1].Sr;
             Sr[9]=field[ii][jj+2].Sr-field[ii][jj].Sr;
             Sr[10]=field[ii+1][jj+2].Sr-field[ii+1][jj].Sr;
             Sr[11]=field[ii+1][jj+1].Sr-field[ii+1][jj-1].Sr;
             Sr[12]=field[ii+1][jj+1].Sr+field[ii-1][jj-1].Sr-field[ii-1][jj+1].Sr-field[ii+1][jj-1].Sr;
             Sr[13]=field[ii+1][jj+2].Sr+field[ii-1][jj].Sr-field[ii-1][jj+2].Sr-field[ii+1][jj].Sr;
             Sr[14]=field[ii+2][jj+2].Sr+field[ii][jj].Sr-field[ii][jj+2].Sr-field[ii+2][jj].Sr;
             Sr[15]=field[ii+2][jj+1].Sr+field[ii][jj-1].Sr-field[ii][jj+1].Sr-field[ii+2][jj-1].Sr;

             Sl[0]=field[ii][jj].Sl;
             Sl[1]=field[ii][jj+1].Sl;
             Sl[2]=field[ii+1][jj+1].Sl;
             Sl[3]=field[ii+1][jj].Sl;
             Sl[4]=field[ii+1][jj].Sl-field[ii-1][jj].Sl;
             Sl[5]=field[ii+1][jj+1].Sl-field[ii-1][jj+1].Sl;
             Sl[6]=field[ii+2][jj+1].Sl-field[ii][jj+1].Sl;
             Sl[7]=field[ii+2][jj].Sl-field[ii][jj].Sl;
             Sl[8]=field[ii][jj+1].Sl-field[ii][jj-1].Sl;
             Sl[9]=field[ii][jj+2].Sl-field[ii][jj].Sl;
             Sl[10]=field[ii+1][jj+2].Sl-field[ii+1][jj].Sl;
             Sl[11]=field[ii+1][jj+1].Sl-field[ii+1][jj-1].Sl;
             Sl[12]=field[ii+1][jj+1].Sl+field[ii-1][jj-1].Sl-field[ii-1][jj+1].Sl-field[ii+1][jj-1].Sl;
             Sl[13]=field[ii+1][jj+2].Sl+field[ii-1][jj].Sl-field[ii-1][jj+2].Sl-field[ii+1][jj].Sl;
             Sl[14]=field[ii+2][jj+2].Sl+field[ii][jj].Sl-field[ii][jj+2].Sl-field[ii+2][jj].Sl;
             Sl[15]=field[ii+2][jj+1].Sl+field[ii][jj-1].Sl-field[ii][jj+1].Sl-field[ii+2][jj-1].Sl;

             E1[0]=field[i][j].E1;
             E1[1]=field[i][j+1].E1;
             E1[2]=field[i+1][j+1].E1;
             E1[3]=field[i+1][j].E1;
             E1[4]=field[i+1][j].E1-field[i-1][j].E1;
             E1[5]=field[i+1][j+1].E1-field[i-1][j+1].E1;
             E1[6]=field[i+2][j+1].E1-field[i][j+1].E1;
             E1[7]=field[i+2][j].E1-field[i][j].E1;
             E1[8]=field[i][j+1].E1-field[i][j-1].E1;
             E1[9]=field[i][j+2].E1-field[i][j].E1;
             E1[10]=field[i+1][j+2].E1-field[i+1][j].E1;
             E1[11]=field[i+1][j+1].E1-field[i+1][j-1].E1;
             E1[12]=field[i+1][j+1].E1+field[i-1][j-1].E1-field[i-1][j+1].E1-field[i+1][j-1].E1;
             E1[13]=field[i+1][j+2].E1+field[i-1][j].E1-field[i-1][j+2].E1-field[i+1][j].E1;
             E1[14]=field[i+2][j+2].E1+field[i][j].E1-field[i][j+2].E1-field[i+2][j].E1;
             E1[15]=field[i+2][j+1].E1+field[i][j-1].E1-field[i][j+1].E1-field[i+2][j-1].E1;

             B1[0]=field[i][j].B1;
             B1[1]=field[i][j+1].B1;
             B1[2]=field[i+1][j+1].B1;
             B1[3]=field[i+1][j].B1;
             B1[4]=field[i+1][j].B1-field[i-1][j].B1;
             B1[5]=field[i+1][j+1].B1-field[i-1][j+1].B1;
             B1[6]=field[i+2][j+1].B1-field[i][j+1].B1;
             B1[7]=field[i+2][j].B1-field[i][j].B1;
             B1[8]=field[i][j+1].B1-field[i][j-1].B1;
             B1[9]=field[i][j+2].B1-field[i][j].B1;
             B1[10]=field[i+1][j+2].B1-field[i+1][j].B1;
             B1[11]=field[i+1][j+1].B1-field[i+1][j-1].B1;
             B1[12]=field[i+1][j+1].B1+field[i-1][j-1].B1-field[i-1][j+1].B1-field[i+1][j-1].B1;
             B1[13]=field[i+1][j+2].B1+field[i-1][j].B1-field[i-1][j+2].B1-field[i+1][j].B1;
             B1[14]=field[i+2][j+2].B1+field[i][j].B1-field[i][j+2].B1-field[i+2][j].B1;
             B1[15]=field[i+2][j+1].B1+field[i][j-1].B1-field[i][j+1].B1-field[i+2][j-1].B1;

             p->Pr=bicubic(Pr,wt,xx,yy)+extPr;
             p->Pl=bicubic(Pl,wt,xx,yy)+extPl;
             p->Sr=bicubic(Sr,wt,xx,yy)+extSr;
             p->Sl=bicubic(Sl,wt,xx,yy)+extSl;
             p->E1=bicubic(E1,wt,x,y)+extE1;
             p->B1=bicubic(B1,wt,x,y)+extB1;

//             p->E1=E1+extE1; p->Pr=Pr+extPr; p->Pl=Pl+extPl;

             ii=(int)(x+i+0.5);  jj=(int)(y+j+0.5);
             alpha=x-((int)(x+0.5));           beta=y-((int)(y+0.5));             
             Wx1=0.5*(0.5-alpha)*(0.5-alpha);
             Wx2=0.75-alpha*alpha;
             Wx3=0.5*(0.5+alpha)*(0.5+alpha);
             Wy1=0.5*(0.5-beta)*(0.5-beta);
             Wy2=0.75-beta*beta;
             Wy3=0.5*(0.5+beta)*(0.5+beta);
/*
               
             alpha=x+i;  beta=y+j; 
             ii=(int)(alpha+0.5);  jj=(int)(beta+0.5);
             x=alpha-ii;           y=beta-jj;             
             Wx1=0.5*(0.5-x)*(0.5-x);
             Wx2=0.75-x*x;
             Wx3=0.5*(0.5+x)*(0.5+x);
             Wy1=0.5*(0.5-y)*(0.5-y);
             Wy2=0.75-y*y;
             Wy3=0.5*(0.5+y)*(0.5+y);
>>>>>>> dev
            
             E1=Wx1*Wy1*field[ii-1][jj-1].E1+Wx2*Wy1*field[ii+0][jj-1].E1+Wx3*Wy1*field[ii+1][jj-1].E1
               +Wx1*Wy2*field[ii-1][jj+0].E1+Wx2*Wy2*field[ii+0][jj+0].E1+Wx3*Wy2*field[ii+1][jj+0].E1
               +Wx1*Wy3*field[ii-1][jj+1].E1+Wx2*Wy3*field[ii+0][jj+1].E1+Wx3*Wy3*field[ii+1][jj+1].E1;

             B1=Wx1*Wy1*field[ii-1][jj-1].B1+Wx2*Wy1*field[ii+0][jj-1].B1+Wx3*Wy1*field[ii+1][jj-1].B1
               +Wx1*Wy2*field[ii-1][jj+0].B1+Wx2*Wy2*field[ii+0][jj+0].B1+Wx3*Wy2*field[ii+1][jj+0].B1
               +Wx1*Wy3*field[ii-1][jj+1].B1+Wx2*Wy3*field[ii+0][jj+1].B1+Wx3*Wy3*field[ii+1][jj+1].B1;

*/

             p=p->next;
             cnt++;
           }
         }		//for(s)        
       }		   //for(i,j)
   }           //End of fieldType=1

}

void interpolation2D_1st(Domain *D,External *Ext)
{
   int i,j,i1,j1,istart,iend,jstart,jend,s,cnt;
   float E1,Pr,Pl,B1,Sr,Sl,extE1,extB1,extPr,extPl,extSr,extSl,x,y,x1,y1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   Particle **particle;
   particle=D->particle;

   extE1=Ext->E1;
   extB1=Ext->B1;
   extPr=Ext->Pr;
   extPl=Ext->Pl;
   extSr=Ext->Sr;
   extSl=Ext->Sl;
   if(D->fieldType==1)
   {   
     FieldDSX **field;
     field=D->fieldDSX;

     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++)
       {
         for(s=0; s<D->nSpecies; s++)
         {
           cnt=0;
           p=particle[i][j].head[s]->pt;
           while(p)
           {
             x=p->x;  y=p->y;
             i1=(int)(i+x+0.5)-1;
             j1=(int)(j+y+0.5)-1;
             x1=x+0.5-((int)(x+0.5));
             y1=y+0.5-((int)(y+0.5));

             Pr=(1-x1)*(1-y1)*field[i1][j1].Pr+(x1)*(y1)*field[i1+1][j1+1].Pr+(1-x1)*(y1)*field[i1][j1+1].Pr+(x1)*(1-y1)*field[i1+1][j1].Pr;
             Pl=(1-x1)*(1-y1)*field[i1][j1].Pl+(x1)*(y1)*field[i1+1][j1+1].Pl+(1-x1)*(y1)*field[i1][j1+1].Pl+(x1)*(1-y1)*field[i1+1][j1].Pl;
             Sr=(1-x1)*(1-y1)*field[i1][j1].Sr+(x1)*(y1)*field[i1+1][j1+1].Sr+(1-x1)*(y1)*field[i1][j1+1].Sr+(x1)*(1-y1)*field[i1+1][j1].Sr;
             Sl=(1-x1)*(1-y1)*field[i1][j1].Sl+(x1)*(y1)*field[i1+1][j1+1].Sl+(1-x1)*(y1)*field[i1][j1+1].Sl+(x1)*(1-y1)*field[i1+1][j1].Sl;
             i1=(int)(i+x+0.5);
             j1=(int)(j+y+0.5);
             x1=x+0.5-((int)(x+0.5));
             y1=y+0.5-((int)(y+0.5));
//             x1=x+i-i1;  y1=y+j-j1;

             E1=(1-x)*(1-y)*field[i][j].E1+x*y*field[i+1][j+1].E1+x*(1-y)*field[i+1][j].E1+(1-x)*y*field[i][j+1].E1;
             B1=(1-x)*(1-y)*field[i][j].B1+x*y*field[i+1][j+1].B1+x*(1-y)*field[i+1][j].B1+(1-x)*y*field[i][j+1].B1;

             p->E1=E1+extE1; p->Pr=Pr+extPr; p->Pl=Pl+extPl;
             p->B1=B1+extB1; p->Sr=Sr+extSr; p->Sl=Sl+extSl;

             p=p->next;
             cnt++;
           }
         }		//for(s)        
       }		   //for(i,j)
   }           //End of fieldType=1

}


double bicubic(double data[], int wt[][16], double x, double y)
{
   int m,n,cnt;
   float coef[4][4] = {0,0,0,0,
                       0,0,0,0,
                       0,0,0,0,
                       0,0,0,0};
   float cl[16];
   float xx;
   double result;

   for(m=0; m<16; m++)
   {
     xx=0.0;
     for(n=0; n<16; n++) xx += wt[m][n]*data[n];
     cl[m]=xx;
   }

   cnt=0;
   for(m=0; m<4; m++)
     for(n=0; n<4; n++)
       coef[m][n]=cl[cnt++];

   result=0.0;
   for(m=3; m>=0; m--)
     result=x*result + ((coef[m][3]*y+coef[m][2])*y+coef[m][1])*y + coef[m][0];
   return result;

}

