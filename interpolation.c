#include <stdio.h>
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

void interpolation2D_2nd(Domain *D,External *Ext)
{
   int i,j,ii,jj,istart,iend,jstart,jend,s,cnt;
   float E1,Pr,Pl,B1,Sr,Sl,x,y,extE1,extB1,extPr,extPl,extSr,extSl,alpha,beta;
   float Wx1,Wx2,Wx3,Wy1,Wy2,Wy3;
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
             Wx1=0.5*(1-x)*(1-x);
             Wx2=0.75-(0.5-x)*(0.5-x);
             Wx3=0.5*x*x;
             Wy1=0.5*(1-y)*(1-y);
             Wy2=0.75-(0.5-y)*(0.5-y);
             Wy3=0.5*y*y;

             Pr=Wx1*Wy1*field[i-1][j-1].Pr+Wx2*Wy1*field[i+0][j-1].Pr+Wx3*Wy1*field[i+1][j-1].Pr
               +Wx1*Wy2*field[i-1][j+0].Pr+Wx2*Wy2*field[i+0][j+0].Pr+Wx3*Wy2*field[i+1][j+0].Pr
               +Wx1*Wy3*field[i-1][j+1].Pr+Wx2*Wy3*field[i+0][j+1].Pr+Wx3*Wy3*field[i+1][j+1].Pr;

             Pl=Wx1*Wy1*field[i-1][j-1].Pl+Wx2*Wy1*field[i+0][j-1].Pl+Wx3*Wy1*field[i+1][j-1].Pl
               +Wx1*Wy2*field[i-1][j+0].Pl+Wx2*Wy2*field[i+0][j+0].Pl+Wx3*Wy2*field[i+1][j+0].Pl
               +Wx1*Wy3*field[i-1][j+1].Pl+Wx2*Wy3*field[i+0][j+1].Pl+Wx3*Wy3*field[i+1][j+1].Pl;

             Sr=Wx1*Wy1*field[i-1][j-1].Sr+Wx2*Wy1*field[i+0][j-1].Sr+Wx3*Wy1*field[i+1][j-1].Sr
               +Wx1*Wy2*field[i-1][j+0].Sr+Wx2*Wy2*field[i+0][j+0].Sr+Wx3*Wy2*field[i+1][j+0].Sr
               +Wx1*Wy3*field[i-1][j+1].Sr+Wx2*Wy3*field[i+0][j+1].Sr+Wx3*Wy3*field[i+1][j+1].Sr;

             Sl=Wx1*Wy1*field[i-1][j-1].Sl+Wx2*Wy1*field[i+0][j-1].Sl+Wx3*Wy1*field[i+1][j-1].Sl
               +Wx1*Wy2*field[i-1][j+0].Sl+Wx2*Wy2*field[i+0][j+0].Sl+Wx3*Wy2*field[i+1][j+0].Sl
               +Wx1*Wy3*field[i-1][j+1].Sl+Wx2*Wy3*field[i+0][j+1].Sl+Wx3*Wy3*field[i+1][j+1].Sl;

             alpha=x+i;  beta=y+j; 
             ii=(int)(alpha+0.5);  jj=(int)(beta+0.5);
             x=alpha-ii;           y=beta-jj;             
             Wx1=0.5*(0.5-x)*(0.5-x);
             Wx2=0.75-x*x;
             Wx3=0.5*(0.5+x)*(0.5+x);
             Wy1=0.5*(0.5-y)*(0.5-y);
             Wy2=0.75-y*y;
             Wy3=0.5*(0.5+y)*(0.5+y);
            
             E1=Wx1*Wy1*field[ii-1][jj-1].E1+Wx2*Wy1*field[ii+0][jj-1].E1+Wx3*Wy1*field[ii+1][jj-1].E1
               +Wx1*Wy2*field[ii-1][jj+0].E1+Wx2*Wy2*field[ii+0][jj+0].E1+Wx3*Wy2*field[ii+1][jj+0].E1
               +Wx1*Wy3*field[ii-1][jj+1].E1+Wx2*Wy3*field[ii+0][jj+1].E1+Wx3*Wy3*field[ii+1][jj+1].E1;

             B1=Wx1*Wy1*field[ii-1][jj-1].B1+Wx2*Wy1*field[ii+0][jj-1].B1+Wx3*Wy1*field[ii+1][jj-1].B1
               +Wx1*Wy2*field[ii-1][jj+0].B1+Wx2*Wy2*field[ii+0][jj+0].B1+Wx3*Wy2*field[ii+1][jj+0].B1
               +Wx1*Wy3*field[ii-1][jj+1].B1+Wx2*Wy3*field[ii+0][jj+1].B1+Wx3*Wy3*field[ii+1][jj+1].B1;

             p->E1=E1+extE1; p->Pr=Pr+extPr; p->Pl=Pl+extPl;
             p->B1=B1+extB1; p->Sr=Sr+extSr; p->Sl=Sl+extSl;

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
             i1=(int)(i+x+0.5);
             j1=(int)(j+y+0.5);
             x1=x+i-i1;  y1=y+j-j1;

             E1=(1-x)*(1-y)*field[i][j].E1+x*y*field[i+1][j+1].E1+x*(1-y)*field[i+1][j].E1+(1-x)*y*field[i][j+1].E1;
             B1=(1-x)*(1-y)*field[i][j].B1+x*y*field[i+1][j+1].B1+x*(1-y)*field[i+1][j].B1+(1-x)*y*field[i][j+1].B1;
             Pr=(0.5-x1)*(0.5-y1)*field[i1-1][j1-1].Pr+(0.5+x1)*(0.5+y1)*field[i1][j1].Pr+(0.5-x1)*(0.5+y1)*field[i1-1][j1].Pr+(0.5+x1)*(0.5-y1)*field[i1][j1-1].Pr;
             Pl=(0.5-x1)*(0.5-y1)*field[i1-1][j1-1].Pl+(0.5+x1)*(0.5+y1)*field[i1][j1].Pl+(0.5-x1)*(0.5+y1)*field[i1-1][j1].Pl+(0.5+x1)*(0.5-y1)*field[i1][j1-1].Pl;
             Sr=(0.5-x1)*(0.5-y1)*field[i1-1][j1-1].Sr+(0.5+x1)*(0.5+y1)*field[i1][j1].Sr+(0.5-x1)*(0.5+y1)*field[i1-1][j1].Sr+(0.5+x1)*(0.5-y1)*field[i1][j1-1].Sr;
             Sl=(0.5-x1)*(0.5-y1)*field[i1-1][j1-1].Sl+(0.5+x1)*(0.5+y1)*field[i1][j1].Sl+(0.5-x1)*(0.5+y1)*field[i1-1][j1].Sl+(0.5+x1)*(0.5-y1)*field[i1][j1-1].Sl;

             p->E1=E1+extE1; p->Pr=Pr+extPr; p->Pl=Pl+extPl;
             p->B1=B1+extB1; p->Sr=Sr+extSr; p->Sl=Sl+extSl;

             p=p->next;
             cnt++;
           }
         }		//for(s)        
       }		   //for(i,j)
   }           //End of fieldType=1

}
