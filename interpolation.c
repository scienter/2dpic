#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void interpolation2D_2nd(Domain *D,External *Ext)  //bicubic
{
   int i,j,ii,jj,i1,j1,istart,iend,jstart,jend,s;
   float x,y,Pr,Pl,Sr,Sl,E1,B1;
   float totalPr,totalPl,totalSr,totalSl,totalE1,totalB1;
   float Wx[3],Wy[3];
   float extPr,extPl,extSr,extSl,extE1,extB1;
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
           p=particle[i][j].head[s]->pt;
           while(p)
           {
             x=p->x;  y=p->y;
             Wx[0]=0.5*(1-x)*(1-x);
             Wx[1]=0.75-(0.5-x)*(0.5-x);
             Wx[2]=0.5*x*x;
             Wy[0]=0.5*(1-y)*(1-y);
             Wy[1]=0.75-(0.5-y)*(0.5-y);
             Wy[2]=0.5*y*y;

             totalPr=totalPl=totalSr=totalSl=0;
             for(jj=0; jj<3; jj++)
             {
               Pr=Pl=Sr=Sl=0;
               for(ii=0; ii<3; ii++)
               {
                 Pr+=field[i-1+ii][j-1+jj].Pr*Wx[ii];
                 Pl+=field[i-1+ii][j-1+jj].Pl*Wx[ii];
                 Sr+=field[i-1+ii][j-1+jj].Sr*Wx[ii];
                 Sl+=field[i-1+ii][j-1+jj].Sl*Wx[ii];
               }
               totalPr+=Pr*Wy[jj];
               totalPl+=Pl*Wy[jj];
               totalSr+=Sr*Wy[jj];
               totalSl+=Sl*Wy[jj];
             }

             i1=(int)(i+x+0.5);
             j1=(int)(j+y+0.5);
             x=i+x-i1;
             y=j+y-j1;
             Wx[0]=0.5*(0.5-x)*(0.5-x);
             Wx[1]=0.75-x*x;
             Wx[2]=0.5*(x+0.5)*(x+0.5);
             Wy[0]=0.5*(0.5-y)*(0.5-y);
             Wy[1]=0.75-y*y;
             Wy[2]=0.5*(y+0.5)*(y+0.5);

             totalE1=totalB1=0;
             for(jj=0; jj<3; jj++)
             {
               E1=B1=0;
               for(ii=0; ii<3; ii++)
               {
                 E1+=field[i1-1+ii][j1-1+jj].E1*Wx[ii];
                 B1+=field[i1-1+ii][j1-1+jj].B1*Wx[ii];
               }
               totalE1+=E1*Wy[jj];
               totalB1+=B1*Wy[jj];
             }


             p->E1=totalE1+extE1; p->Pr=totalPr+extPr; p->Pl=totalPl+extPl;
             p->B1=totalB1+extB1; p->Sr=totalSr+extSr; p->Sl=totalSl+extSl;
         
             p=p->next;
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
             x1=x+0.5-((int)(x+0.5));
             y1=y+0.5-((int)(y+0.5));
//             x1=x+i-i1;  y1=y+j-j1;

             E1=(1-x)*(1-y)*field[i][j].E1+x*y*field[i+1][j+1].E1+x*(1-y)*field[i+1][j].E1+(1-x)*y*field[i][j+1].E1;
             B1=(1-x)*(1-y)*field[i][j].B1+x*y*field[i+1][j+1].B1+x*(1-y)*field[i+1][j].B1+(1-x)*y*field[i][j+1].B1;
             Pr=(1.0-x1)*(1-y1)*field[i1-1][j1-1].Pr+(0.0+x1)*(0.0+y1)*field[i1][j1].Pr+(1.0-x1)*(0.0+y1)*field[i1-1][j1].Pr+(0.0+x1)*(1.0-y1)*field[i1][j1-1].Pr;
             Pl=(1.0-x1)*(1.0-y1)*field[i1-1][j1-1].Pl+(0.0+x1)*(0.0+y1)*field[i1][j1].Pl+(1.0-x1)*(0.0+y1)*field[i1-1][j1].Pl+(0.0+x1)*(1.0-y1)*field[i1][j1-1].Pl;
             Sr=(1.0-x1)*(1.0-y1)*field[i1-1][j1-1].Sr+(0.0+x1)*(0.0+y1)*field[i1][j1].Sr+(1.0-x1)*(0.0+y1)*field[i1-1][j1].Sr+(0.0+x1)*(1.0-y1)*field[i1][j1-1].Sr;
             Sl=(1.0-x1)*(1.0-y1)*field[i1-1][j1-1].Sl+(0.0+x1)*(0.0+y1)*field[i1][j1].Sl+(1.0-x1)*(0.0+y1)*field[i1-1][j1].Sl+(0.0+x1)*(1.0-y1)*field[i1][j1-1].Sl;

             p->E1=E1+extE1; p->Pr=Pr+extPr; p->Pl=Pl+extPl;
             p->B1=B1+extB1; p->Sr=Sr+extSr; p->Sl=Sl+extSl;

             p=p->next;
             cnt++;
           }
         }		//for(s)        
       }		   //for(i,j)
   }           //End of fieldType=1

}


