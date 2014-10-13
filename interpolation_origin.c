#include <stdio.h>
#include "mesh.h"
#include "mpi.h"

//void linearInterpolation2D(Domain *D)
void interpolation2D_1st(Domain *D)
{
   int i,j,s,cnt;
   float E1,Pr,Pl,B1,Sr,Sl,x,y;
   float Pr_ij, Pr_ipjp, Pr_ipj, Pr_ijp;
   float Pl_ij, Pl_ipjp, Pl_ipj, Pl_ijp;
   float Sr_ij, Sr_ipjp, Sr_ipj, Sr_ijp;
   float Sl_ij, Sl_ipjp, Sl_ipj, Sl_ijp;
   ptclList *p;
   FieldDSX **field;
   field=D->fieldDSX;
   Particle **particle;
   particle=D->particle;

    int myrank, nprocs;    

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   
   for(i=D->istart; i<D->iend; i++)
     for(j=D->jstart; j<D->jend; j++)
     {
        Pr_ij=0.25*(field[i][j].Pr+field[i-1][j].Pr+field[i-1][j-1].Pr+field[i][j-1].Pr);	//i,j
        Pr_ipjp=0.25*(field[i+1][j+1].Pr+field[i][j+1].Pr+field[i][j].Pr+field[i+1][j].Pr);		//i+1,j+1
        Pr_ipj=0.25*(field[i+1][j].Pr+field[i][j].Pr+field[i][j-1].Pr+field[i+1][j-1].Pr);	//i+1,j
        Pr_ijp=0.25*(field[i][j+1].Pr+field[i-1][j+1].Pr+field[i-1][j].Pr+field[i][j].Pr);		//i,j+1
        Pl_ij=0.25*(field[i][j].Pl+field[i-1][j].Pl+field[i-1][j-1].Pl+field[i][j-1].Pl);	//i,j
        Pl_ipjp=0.25*(field[i+1][j+1].Pl+field[i][j+1].Pl+field[i][j].Pl+field[i+1][j].Pl);		//i+1,j+1
        Pl_ipj=0.25*(field[i+1][j].Pl+field[i][j].Pl+field[i][j-1].Pl+field[i+1][j-1].Pl);	//i+1,j
        Pl_ijp=0.25*(field[i][j+1].Pl+field[i-1][j+1].Pl+field[i-1][j].Pl+field[i][j].Pl);		//i,j+1
        Sr_ij=0.25*(field[i][j].Sr+field[i-1][j].Sr+field[i-1][j-1].Sr+field[i][j-1].Sr);	//i,j
        Sr_ipjp=0.25*(field[i+1][j+1].Sr+field[i][j+1].Sr+field[i][j].Sr+field[i+1][j].Sr);		//i+1,j+1
        Sr_ipj=0.25*(field[i+1][j].Sr+field[i][j].Sr+field[i][j-1].Sr+field[i+1][j-1].Sr);	//i+1,j
        Sr_ijp=0.25*(field[i][j+1].Sr+field[i-1][j+1].Sr+field[i-1][j].Sr+field[i][j].Sr);		//i,j+1
        Sl_ij=0.25*(field[i][j].Sl+field[i-1][j].Sl+field[i-1][j-1].Sl+field[i][j-1].Sl);	//i,j
        Sl_ipjp=0.25*(field[i+1][j+1].Sl+field[i][j+1].Sl+field[i][j].Sl+field[i+1][j].Sl);		//i+1,j+1
        Sl_ipj=0.25*(field[i+1][j].Sl+field[i][j].Sl+field[i][j-1].Sl+field[i+1][j-1].Sl);	//i+1,j
        Sl_ijp=0.25*(field[i][j+1].Sl+field[i-1][j+1].Sl+field[i-1][j].Sl+field[i][j].Sl);		//i,j+1
       for(s=0; s<D->nSpecies; s++)
       {
         cnt=0;
         p=particle[i][j].head[s]->pt;
         while(p)
         {
            x=p->x;  y=p->y;

            E1=(1-x)*(1-y)*field[i][j].E1+x*y*field[i+1][j+1].E1+x*(1-y)*field[i+1][j].E1+(1-x)*y*field[i][j+1].E1;
            B1=(1-x)*(1-y)*field[i][j].B1+x*y*field[i+1][j+1].B1+x*(1-y)*field[i+1][j].B1+(1-x)*y*field[i][j+1].B1;
            Pr=(1-x)*(1-y)*Pr_ij+x*y*Pr_ipjp+x*(1-y)*Pr_ipj+(1-x)*y*Pr_ijp;
            Pl=(1-x)*(1-y)*Pl_ij+x*y*Pl_ipjp+x*(1-y)*Pl_ipj+(1-x)*y*Pl_ijp;
            Sr=(1-x)*(1-y)*Sr_ij+x*y*Sr_ipjp+x*(1-y)*Sr_ipj+(1-x)*y*Sr_ijp;
            Sl=(1-x)*(1-y)*Sl_ij+x*y*Sl_ipjp+x*(1-y)*Sl_ipj+(1-x)*y*Sl_ijp;

            p->E1=E1; p->Pr=Pr; p->Pl=Pl;
            p->B1=B1; p->Sr=Sr; p->Sl=Sl;

            p=p->next;
            cnt++;
         }
       }		//for(s)        

     }		//for(i,j)

}
