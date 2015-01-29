#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "plasma.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

// Temporary routine to dump and restore data
void saveDump2D(Domain D,int iteration)
{
   FILE *out;
   char name[100];
   int i,j,istart,iend,jstart,jend,n,s,cntParticle,cntField;
   Particle **particle;
   particle=D.particle;
   Probe **probe;
   probe=D.probe;
   ptclList *p;
   LoadList *LL;
   int myrank, nprocs;       

   istart=D.istart;
   iend=D.iend;
   jstart=D.jstart;
   jend=D.jend;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(name,"dump%d_%d",iteration,myrank);
   out = fopen(name, "w");   

   // Save simulation Domain information
   fwrite(&(D.minXSub),sizeof(int),1,out);
   fwrite(&(D.maxXSub),sizeof(int),1,out);
   fwrite(&(D.minYSub),sizeof(int),1,out);
   fwrite(&(D.maxYSub),sizeof(int),1,out);
   fwrite(&(D.maxStep),sizeof(int),1,out);
   fwrite(&(D.dx),sizeof(float),1,out);
   fwrite(&(D.dy),sizeof(float),1,out);
   fwrite(&(D.istart),sizeof(int),1,out);
   fwrite(&(D.iend),sizeof(int),1,out);
   fwrite(&(D.jstart),sizeof(int),1,out);
   fwrite(&(D.jend),sizeof(int),1,out);
   fwrite(&(D.lambda),sizeof(float),1,out);
   fwrite(&(D.nxSub),sizeof(int),1,out);
   fwrite(&(D.nySub),sizeof(int),1,out);
   
   LL=D.loadList;
   while(LL->next)
   {
     fwrite(&(LL->numberInCell),sizeof(float),1,out);
     LL=LL->next;
   }

   // Save informations of particles inside the domain
   cntParticle=0;
   for(s=0; s<D.nSpecies; s++)
     for (i=istart; i<iend; i++) 
       for (j=jstart; j<jend; j++)
       {
         p = particle[i][j].head[s]->pt;
         while(p)  { 
           p=p->next; 
           cntParticle++;    
         }
       }
       fwrite(&cntParticle,sizeof(int),1,out);

   for(s=0; s<D.nSpecies; s++)
     for (i=istart; i<iend; i++) 
       for (j=jstart; j<jend; j++)
       {
         p = particle[i][j].head[s]->pt;
         while(p)  
         { 
           fwrite(&(s),sizeof(int),1,out);   
           fwrite(&(i),sizeof(int),1,out);   
           fwrite(&(j),sizeof(int),1,out);   
           fwrite(&(p->x),sizeof(float),1,out);   
           fwrite(&(p->y),sizeof(float),1,out);   
           fwrite(&(p->p1),sizeof(float),1,out);   
           fwrite(&(p->p2),sizeof(float),1,out);   
           fwrite(&(p->p3),sizeof(float),1,out);   
           fwrite(&(p->index),sizeof(float),1,out);   

           p = p->next; 
         }            
       }		//End of for(i,j)

   if(D.fieldType==1)
   {
     FieldDSX **field;
     field=D.fieldDSX;

     cntField=(D.nxSub+5)*(D.nySub+5);
     for(i=0; i<D.nxSub+5; i++) 
       for(j=0; j<D.nySub+5; j++)
       {
         fwrite(&(i),sizeof(int),1,out);   
         fwrite(&(j),sizeof(int),1,out);   
         fwrite(&(field[i][j].E1),sizeof(float),1,out);   
         fwrite(&(field[i][j].Pr),sizeof(float),1,out);   
         fwrite(&(field[i][j].Pl),sizeof(float),1,out);   
         fwrite(&(field[i][j].B1),sizeof(float),1,out);   
         fwrite(&(field[i][j].Sr),sizeof(float),1,out);   
         fwrite(&(field[i][j].Sl),sizeof(float),1,out);   
         fwrite(&(field[i][j].E1C),sizeof(float),1,out);   
         fwrite(&(field[i][j].PrC),sizeof(float),1,out);   
         fwrite(&(field[i][j].PlC),sizeof(float),1,out);   
         fwrite(&(field[i][j].B1C),sizeof(float),1,out);   
         fwrite(&(field[i][j].SrC),sizeof(float),1,out);   
         fwrite(&(field[i][j].SlC),sizeof(float),1,out);   
         fwrite(&(field[i][j].J1),sizeof(float),1,out);   
         fwrite(&(field[i][j].J2),sizeof(float),1,out);   
         fwrite(&(field[i][j].J3),sizeof(float),1,out);   
         fwrite(&(field[i][j].J1Old),sizeof(float),1,out);   
         fwrite(&(field[i][j].J2Old),sizeof(float),1,out);   
         fwrite(&(field[i][j].J3Old),sizeof(float),1,out);   
       }
   }
/*  
   //save Probe data
   if(D.probeNum>0)
   {
     for(n=0; n<D.probeNum; n++)
       for(i=0; i<=D.maxStep; i++)
       { 
         fwrite(&(probe[n][i].Pr),sizeof(float),1,out);   
         fwrite(&(probe[n][i].Pl),sizeof(float),1,out);   
         fwrite(&(probe[n][i].Sr),sizeof(float),1,out);   
         fwrite(&(probe[n][i].Sl),sizeof(float),1,out);   
         fwrite(&(probe[n][i].E1),sizeof(float),1,out);   
         fwrite(&(probe[n][i].B1),sizeof(float),1,out);   
       }
   }
*/
   fclose(out);
}


// Temporary routine to dump and restore data
void restoreData2D(Domain *D, int iteration)
{
   FILE *in;
   char name[100];
   int i,j,istart,iend,jstart,jend,n,s,cntParticle,cntField;
   float tmp;
   double tmp1;
   ptclList *p;
   Particle **particle;
   particle=D->particle;
   LoadList *LL;
   Probe **probe;
   probe=D->probe;
   int myrank, nprocs;    

   istart=D->istart;
   iend=D->iend;
   jstart=D->jstart;
   jend=D->jend;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(name,"dump%d_%d",iteration,myrank);
   in = fopen(name, "r");   

   // restore simulation Domain information
   fread(&(D->minXSub),sizeof(int),1,in);
   fread(&(D->maxXSub),sizeof(int),1,in);
   fread(&(D->minYSub),sizeof(int),1,in);
   fread(&(D->maxYSub),sizeof(int),1,in);
   fread(&(D->maxStep),sizeof(int),1,in);
   fread(&(D->dx),sizeof(float),1,in);
   fread(&(D->dy),sizeof(float),1,in);
   fread(&(D->istart),sizeof(int),1,in);
   fread(&(D->iend),sizeof(int),1,in);
   fread(&(D->jstart),sizeof(int),1,in);
   fread(&(D->jend),sizeof(int),1,in);
   fread(&(D->lambda),sizeof(float),1,in);
   fread(&(D->nxSub),sizeof(int),1,in);
   fread(&(D->nySub),sizeof(int),1,in);

   LL=D->loadList;
   while(LL->next)
   {
     fread(&(LL->numberInCell),sizeof(float),1,in);
     LL=LL->next;
   }

   // restore informations of particles inside the domain
   fread(&cntParticle,sizeof(int),1,in);

   for(n=0; n<cntParticle; n++)  
   { 
     fread(&(s),sizeof(int),1,in);
     fread(&(i),sizeof(int),1,in);
     fread(&(j),sizeof(int),1,in);
     p = (ptclList *)malloc(sizeof(ptclList)); 
     p->next = particle[i][j].head[s]->pt;
     particle[i][j].head[s]->pt = p;
     fread(&(p->x),sizeof(float),1,in);   
     fread(&(p->y),sizeof(float),1,in);   
     fread(&(p->p1),sizeof(float),1,in);   
     fread(&(p->p2),sizeof(float),1,in);   
     fread(&(p->p3),sizeof(float),1,in);   
     fread(&(p->index),sizeof(float),1,in); 
   }

   if(D->fieldType==1)
   {
     FieldDSX **field;
     field=D->fieldDSX;

     cntField=(D->nxSub+5)*(D->nySub+5);
     for(n=0; n<cntField; n++)  
     { 
       fread(&(i),sizeof(int),1,in);
       fread(&(j),sizeof(int),1,in);
       fread(&(field[i][j].E1),sizeof(float),1,in);
       fread(&(field[i][j].Pr),sizeof(float),1,in);
       fread(&(field[i][j].Pl),sizeof(float),1,in);
       fread(&(field[i][j].B1),sizeof(float),1,in);
       fread(&(field[i][j].Sr),sizeof(float),1,in);
       fread(&(field[i][j].Sl),sizeof(float),1,in);
       fread(&(field[i][j].E1C),sizeof(float),1,in);
       fread(&(field[i][j].PrC),sizeof(float),1,in);
       fread(&(field[i][j].PlC),sizeof(float),1,in);
       fread(&(field[i][j].B1C),sizeof(float),1,in);
       fread(&(field[i][j].SrC),sizeof(float),1,in);
       fread(&(field[i][j].SlC),sizeof(float),1,in);
       fread(&(field[i][j].J1),sizeof(float),1,in);
       fread(&(field[i][j].J2),sizeof(float),1,in);
       fread(&(field[i][j].J3),sizeof(float),1,in);
       fread(&(field[i][j].J1Old),sizeof(float),1,in);
       fread(&(field[i][j].J2Old),sizeof(float),1,in);
       fread(&(field[i][j].J3Old),sizeof(float),1,in);
     }
   }       //End of fieldType=1
/*
   //restore Probe data
   if(D->probeNum>0)
   {
     for(n=0; n<D->probeNum; n++)
       for(i=0; i<=maxStep; i++)
       { 
         fread(&(probe[n][i].Pr),sizeof(float),1,in);   
         fread(&(probe[n][i].Pl),sizeof(float),1,in);   
         fread(&(probe[n][i].Sr),sizeof(float),1,in);   
         fread(&(probe[n][i].Sl),sizeof(float),1,in);   
         fread(&(probe[n][i].E1),sizeof(float),1,in);   
         fread(&(probe[n][i].B1),sizeof(float),1,in);   
       }
   }
*/   
   fclose(in);
}


// Temporary routine to dump and restore data
void restoreDataInter2D(Domain *D, int iteration)
{
   FILE *in;
   char name[100];
   int i,j,ii,jj,istart,iend,jstart,jend,n,s,cntParticle,cntField,cnt;
   int minXSub,maxXSub,minYSub,maxYSub,maxStep,nxSub,nySub;
   float tmp,dx,dy,lambda,xx,yy,x,y,p1,p2,p3,index,numberInCell;
   double tmp1;
   FieldDSX **mesh;
   ptclList *p;
   Particle **particle;
   particle=D->particle;
   Probe **probe;
   probe=D->probe;
   LoadList *LL;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   sprintf(name,"dump%d_%d",iteration,myrank);
   in = fopen(name, "r");   

   // restore simulation Domain information
   fread(&(minXSub),sizeof(int),1,in);
   fread(&(maxXSub),sizeof(int),1,in);
   fread(&(minYSub),sizeof(int),1,in);
   fread(&(maxYSub),sizeof(int),1,in);
   fread(&(maxStep),sizeof(int),1,in);
   fread(&(dx),sizeof(float),1,in);
   fread(&(dy),sizeof(float),1,in);
   fread(&(istart),sizeof(int),1,in);
   fread(&(iend),sizeof(int),1,in);
   fread(&(jstart),sizeof(int),1,in);
   fread(&(jend),sizeof(int),1,in);
   fread(&(lambda),sizeof(float),1,in);
   fread(&(nxSub),sizeof(int),1,in);
   fread(&(nySub),sizeof(int),1,in);

   LL=D->loadList;
   while(LL->next)
   {
     fread(&(numberInCell),sizeof(float),1,in);
//     LL->numberInCell=numberInCell*D->dx*D->lambda*D->dy*D->lambda/dx/lambda/dy/lambda;
     LL->numberInCell=numberInCell;
     LL=LL->next;
   }

   // restore informations of particles inside the domain
   fread(&cntParticle,sizeof(int),1,in);

   cnt=0;
   for(n=0; n<cntParticle; n++)  
   { 
     fread(&(s),sizeof(int),1,in);
     fread(&(i),sizeof(int),1,in);
     fread(&(j),sizeof(int),1,in);
     fread(&(xx),sizeof(float),1,in);   
     fread(&(yy),sizeof(float),1,in);   
     xx=(i-istart+minXSub+xx)*dx*lambda;   //physical position
     yy=(j-jstart+minYSub+yy)*dy*lambda;
     x=xx/D->lambda/D->dx-(int)(xx/D->lambda/D->dx);
     y=yy/D->lambda/D->dy-(int)(yy/D->lambda/D->dy);
     i=((int)(xx/D->lambda/D->dx))-D->minXSub+D->istart;
     j=((int)(yy/D->lambda/D->dy))-D->minYSub+D->jstart;
//if(i<D->istart || i>D->iend || j<D->jstart-1 || j>D->jend)
//printf("i=%d, j=%d\n",i,j);
     fread(&(p1),sizeof(float),1,in);   
     fread(&(p2),sizeof(float),1,in);   
     fread(&(p3),sizeof(float),1,in);   
     fread(&(index),sizeof(float),1,in); 
     if(i>=D->istart && i<D->iend && j>=D->jstart && j<D->jend)
     {
       p = (ptclList *)malloc(sizeof(ptclList)); 
       p->next = particle[i][j].head[s]->pt;
       particle[i][j].head[s]->pt = p;
       p->x=x;
       p->y=y;
       p->p1=p1;
       p->p2=p2;
       p->p3=p3;
       p->index=index;
     }
     cnt++;
//printf("cnt=%d, cntParticle=%d\n",cnt,cntParticle);
   }

   mesh=(FieldDSX **)malloc((nxSub+5)*sizeof(FieldDSX *));
   for(i=0; i<nxSub+5; i++)
     mesh[i]=(FieldDSX *)malloc((nySub+5)*sizeof(FieldDSX ));

   if(D->fieldType==1)
   {
     FieldDSX **field;
     field=D->fieldDSX;

     cntField=(nxSub+5)*(nySub+5);
     for(n=0; n<cntField; n++)  
     { 
       fread(&(i),sizeof(int),1,in);
       fread(&(j),sizeof(int),1,in);
       fread(&(mesh[i][j].E1),sizeof(float),1,in);
       fread(&(mesh[i][j].Pr),sizeof(float),1,in);
       fread(&(mesh[i][j].Pl),sizeof(float),1,in);
       fread(&(mesh[i][j].B1),sizeof(float),1,in);
       fread(&(mesh[i][j].Sr),sizeof(float),1,in);
       fread(&(mesh[i][j].Sl),sizeof(float),1,in);
       fread(&(mesh[i][j].E1C),sizeof(float),1,in);
       fread(&(mesh[i][j].PrC),sizeof(float),1,in);
       fread(&(mesh[i][j].PlC),sizeof(float),1,in);
       fread(&(mesh[i][j].B1C),sizeof(float),1,in);
       fread(&(mesh[i][j].SrC),sizeof(float),1,in);
       fread(&(mesh[i][j].SlC),sizeof(float),1,in);
       fread(&(mesh[i][j].J1),sizeof(float),1,in);
       fread(&(mesh[i][j].J2),sizeof(float),1,in);
       fread(&(mesh[i][j].J3),sizeof(float),1,in);
       fread(&(mesh[i][j].J1Old),sizeof(float),1,in);
       fread(&(mesh[i][j].J2Old),sizeof(float),1,in);
       fread(&(mesh[i][j].J3Old),sizeof(float),1,in);
     }
     for(i=D->istart; i<D->iend; i++)
       for(j=D->jstart; j<D->jend; j++)
       {
         xx=(i-D->istart+D->minXSub)*D->dx*D->lambda;  //physical position
         yy=(j-D->jstart+D->minYSub)*D->dy*D->lambda;  //physical position
         ii=((int)(xx/lambda/dx))-minXSub+istart;	//data's index
         jj=((int)(yy/lambda/dy))-minYSub+jstart;	//data's index
         if(ii>istart && ii<iend && jj>jstart && jj<jend)
         { 
//printf("ii=%d, jj=%d \n",ii,jj);
           x=xx/lambda/dx-((int)(xx/lambda/dx));	//x weight
           y=yy/lambda/dy-((int)(yy/lambda/dy));	//y weight
           field[i][j].E1=(1-x)*(1-y)*mesh[ii][jj].E1
                         +(1-x)*    y*mesh[ii][jj+1].E1
                         +    x*(1-y)*mesh[ii+1][jj].E1
                         +    x*    y*mesh[ii+1][jj+1].E1;
           field[i][j].Pr=(1-x)*(1-y)*mesh[ii][jj].Pr
                         +(1-x)*    y*mesh[ii][jj+1].Pr
                         +    x*(1-y)*mesh[ii+1][jj].Pr
                         +    x*    y*mesh[ii+1][jj+1].Pr;
           field[i][j].Pl=(1-x)*(1-y)*mesh[ii][jj].Pl
                         +(1-x)*    y*mesh[ii][jj+1].Pl
                         +    x*(1-y)*mesh[ii+1][jj].Pl
                         +    x*    y*mesh[ii+1][jj+1].Pl;
           field[i][j].B1=(1-x)*(1-y)*mesh[ii][jj].B1
                         +(1-x)*    y*mesh[ii][jj+1].B1
                         +    x*(1-y)*mesh[ii+1][jj].B1
                         +    x*    y*mesh[ii+1][jj+1].B1;
           field[i][j].Sr=(1-x)*(1-y)*mesh[ii][jj].Sr
                         +(1-x)*    y*mesh[ii][jj+1].Sr
                         +    x*(1-y)*mesh[ii+1][jj].Sr
                         +    x*    y*mesh[ii+1][jj+1].Sr;
           field[i][j].Sl=(1-x)*(1-y)*mesh[ii][jj].Sl
                         +(1-x)*    y*mesh[ii][jj+1].Sl
                         +    x*(1-y)*mesh[ii+1][jj].Sl
                         +    x*    y*mesh[ii+1][jj+1].Sl;
           field[i][j].E1C=(1-x)*(1-y)*mesh[ii][jj].E1C
                          +(1-x)*    y*mesh[ii][jj+1].E1C
                          +    x*(1-y)*mesh[ii+1][jj].E1C
                          +    x*    y*mesh[ii+1][jj+1].E1C;
           field[i][j].PrC=(1-x)*(1-y)*mesh[ii][jj].PrC
                          +(1-x)*    y*mesh[ii][jj+1].PrC
                          +    x*(1-y)*mesh[ii+1][jj].PrC
                          +    x*    y*mesh[ii+1][jj+1].PrC;
           field[i][j].PlC=(1-x)*(1-y)*mesh[ii][jj].PlC
                          +(1-x)*    y*mesh[ii][jj+1].PlC
                          +    x*(1-y)*mesh[ii+1][jj].PlC
                          +    x*    y*mesh[ii+1][jj+1].PlC;
           field[i][j].B1C=(1-x)*(1-y)*mesh[ii][jj].B1C
                          +(1-x)*    y*mesh[ii][jj+1].B1C
                          +    x*(1-y)*mesh[ii+1][jj].B1C
                          +    x*    y*mesh[ii+1][jj+1].B1C;
           field[i][j].SrC=(1-x)*(1-y)*mesh[ii][jj].SrC
                          +(1-x)*    y*mesh[ii][jj+1].SrC
                          +    x*(1-y)*mesh[ii+1][jj].SrC
                          +    x*    y*mesh[ii+1][jj+1].SrC;
           field[i][j].SlC=(1-x)*(1-y)*mesh[ii][jj].SlC
                          +(1-x)*    y*mesh[ii][jj+1].SlC
                          +    x*(1-y)*mesh[ii+1][jj].SlC
                          +    x*    y*mesh[ii+1][jj+1].SlC;
           field[i][j].J1=(1-x)*(1-y)*mesh[ii][jj].J1
                         +(1-x)*    y*mesh[ii][jj+1].J1
                         +    x*(1-y)*mesh[ii+1][jj].J1
                         +    x*    y*mesh[ii+1][jj+1].J1;
           field[i][j].J2=(1-x)*(1-y)*mesh[ii][jj].J2
                         +(1-x)*    y*mesh[ii][jj+1].J2
                         +    x*(1-y)*mesh[ii+1][jj].J2
                         +    x*    y*mesh[ii+1][jj+1].J2;
           field[i][j].J3=(1-x)*(1-y)*mesh[ii][jj].J3
                         +(1-x)*    y*mesh[ii][jj+1].J3
                         +    x*(1-y)*mesh[ii+1][jj].J3
                         +    x*    y*mesh[ii+1][jj+1].J3;
           field[i][j].J1Old=(1-x)*(1-y)*mesh[ii][jj].J1Old
                            +(1-x)*    y*mesh[ii][jj+1].J1Old
                            +    x*(1-y)*mesh[ii+1][jj].J1Old
                            +    x*    y*mesh[ii+1][jj+1].J1Old;
           field[i][j].J2Old=(1-x)*(1-y)*mesh[ii][jj].J2Old
                            +(1-x)*    y*mesh[ii][jj+1].J2Old
                            +    x*(1-y)*mesh[ii+1][jj].J2Old
                            +    x*    y*mesh[ii+1][jj+1].J2Old;
           field[i][j].J3Old=(1-x)*(1-y)*mesh[ii][jj].J3Old
                            +(1-x)*    y*mesh[ii][jj+1].J3Old
                            +    x*(1-y)*mesh[ii+1][jj].J3Old
                            +    x*    y*mesh[ii+1][jj+1].J3Old;
         }		//End of (if ii and jj is in the range)
       }			//End of i and j
printf("done\n");
   }       //End of fieldType=1


   //restore Probe data
//   if(D->probeNum>0)
//   {
//     for(n=0; n<D->probeNum; n++)
//       for(i=0; i<=maxStep; i++)
//       { 
//         fread(&(probe[n][i].Pr),sizeof(float),1,in);   
//         fread(&(probe[n][i].Pl),sizeof(float),1,in);   
//         fread(&(probe[n][i].Sr),sizeof(float),1,in);   
//         fread(&(probe[n][i].Sl),sizeof(float),1,in);   
//         fread(&(probe[n][i].E1),sizeof(float),1,in);   
//         fread(&(probe[n][i].B1),sizeof(float),1,in);   
//       }
//   }
   
   fclose(in);
}

