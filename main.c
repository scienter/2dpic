#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "plasma.h"
#include "mpi.h"
#include <time.h>

int main(int argc, char *argv[])
{
    int i,j,k,n,iteration=0,filter,boost,filterStep,labSaveStep;
    float factor,time_spent;
    clock_t begin,end;
    double t;
    char name[100];
    FILE *out;
    Domain D;  
    LaserList *L;
    External Ext;
    int myrank, nTasks;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    void parameterSetting();
    void boundaray();    
    void saveField2D();  
    void saveRaman2D();  
    void saveParticle();
    void saveRho2D();
    void saveProbe();
    void saveDump2D();
//    void boostShot();
    void solveField1D();
    void solveField2DC_DSX();
    void solveField2D_DSX();
    void absorpbing();
    void absorpbingC();
    void periodY();
    void periodYC();
    void periodY1core();
    void periodY1coreC();
    void filterField();
    void loadLaser1D();
    void loadLaser2D();
    void boostLoadLaser2D();
    void loadPlasma2D();    
    void loadMovingPlasma2DBoost();    
    void restoreData2D();    
//    void loadPoint();    
    void interpolation2D_1st();
    void particlePush();       
    void updateCurrent2D_DSX_1st(); 
    void updateCurrent2D_DSX_2nd(); 
    void updateCurrent2D_DSX_3rd(); 
    void removeEdge2D();
    void removeEdge2DBoost();
    void removeField();
    void movingDomain2D();
    void loadMovingPlasma2D();
    void probe();
    void MPI_TransferF_DSX_Yminus();
    void MPI_TransferF_DSX_YminusC();
    void MPI_TransferF_DSX_Yplus();
    void MPI_TransferF_DSX_YplusC();
//    void MPI_TransferF_XplusFilter();
    void MPI_TransferP_Xplus();
    void MPI_TransferP_Xminus();
//    void MPI_TransferP_Moving();
    void MPI_TransferJ_DSX_Yminus();
    void MPI_TransferJ_DSX_Yplus();
    void rearrangeParticles();

//    begin=MPI_Wtime();
    begin=clock();

    if(argc < 2) 
    {  
      printf("mpirun -np N show [inputFile] [dumpNum]\n"); 
      exit(0); 
    }
    if(FindParameters("Domain",1,"filter",argv[1],name)) filter=atoi(name);
    else  filter=0;
    if(FindParameters("Domain",1,"filterStep",argv[1],name)) filterStep=atoi(name);
    else  filterStep=10;

    //parameter setting
    parameterSetting(&D,&Ext,argv[1]);

    //create mesh
    boundary(&D,&Ext);

    //load plasma or load dump file
    if(argc >= 3 && D.dumpInter==0)
    {   
      iteration = atoi(argv[2]);
      restoreData2D(&D,iteration);
      t=D.dt*iteration; 
    }
    else if(argc >= 3 && D.dumpInter==1)
    {   
      iteration = atoi(argv[2]);
      restoreDataInter2D(&D,iteration);
      t=D.dt*iteration;
      iteration++; 
      t=D.dt*iteration;
    }
    else 
    {
      loadPlasma2D(&D);
      t=0;
    }

    if(D.boostOn==1)	{
      L=D.laserList;
      while(L->next)  {
        boostLoadLaser2D(&D,L); 
        L=L->next;
      }
      MPI_TransferF_DSX_Yminus(&D,D.numShareDn);
      MPI_TransferF_DSX_Yplus(&D,D.numShareUp);
   }

    //rooping time 
   labSaveStep=D.boostSaveStep;
   factor=D.gamma*(1+D.beta);

    while(iteration<=D.maxStep)
    {
   
       //probe data
       probe(&D,iteration);
       
       
       if(D.boostOn==1 && D.boostSave==1)
       {
          if(filter==1 && iteration%filterStep==0)
            filterField(&D);       

          boostShot(&D,iteration);    
          if(iteration>=D.maxT)   
             iteration=D.maxStep+1;
       }

       //save File      
       if(iteration%D.saveStep==0 && iteration>=D.saveStart)   
       {

          if(D.fieldSave==1) { 
            saveField2D(&D,iteration);
            saveRaman2D(&D,iteration);
            if(myrank==0)
              printf("field%d is made.\n",iteration);  
          }
          if(D.particleSave==1) { 
            saveParticle(&D,iteration);
            if(myrank==0)
              printf("particle%d is made.\n",iteration);              
          }
          if(D.rhoSave==1) { 
            saveRho2D(&D,iteration);
            if(myrank==0)
              printf("rho%d is made.\n",iteration);              
          }
          if(D.dumpSave==1 && iteration>=D.dumpStart) { 
            saveDump2D(D,iteration);
            if(myrank==0)
              printf("dump%d is made.\n",iteration);              
          }

          if(D.probeNum>0) { 
            saveProbe(&D,iteration);
            if(myrank==0)
              printf("probe%d is made.\n",iteration);              
          }
       }

       //load laser
       if(D.boostOn==0)	{
         L=D.laserList;
         while(L->next)  {
           if(L->direction==1)     loadLaser2D(&D,L,t); 
           else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
           L=L->next;
         }
       }
//       else if(D.boostOn==1)	{
//         L=D.laserList;
//         while(L->next)  {
//           boostLoadLaser2D(&D,L,t); 
//           L=L->next;
//         }
//         MPI_TransferF_DSX_Yminus(&D,D.numShareDn);
//         MPI_TransferF_DSX_Yplus(&D,D.numShareUp);
//       }


//       if(nTasks==1)  periodY1core(&D);   
//       else           periodY(&D);
       if(D.fieldType==1)
       {
         if(D.pmlOn==1)   absorpbing(&D);
         solveField2DC_DSX(&D);
         MPI_TransferF_DSX_YminusC(&D,D.numShareDn);
         MPI_TransferF_DSX_YplusC(&D,D.numShareUp);
//         if()  periodY1coreC(&D);   

         if(D.pmlOn==1)   absorpbingC(&D);
         solveField2D_DSX(&D);
         MPI_TransferF_DSX_Yminus(&D,D.numShareDn);
         MPI_TransferF_DSX_Yplus(&D,D.numShareUp);
       }

       if(D.interpolationType==1)
         interpolation2D_1st(&D,&Ext);
       else if(D.interpolationType==2)
         interpolation2D_2nd(&D,&Ext);
//       else if(D.currentType==3)
//         interpolation2D_3rd(&D,&Ext);


       particlePush2D(&D);

       if(D.fieldType==1)
       {

         if(D.currentType==1)
           updateCurrent2D_DSX_1st(&D);
         else if(D.currentType==2)
           updateCurrent2D_DSX_2nd(&D);
         else if(D.currentType==3)
           updateCurrent2D_DSX_3rd(&D);
         MPI_TransferJ_DSX_Yplus(&D);
         MPI_TransferJ_DSX_Yminus(&D);

       }

       if (iteration>=D.nx && D.moving==1 && D.boostOn==0)
       {
          movingDomain2D(&D);
          loadMovingPlasma2D(&D);          
          rearrangeParticles2D(&D);
          MPI_TransferP_Yminus(&D);
          MPI_TransferP_Yplus(&D);
          removeEdge2D(&D);
       }
       else if(D.boostOn==1) {
          loadMovingPlasma2DBoost(&D);          
          movingDomain2D(&D);
          rearrangeParticles2D(&D);
          MPI_TransferP_Yminus(&D);
          MPI_TransferP_Yplus(&D);
          removeEdge2DBoost(&D);
       }
       else
       {
          rearrangeParticles2D(&D);
          MPI_TransferP_Yminus(&D);
          MPI_TransferP_Yplus(&D);
          removeEdge2D(&D);
       }

       //time update
       t+=D.dt;  
       if(iteration%10==0 && myrank==0)  
          printf("iteration = %d\n",iteration);           
       iteration+=1;

    }     //end of time roop                  

    end=clock();

    if(myrank==0)
    {
      time_spent=(end-begin)/CLOCKS_PER_SEC;

      //make 'report' file
      sprintf(name,"report");
      out = fopen(name,"w");
      fprintf(out,"nx=%d\n",D.nx);
      fprintf(out,"ny=%d\n",D.ny);
      fprintf(out,"cores=%d\n",nTasks);
      fprintf(out,"running time=%gm\n",time_spent/60.0);
      fclose(out);
    }


    clean2D(&D);
  
    MPI_Finalize();

    return 0;
}
