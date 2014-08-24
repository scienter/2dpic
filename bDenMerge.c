#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"

main(int argc, char *argv[])
{
   int i,j,rank,time;
   float *x,*y,**rho;
   float xx,yy,density;
   int nx,ny,initial,final,timeStep,cores;
   char name[100];
   FILE *in,*out;

   if(argc < 6)
   {
      printf("bDenMerger initial final timeStep cores nx ny\n");
      exit(0);
   }

   initial = atoi(argv[1]);
   final = atoi(argv[2]);
   timeStep = atoi(argv[3]);
   cores = atoi(argv[4]);
   nx = atoi(argv[5]);
   ny = atoi(argv[6]);

   x=(float *)malloc(nx*sizeof(float ));
   y=(float *)malloc(ny*sizeof(float ));
   rho=(float **)malloc(nx*sizeof(float *));
   for(i=0; i<nx; i++)
     rho[i]=(float *)malloc(ny*sizeof(float ));

   for(time=initial; time <= final; time+=timeStep) 
   {
     for(i=0; i<nx; i++)
       for(j=0; j<ny; j++)
       {
         x[i]=0.0;
         y[j]=0.0;    
         rho[i][j]=0.0;
       }

     for(i=0; i<nx; i++)
     {     
       j=0;
       for(rank=0; rank<cores; rank++)
       {
         sprintf(name,"bDensity%d_%d_%d",time,i,rank);
         in = fopen(name,"r");
         if(in!=NULL)
         {
           while(fscanf(in,"%g %g %g",&xx,&yy,&density)!=EOF)
           {
             x[i]=xx;
             y[j]=yy;    
             rho[i][j]=density;
             j++;
           } 
           fclose(in);
         }
//         else fclose(in);	
       }	//End of cores
     }		//End of i

     sprintf(name,"bDensity%d",time);
     remove(name);
     out=fopen(name,"w");
     for(i=0; i<nx; i++)
     {
       if(x[i]>0)
       {
         for(j=0; j<ny; j++)
         {
            fprintf(out,"%g %g %g\n",x[i],y[j],rho[i][j]);
         }
         fprintf(out,"\n");
       }
     }
     fclose(out);
     printf("bDensity%d is made.\n",time);

   }	//End of for(time)

   for(i=0; i<nx; i++)
   {
     free(rho[i]);
   }
   free(x);
   free(y);
}

