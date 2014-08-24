#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"

main(int argc, char *argv[])
{
   int i,j,rank,time;
   float *x,*y,**E1,**E2,**E3,**B1,**B2,**B3;
   float xx,yy,e1,e2,e3,b1,b2,b3;
   int nx,ny,initial,final,timeStep,cores;
   char name[100];
   FILE *in,*out;

   if(argc < 6)
   {
      printf("bFeildMerger initial final timeStep cores nx ny\n");
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
   E1=(float **)malloc(nx*sizeof(float *));
   E2=(float **)malloc(nx*sizeof(float *));
   E3=(float **)malloc(nx*sizeof(float *));
   B1=(float **)malloc(nx*sizeof(float *));
   B2=(float **)malloc(nx*sizeof(float *));
   B3=(float **)malloc(nx*sizeof(float *));
   for(i=0; i<nx; i++)
   {
     E1[i]=(float *)malloc(ny*sizeof(float ));
     E2[i]=(float *)malloc(ny*sizeof(float ));
     E3[i]=(float *)malloc(ny*sizeof(float ));
     B1[i]=(float *)malloc(ny*sizeof(float ));
     B2[i]=(float *)malloc(ny*sizeof(float ));
     B3[i]=(float *)malloc(ny*sizeof(float ));
   }

   for(time=initial; time <= final; time+=timeStep) 
   {
     for(i=0; i<nx; i++)
       for(j=0; j<ny; j++)
       {
         x[i]=0.0;
         y[j]=0.0;    
         E1[i][j]=0.0;
         E2[i][j]=0.0;
         E3[i][j]=0.0;
         B1[i][j]=0.0;
         B2[i][j]=0.0;
         B3[i][j]=0.0;
       }

     for(i=0; i<nx; i++)
     {     
       j=0;
       for(rank=0; rank<cores; rank++)
       {
         sprintf(name,"bField%d_%d_%d",time,i,rank);
         in = fopen(name,"r");
         if(in!=NULL)
         {
           while(fscanf(in,"%g %g %g %g %g %g %g %g",&xx,&yy,&e1,&e2,&e3,&b1,&b2,&b3)!=EOF)
           {
             x[i]=xx;
             y[j]=yy;    
             E1[i][j]=e1;
             E2[i][j]=e2;
             E3[i][j]=e3;
             B1[i][j]=b1;
             B2[i][j]=b2;
             B3[i][j]=b3;
             j++;
           } 
           fclose(in);
         }
//         else fclose(in);	
       }	//End of cores
     }		//End of i

     sprintf(name,"bField%d",time);
     remove(name);
     out=fopen(name,"w");
     for(i=0; i<nx; i++)
     {
       if(x[i]>0)
       {
         for(j=0; j<ny; j++)
         {
            fprintf(out,"%g %g %g %g %g %g %g %g\n",x[i],y[j],
                       E1[i][j],E2[i][j],E3[i][j],
                       B1[i][j],B2[i][j],B3[i][j]);
         }
         fprintf(out,"\n");
       }
     }
     fclose(out);
     printf("bField%d is made.\n",time);

   }	//End of for(time)

   for(i=0; i<nx; i++)
   {
     free(E1[i]);
     free(E2[i]);
     free(E3[i]);
     free(B1[i]);
     free(B2[i]);
     free(B3[i]);
   }
   free(x);
   free(y);
}

char* GetFileName(char file_path[])
{
   char *file_name;

   while(*file_path)
   {
     if(*file_path == '/' && (file_path+1) != NULL)
     {
       file_name = file_path+1;
     }
     else
     {}
     file_path++;
   }
   return file_name;
}
