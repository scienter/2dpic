#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
   float minX,maxX,dx,tmp,xx,u1,u2,u3,gamma,index,superP;
   char name[100];
   int i,cnt,NX,countSp,s,step,initial,final,timeStep;
   FILE *in, *out;
   
   if (argc < 8) 
   {
      printf("density initial final timeStep minX maxX NX superP countSpecies\n");
      exit(0);
   }
   initial = atoi(argv[1]);
   final = atoi(argv[2]);
   timeStep = atoi(argv[3]);
   minX = atof(argv[4]);
   maxX = atof(argv[5]);
   NX = atoi(argv[6]);
   superP = atof(argv[7]);
   countSp = atoi(argv[8]);
   dx=(maxX-minX)/((float)(NX));

   for(s=0; s<countSp; s++)
   {
     for(step=initial; step<=final; step+=timeStep)
     {
      sprintf(name,"%dParticle%d",s,step);
      in = fopen(name,"r");
      cnt=0; 
      while(fscanf(in,"%g %g %g %g %g %g",&tmp,&tmp,&tmp,&tmp,&tmp,&tmp)!=EOF) {
         cnt=cnt+1;
      }
      fclose(in);
      float x[cnt],ne[cnt];
    

      for(i=0; i<NX; i++)
      {
         x[i]=minX+i*dx;
         ne[i]=0;
      }

      sprintf(name,"%dParticle%d",s,step);
      in = fopen(name,"r");
      while(fscanf(in,"%g %g %g %g %g %g",&xx,&u1,&u2,&u3,&gamma,&index)!=EOF)
      {
         i=((int)((xx-minX)/dx));
         ne[i]+=1;
      }
      fclose(in);

      
      sprintf(name,"%dDensity%d",s,step);
      out = fopen(name,"w");
      for (i=0; i<NX; i++) 
         fprintf(out,"%g %g\n",x[i],ne[i]*superP/dx);
      fclose(out);
      printf("%dDensity%d is saved.\n",s,step);
    }
   }
}
