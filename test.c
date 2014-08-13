#include "stdio.h"
#include "stdlib.h"
#include "math.h"

main(int argc, char *argv[]) {
  float x,absX;
  int intX,result,i;
  
  x = atof(argv[1]);
  absX=fabs(x);
  float *a;

  a = (float *)malloc(2*sizeof(float ));
  for(i=0; i<2; i++)
  {
    a[i]=i*0.1;
    printf("a[%d]=%g \n",i,a[i]);
  }


  result=((int)(x-1))-x/fabs(x);
//  result=((int)(x-1.0))+1;

  free(a);
} 
        
