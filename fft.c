#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/*
 this computes an in-place complex to complex FFT
 x and y are the real and imaginary arrays of 2^m points.
 dir = 1 gives forward transform
 dir = -1 gives reverse transform
*/


main(int argc, char *argv[]) {
  float tmp, norEx, norEy, norEz, f, lowFre, range;
  int i, cnt, dataNum, binOrder;
  char str[100];
  FILE *in;
  void four1();
  
  if(argc < 2) {
    printf("[fft] [file] [order of 2]\n");
    exit(0);
  }

  binOrder = atof(argv[2]);

  dataNum=1;
  for (i=1; i<= binOrder; i++) {
    dataNum=dataNum*2;
  }

  in = fopen(argv[1],"r");
  cnt=0;
  while(fscanf(in,"%g %g %g %g %g", &tmp, &tmp, &tmp, &tmp, &tmp)!=EOF) {
//  while(fscanf(in,"%g %g", &tmp, &tmp)!=EOF) {
     cnt = cnt +1;
  }   
  fclose(in);

  if (dataNum != cnt) {
    printf("data number is not order of 2. \n");
    printf("present data number is %d!. \n",cnt);
    exit(0);
  }

  float x[dataNum*2], Pr[dataNum*2], Pl[dataNum*2], Sr[dataNum*2], Sl[dataNum*2];

  in = fopen(argv[1],"r");
//  fgets(str,100,in);
  for (i=0; i< dataNum ; i++) {
    fscanf(in,"%g %g %g %g %g", &x[i*2],&Pr[i*2],&Pl[i*2],&Sr[i*2],&Sl[i*2]);
//    fscanf(in,"%g %g",&t[i*2], &lamF[i*2]);
    x[i*2+1]=0.0;
    Pr[i*2+1]=0.0;
    Pl[i*2+1]=0.0;
    Sr[i*2+1]=0.0;
    Sl[i*2+1]=0.0;
  }   
  fclose(in);

  four1(Pr, dataNum, 1);
  four1(Pl, dataNum, 1);
  four1(Sr, dataNum, 1);
  four1(Sl, dataNum, 1);

  lowFre=0.0;
  range=x[(dataNum-1)*2]-x[0];
  for (i=0; i< dataNum; i++) {
    norEx=sqrt(Ex[i*2]*Ex[i*2]+Ex[i*2+1]*Ex[i*2+1]);
    norEy=sqrt(Ey[i*2]*Ey[i*2]+Ey[i*2+1]*Ey[i*2+1]);    
    norEz=sqrt(Ez[i*2]*Ez[i*2]+Ez[i*2+1]*Ez[i*2+1]);
    f=lowFre/range;
    lowFre=lowFre+1.0;
    printf("%g %g %g %g\n", f,norEx,norEy,norEz);
  }
}

void four1(float *data, int n, int isign)
{
  int nn, mmax,m,j,istep,i;
  float wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;

 /* This is the bit-reversal section of the routine.
    Exchange the two complex numbers. */
  nn = n << 1;
  j = 1;
  for (i=1; i<nn; i+=2) {
    if (j > i) {
      SWAP(data[j-1], data[i-1]);
      SWAP(data[j], data[i]);
    }
    m = n;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }  
    j += m; 
  }

/* Here begins the Danielson-Lanczos section of the routine.*/
  mmax = 2;
  while (nn > mmax) {
    istep = mmax << 1;
    theta = isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m=1; m<mmax; m+=2) {
      for (i=m; i<=nn; i+=istep) {
        j=i+mmax;
        tempr = wr*data[j-1]-wi*data[j];
        tempi = wr*data[j] + wi*data[j-1];
        data[j-1] = data[i-1] -tempr;
        data[j] = data[i] - tempi;
        data[i-1] += tempr;
        data[i] += tempi;
      }
      wr = (wtemp=wr)*wpr-wi*wpi+wr;
      wi = wi*wpr+wtemp*wpi+wi;
    }
    mmax = istep;
  }
} 
        
