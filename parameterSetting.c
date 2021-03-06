#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "plasma.h"
#include "mesh.h"
#include "constants.h"
#include <math.h>


void parameterSetting(Domain *D,External *Ext, char *input)
{
   int FindParameters();
   int findLoadParameters();
   int findLaserParameters();
   float minX,maxX,minY,maxY,positionX,positionY,factor,pMinX,pMaxX,pPosition;
   float normalB,normalE,E1,E2,E3,B1,B2,B3,dyoverdx;
   char str[100],name[100];
   int rank,minT,maxT,tmpInt,fail=0,cnt;
   int i,j,n,numProbeX,numProbeY,probeType;
   float lambda,tmpFloat,probeDx,probeDy,maxProbeX,minProbeX,maxProbeY,minProbeY;
   LoadList *LL, *New;
   LaserList *L, *LNew;

   //Field Type
   if(FindParameters("Domain",1,"FieldType",input,str)) D->fieldType=atoi(str);
   else  {
      printf("in [Domain], FieldType=?  (1:DSX, 12:DSXY, 123:DSXYZ)\n");
      fail=1;
   }
   //Current Type
   if(FindParameters("Domain",1,"CurrentType",input,str)) D->currentType=atoi(str);
   else  {
      printf("in [Domain], CurrentType=?  (1:1st, 2:2nd, 3:3rd Order)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"InterpolationType",input,str)) D->interpolationType=atoi(str);
   else 
      D->interpolationType=1;
   if(D->interpolationType==1)  { D->numShareUp=1; D->numShareDn=1; }
   else if(D->interpolationType==2)  { D->numShareUp=2; D->numShareDn=3; }


   //Boost frame
   if(FindParameters("Domain",1,"boostGamma",input,str)) D->gamma=atof(str);
   else D->gamma=1;
   if(FindParameters("Domain",1,"boostSave",input,str)) D->boostSave=atoi(str);
   else D->boostSave=0;
   if(D->gamma>1)   D->boostOn=1;
   else             D->boostOn=0;
   D->beta=sqrt(1-1.0/D->gamma/D->gamma);



   //Domain parameter setting
   if(FindParameters("Domain",1,"maxStep",input,str)) D->maxStep=atoi(str);
   else  {
      printf("Domain maxStep must be defined.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"saveStep",input,str)) D->saveStep=atoi(str);
   else  {
      printf("Domain saveStep must be defined.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"saveStart",input,str)) D->saveStart=atoi(str);
   else  {
      printf("Domain saveStart must be defined.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"dumpSave",input,str)) D->dumpSave=atoi(str);
   else  {
      printf("in [Domain], dumpSave=?  [0:off, 1:on]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"dumpStart",input,name)) D->dumpStart=atoi(name);
   else  D->dumpStart=D->saveStart;
   if(FindParameters("Domain",1,"dumpInter",input,name)) D->dumpInter=atoi(name);
   else  D->dumpInter=0;
   if(FindParameters("Domain",1,"fieldSave",input,name)) D->fieldSave=atoi(name);
   else  D->fieldSave=1;
   if(FindParameters("Domain",1,"particleSave",input,name)) D->particleSave=atoi(name);
   else  D->particleSave=1;
   if(FindParameters("Domain",1,"rhoSave",input,name)) D->rhoSave=atoi(name);
   else  D->rhoSave=1;

   if(FindParameters("Domain",1,"minX",input,str)) 
   {
      D->minX=atof(str);
      D->minX*=D->gamma*(1+D->beta);
   }
   else  {
      printf("minX value must be defined.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"maxX",input,str)) 
   {
      maxX=atof(str);
      maxX*=D->gamma*(1+D->beta);
   }
   else  {
      printf("maxX value must be defined.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"minY",input,str)) 
      D->minY=atof(str);
   else  {
      printf("minY value must be defined.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"maxY",input,str)) 
      maxY=atof(str);
   else  {
      printf("maxY value must be defined.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"dy_over_dx",input,str)) 
      dyoverdx=atof(str);
   else  {
      printf("in [Domain], dy_over_dx=?  [dy/dx]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"divisionDirection",input,str)) D->divDirection=atof(str);
   else  {
      printf("in [Domain], divisionDirection=? How many cores do you use?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"moving",input,str)) D->moving=atoi(str);
   else  {
      printf("moving must be defined.\n");
      printf("1:moving domain on,  0:moving domain off\n");
      fail=1;
   }   
   if(FindParameters("Domain",1,"lambda",input,str)) 
   {
      D->lambda=atof(str);
      D->lambda*=D->gamma*(1+D->beta);
   }
   else  {
      printf("in [Domain], lambda=? [m]\n");
      printf("basic parameter for dx, dt.\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"divisionLambda",input,str)) 
   {
      D->divisionLambda=atof(str);
   }
   else  {
      printf("In [Domain], divisionLambda=? [number of devided wavelength]\n");
      fail=1;
   }

   //pml
   if(FindParameters("Domain",1,"pmlOn",input,str)) D->pmlOn=atoi(str);
   else  {
      printf("in [Domain], pmlOn=? \n");
      fail=1;
   }
   if(D->pmlOn==1)
   {
     if(FindParameters("Domain",1,"pmlCell",input,str)) D->pmlCell=atoi(str);
     else  {
        printf("in [Domain], pmlCell=? [ea]\n");
        fail=1;
     }
   }


   //External field parameter setting
   if(FindParameters("External",1,"E1",input,str)) E1=atof(str);
   else  {
      printf("in [External], E1=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"E2",input,str)) E2=atof(str);
   else  {
      printf("in [External], E2=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"E3",input,str)) E3=atof(str);
   else  {
      printf("in [External], E3=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"B1",input,str)) B1=atof(str);
   else  {
      printf("in [External], B1=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"B2",input,str)) B2=atof(str);
   else  {
      printf("in [External], B2=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"B3",input,str)) B3=atof(str);
   else  {
      printf("in [External], B3=? [Tesla]\n");
      fail=1;
   }


   //additional Domain parameters  
   D->dx=1.0/D->divisionLambda;
   D->dt=D->dx;   
   D->dy=D->dx*dyoverdx;
/*
   if(D->gamma>1) 
   {
     D->dy=D->dy*D->divisionLambda*0.5/dyoverdx/D->gamma/(1.0+D->beta);
     D->dx=D->dy*0.5;
     printf("labframe, dx=%g, dy=%g\n",D->dx,D->dy);
   }
*/

   tmpFloat=D->dx/(D->gamma*(1+D->beta));
   D->dy=tmpFloat*dyoverdx;
   printf("dx=%g, dy=%g\n",D->dx,D->dy);
   if(D->dy<=D->dx*1.5)   {
      printf("dyOverDx is too low. It must be over than %g.\n", D->gamma*(1+D->beta)*1.5);
      fail=1;
   }   

   D->nx=((int)((maxX-D->minX)/D->lambda/D->dx));
   D->ny=((int)((maxY-D->minY)/D->lambda/D->dy));
   D->omega=2*pi*velocityC/D->lambda;
   if(D->boostOn==1)   {
      D->minXSub=-D->nx;
      D->minYSub=0;
   }
   else   {
      D->minXSub=0;
      D->minYSub=(int)(D->minY/D->lambda/D->dy);
   }

   //Probe parameter
   if(FindParameters("Probe",1,"probeType",input,str)) probeType=atoi(str);
   else probeType=0;
   D->probeNum=0;

   if(probeType==0)
   {
     if(FindParameters("Probe",1,"probeNum",input,str)) D->probeNum=atoi(str);
     else  {
       printf("in [Probe], probeNum=?  [ea]\n");
       fail=1;
     }
     if(D->probeNum>0)
     {
       D->probeX=(int *)malloc(D->probeNum*sizeof(int));
       D->probeY=(int *)malloc(D->probeNum*sizeof(int));
       for(i=0; i<D->probeNum; i++)
       {
         sprintf(name,"probeX%d",i);
         if(FindParameters("Probe",1,name,input,str))   
           D->probeX[i]=((int)(atof(str)/D->lambda/D->dx));      
         else  {
           printf("in [Probe], probeX%d=?\n",i);
           fail=1;
         }
         sprintf(name,"probeY%d",i);
         if(FindParameters("Probe",1,name,input,str))      
           D->probeY[i]=((int)(atof(str)/D->lambda/D->dy));      
         else  {
           printf("in [Probe], probeY%d=?\n",i);
           fail=1;
         }
       }
     }  
   }
   else if(probeType==1)
   {
     if(FindParameters("Probe",1,"minProbeX",input,str)) minProbeX=atof(str);
     else  {
       printf("in [Probe], minProbeX=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"maxProbeX",input,str)) maxProbeX=atof(str);
     else  {
       printf("in [Probe], maxProbeX=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"numProbeX",input,str)) numProbeX=atoi(str);
     else 
       numProbeX=1;
     if(FindParameters("Probe",1,"minProbeY",input,str)) minProbeY=atof(str);
     else  {
       printf("in [Probe], minProbeY=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"maxProbeY",input,str)) maxProbeY=atof(str);
     else  {
       printf("in [Probe], maxProbeY=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"numProbeY",input,str)) numProbeY=atoi(str);
     else 
       numProbeY=1;
     if(numProbeX==0 && numProbeY==0)  {
       printf("in [Probe], it must be that numProbeX or numProbeY > 0 !!\n");
       fail=1;
     }
           
     probeDx=(maxProbeX-minProbeX)/((float)numProbeX);
     probeDy=(maxProbeY-minProbeY)/((float)numProbeY);
     D->probeNum=numProbeX*numProbeY;
     D->probeX=(int *)malloc(D->probeNum*sizeof(int));
     D->probeY=(int *)malloc(D->probeNum*sizeof(int));

     n=0;
     for(i=0; i<numProbeX; i++)
       for(j=0; j<numProbeY; j++)
       {       
         tmpFloat=minProbeX+i*probeDx;
         D->probeX[n]=((int)(tmpFloat/D->lambda/D->dx));      
         tmpFloat=minProbeY+j*probeDy;
         D->probeY[n]=((int)(tmpFloat/D->lambda/D->dy));     
         n++;
       } 
   }			//End of else if(probeType=2)
    
   //additional Boost parameters
   factor=D->gamma*(1+D->beta);
   D->minT=(int)(D->maxStep/factor/factor); 	//boost frame iteration
   D->maxT=(int)(D->gamma*(D->maxStep/factor+D->beta*D->nx/factor-factor*D->beta*D->minT));	//boost frame iteration
printf("maxT=%d\n",D->maxT);


   //additional external field parameters
   normalB=eMass*D->omega/(-eCharge);
   normalE=normalB*velocityC;
   Ext->Pr=(E2/normalE+B3/normalB)*0.5;
   Ext->Pl=(E2/normalE-B3/normalB)*0.5;
   Ext->E1=E1/normalE;
   Ext->Sr=(E3/normalE-B2/normalB)*0.5;
   Ext->Sl=(E3/normalE+B2/normalB)*0.5;
   Ext->B1=B1/normalB;

   //Laser parameter setting
   D->laserList = (LaserList *)malloc(sizeof(LaserList));
   D->laserList->next = NULL;
   L = D->laserList;
   rank = 1;
   while(findLaserParameters(rank,L,D,input)) 
   {
      LNew = (LaserList *)malloc(sizeof(LaserList));
      LNew->next = NULL;
      L->next=LNew;
      L=L->next;
      rank ++;
   }
   D->nLaser = rank-1;

   //Plasma parameter setting
   D->loadList = (LoadList *)malloc(sizeof(LoadList));
   D->loadList->next = NULL;
   LL = D->loadList;
   rank = 1;
   while(findLoadParameters(rank, LL, D,input)) 
   {
      New = (LoadList *)malloc(sizeof(LoadList));
      New->next = NULL;
      LL->next=New;
      LL=LL->next;
      rank ++;
   }
   D->nSpecies = rank-1;

   if(fail==1)
      exit(0);

}

int findLaserParameters(int rank, LaserList *L,Domain *D,char *input)
{
   int FindParameters();
   float positionX,positionY;
   char name[100], str[100];
   int fail=0,polarity;

   if(FindParameters("Laser",rank,"polarity",input,str)) polarity=atoi(str);
   else  polarity=0;

   if(polarity)
   {
     if(FindParameters("Laser",rank,"wavelength",input,str)) 
     {
        L->lambda=atof(str);
        L->lambda*=D->gamma*(1.0+D->beta);
     }
     else  L->lambda=D->lambda;
  
     if(FindParameters("Laser",rank,"a0",input,str)) 
        L->amplitude=atof(str);
     else  {
        printf("in [Laser], a0=??\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rU",input,str)) L->rU=atof(str);
     else  {
        printf("in [Laser], rU=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rD",input,str)) L->rD=atof(str);
     else  {
        printf("in [Laser], rD=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"loadPositionX",input,str)) positionX=atof(str);
     else  {
        printf("in [Laser], loadPositionX=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"loadPositionY",input,str)) positionY=atof(str);
     else  {
        printf("in [Laser], loadPositionY=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"beamWaist",input,str)) L->beamWaist=atof(str);
     else  {
        printf("in [Laser], beamWaist=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"focus",input,str)) L->focus=atof(str);
     else  {
        printf("in [Laser], focus=?  [m]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"flat",input,str)) L->flat=atof(str);
     else  L->flat=0.0;
     if(FindParameters("Laser",rank,"direction",input,str)) L->direction=atoi(str);
     else  L->direction=1;

     //additional laser parameters
     L->polarity=polarity;
     L->omega=2*pi*velocityC/L->lambda;
     L->loadPointX=((int)(positionX/D->lambda/D->dx));   
     L->loadPointY=((int)(positionY/D->lambda/D->dy));   
     L->rayleighLength=pi/L->lambda*D->gamma*(1+D->beta)*L->beamWaist*L->beamWaist/D->lambda;
     L->beamWaist=L->beamWaist/D->lambda;
     L->focus=L->focus/D->lambda;
     if(fail==1)
        exit(0);
   }
   return polarity;
}

int findLoadParameters(int rank, LoadList *LL,Domain *D,char *input)
{
   int FindParameters();
   LoadList *New;
   int whatSpecies();
   int whatPlasmaType();
   double whatMass();
   float pointPosition,wp,pDt;
   int whatCharge();
   char name[100], str[100];
   int i,species,fail=0;

   if(FindParameters("Plasma",rank,"type",input,name)) 
      LL->type = whatPlasmaType(name);
   else LL->type=0;

   if(LL->type>0)
   {
      LL->species=species;
      if(FindParameters("Plasma",rank,"density",input,str)) 
      {
         LL->density=atof(str);
         LL->density*=D->gamma;
      }
      else  {
         printf("in [Plasma], density=? [m-3]\n");
         fail=1;
      }

/*    
      //testing optimal dx(divisionLambda) size
      wp=sqrt(LL->density*eCharge*eCharge/eMass/eps0);
      pDt=2*pi/wp;   
      pDt=pDt/(2*pi/D->omega)/20.0;
      if(D->dt>pDt)
      {
         printf("dt must be less then %g!\n",pDt);
         printf("So, divisionLambda>%g!\n",1/pDt);
         fail=1;
      }
*/

      if(FindParameters("Plasma",rank,"species",input,name)) 
         species = whatSpecies(name);
      else  species = 0;
      if(FindParameters("Plasma",rank,"numberInCell",input,str)) 
         LL->numberInCell=atoi(str);
      else  {
         printf("in [Plasma], numberInCell=? \n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"with_next_species",input,str)) 
         LL->withNextSpcs=atoi(str);
      else  {
         printf("in [Plasma], with_next_species=? (0:no, 1:yes)\n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"with_prev_species",input,str)) 
         LL->withPrevSpcs=atoi(str);
      else  {
         printf("in [Plasma], with_prev_species=? (0:no, 1:yes)\n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"startIndex",input,str)) 
         LL->index=atoi(str);
      else  {
         printf("in [Plasma], startIndex=? \n");
         printf("It is the starting particle index which may be 0 as default.\n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"Lnodes",input,str)) LL->lnodes=atoi(str);
      else  {
         printf("in [Plasma], Lnodes=?\n");
         printf("Each nodes indicates the point of plasma density changing.\n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"Tnodes",input,str)) LL->tnodes=atoi(str);
      else  {
         printf("in [Plasma], Tnodes=?\n");
         printf("Each nodes indicates the point of plasma density changing.\n");
         fail=1;
      }
   
      LL->lpoint = (float *)malloc(LL->lnodes*sizeof(float));
      LL->ln = (float *)malloc(LL->lnodes*sizeof(float));   
      for(i=0; i<LL->lnodes; i++)
      {
         sprintf(name,"X%d",i);
         if(FindParameters("Plasma",rank,name,input,str)) 
            LL->lpoint[i] = atof(str)/D->gamma/D->lambda/D->dx;
         else 
         { printf("X%d should be defined.\n",i);  fail=1; }

         sprintf(name,"Ln%d",i);
         if(FindParameters("Plasma",rank,name,input,str)) 
            LL->ln[i] = atof(str);
         else 
         { printf("Ln%d should be defined.\n",i);  fail=1; } 
      }
      LL->tpoint = (float *)malloc(LL->tnodes*sizeof(float));
      LL->tn = (float *)malloc(LL->tnodes*sizeof(float));   
      for(i=0; i<LL->tnodes; i++)
      {
         sprintf(name,"Y%d",i);
         if(FindParameters("Plasma",rank,name,input,str)) {
            LL->tpoint[i] = atof(str)/D->lambda/D->dy;
         }
         else 
         { printf("Y%d should be defined.\n",i);  fail=1; }

         sprintf(name,"Tn%d",i);
         if(FindParameters("Plasma",rank,name,input,str)) 
            LL->tn[i] = atof(str);
         else 
         { printf("Tn%d should be defined.\n",i);  fail=1; } 
      }

      if(FindParameters("Plasma",rank,"cx",input,str))  {
         LL->cx=atof(str)/D->lambda;
      }
      else   LL->cx=D->nx*D->dx; 
      if(FindParameters("Plasma",rank,"cy",input,str))  {
         LL->cy=atof(str)/D->lambda;
      }
      else   LL->cy=D->ny*0.5*D->dy;  

      if(FindParameters("Plasma",rank,"temperature",input,str))  
         LL->temperature=atof(str);
      else   LL->temperature=0.0;	

      LL->mass=whatMass(species);
      LL->charge=whatCharge(species);
      LL->criticalDensity=eps0*eMass*D->omega*D->omega/eCharge/eCharge;
      LL->superP=LL->density*D->lambda*D->dx*D->lambda*D->dy/LL->numberInCell;
      
   }	//end of if(species)

   if(fail==1)
      exit(0);

   return LL->type;
}

int whatSpecies(char *str)
{
   if(strstr(str,"Electron")) 		return Electron;
   else if(strstr(str,"HPlus0"))   	return HPlus0;
   else if(strstr(str,"HPlus1"))   	return HPlus1;
   else if(strstr(str,"HePlus0"))   	return HePlus0;
   else if(strstr(str,"HePlus1"))   	return HePlus1;
   else if(strstr(str,"HePlus2"))   	return HePlus1;
   else if(strstr(str,"CPlus0"))   	return CPlus0;
   else if(strstr(str,"CPlus1"))   	return CPlus1;
   else if(strstr(str,"CPlus2"))   	return CPlus2;
   else if(strstr(str,"CPlus3"))   	return CPlus3;
   else if(strstr(str,"CPlus4"))   	return CPlus4;
   else if(strstr(str,"CPlus5"))   	return CPlus5;
   else if(strstr(str,"CPlus6"))   	return CPlus6;
   else return 0;
}

int whatPlasmaType(char *str)
{
   if(strstr(str,"4Point"))         return Point4;
   else if(strstr(str,"Circle"))   	return Circle;
   else return 0;
}

double whatMass(int species)
{
   if(species == Electron) 		return 1;
   else if(species == HPlus0)  		return 1.00794/eMassU;
   else if(species == HPlus1)  		return (1.00794-1*eMassU)/eMassU;
   else if(species == HePlus0)  	return (4.00260-0*eMassU)/eMassU;
   else if(species == HePlus1)  	return (4.00260-1*eMassU)/eMassU;
   else if(species == HePlus2)  	return (4.00260-2*eMassU)/eMassU;
   else if(species == CPlus0)           return (12.0111-0*eMassU)/eMassU;
   else if(species == CPlus1)           return (12.0111-1*eMassU)/eMassU;
   else if(species == CPlus2)           return (12.0111-2*eMassU)/eMassU;
   else if(species == CPlus3)           return (12.0111-3*eMassU)/eMassU;
   else if(species == CPlus4)           return (12.0111-4*eMassU)/eMassU;
   else if(species == CPlus5)           return (12.0111-5*eMassU)/eMassU;
   else if(species == CPlus6)           return (12.0111-6*eMassU)/eMassU;
   else {  printf("Species' mass not defined\n");  exit(0);  }
}

int whatCharge(int species)
{
   int fail;

   if(species == Electron) 		return -1;
   else if(species == HPlus0)  		return 0;
   else if(species == HPlus1)  		return 1;
   else if(species == HePlus0)  	return 0;
   else if(species == HePlus1)  	return 1;
   else if(species == HePlus2)  	return 2;
   else if(species == CPlus0)           return 0;
   else if(species == CPlus1)           return 1;
   else if(species == CPlus2)           return 2;
   else if(species == CPlus3)           return 3;
   else if(species == CPlus4)           return 4;
   else if(species == CPlus5)           return 5;
   else if(species == CPlus6)           return 6;
   else {  printf("Species' charge not defined\n");  exit(0);  }
}

