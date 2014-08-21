#include "particle.h"
#include "laser.h"
#define FIRST 1
#define SECOND 2
#define THIRD 3


typedef struct _Domain 
{
   int fieldType;
   int currentType;

   int maxStep;
   int saveStep;
   int saveStart;
   int dumpStart;
   int dumpSave;
   int fieldSave;
   int particleSave;
   int rhoSave;

   int nx;             //Total domain
   int ny;             //Total domain
   int nxSub;          //Each core has sub domain        
   int nySub;          //Each core has sub domain        
   int istart;
   int iend;
   int jstart;
   int jend;
   int minXSub;        //Each core has start mesh point in total domain
   int maxXSub;
   int minYSub;        //Each core has start mesh point in total domain
   int maxYSub;
   int numberInCell;
   int moving;         //Moving domain option. 1:on
         
   float lambda;
   float omega;
   float divisionLambda;
   float dt;
   float dx;
   float dy;

   //MPI parameter
   int divDirection;         //how many cores?
   
   double *upJ;
   double *btJ;

   struct _FieldDSX **fieldDSX;    
   struct _UpPML **UpPML;    
   struct _DnPML **DnPML;    
   struct _Particle **particle;    
   struct _Boost **boost;    

   //Plasma load
   struct _LoadList *loadList;
   int nSpecies;

   //Plasma load
   struct _LaserList *laserList;
   int nLaser;

   //Boost
   int boostOn;
   float gamma;
   float beta;
   int minT;	//boost frame's step
   int maxT;	//boost frame's step
   int boostSaveStep;	//lab frame's step
   
   //Probe
   int probeNum;
   int *probeX;
   int *probeY;
   struct _Probe **probe;
   
}  Domain; 

typedef struct _Boost
{
   float x;
   float y;
   float E1;
   float B1;
   float Pr;
   float Pl;
   float Sr;
   float Sl;   
}  Boost;

typedef struct _FieldDSX 
{
   float Pr;
   float Pl;
   float E1;
   float PrC;
   float PlC;
   float E1C;
   float Sr;
   float Sl;
   float B1;
   float SrC;
   float SlC;
   float B1C;
   float J1;
   float J2;
   float J3;
   float J1Old;
   float J2Old;
   float J3Old;

}  FieldDSX;

typedef struct _UpPML 
{
   float Dx;
   float Dy;
   float Dz;
   float Ex;
   float Ey;
   float Ez;
   float Hx;
   float Hy;
   float Hz;
   float Bx;
   float By;
   float Bz;
}  UpPML;

typedef struct _DnPML 
{
   float Bzx;
   float Bzy;
}  DnPML;

typedef struct _Particle 
{
   float rho;
   // Particle List Header
   ptclHead **head;            
}  Particle;

typedef struct _External 
{
   float Pr;
   float Pl;
   float E1;
   float Sr;
   float Sl;
   float B1;
}  External;

typedef struct _Probe
{
   float E1;
   float Pr;
   float Pl;
   float B1;
   float Sr;
   float Sl;
}  Probe;
