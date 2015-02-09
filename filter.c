#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void filterField(Domain *D)
{
    int i,j,a,istart,iend,jstart,jend;  
    float alpha[2];
    float nowE1,nowB1,nowPr,nowPl,nowSr,nowSl;
    float beforeE1,beforeB1,beforePr,beforePl,beforeSr,beforeSl;
    float nowE1C,nowB1C,nowPrC,nowPlC,nowSrC,nowSlC;
    float beforeE1C,beforeB1C,beforePrC,beforePlC,beforeSrC,beforeSlC;
    int myrank,nTasks;  
 
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    alpha[0]=0.5;
    alpha[1]=1.5;

    if(D->fieldType==1)
    {    
      FieldDSX **field;
      field=D->fieldDSX;

      for(j=0; j<jend+3; j++)
      {
        nowE1=field[istart-2][j].E1;
        nowB1=field[istart-2][j].B1;
        nowPr=field[istart-2][j].Pr;
        nowPl=field[istart-2][j].Pl;
        nowSr=field[istart-2][j].Sr;
        nowSl=field[istart-2][j].Pl;
        nowE1C=field[istart-2][j].E1C;
        nowB1C=field[istart-2][j].B1C;
        nowPrC=field[istart-2][j].PrC;
        nowPlC=field[istart-2][j].PlC;
        nowSrC=field[istart-2][j].SrC;
        nowSlC=field[istart-2][j].PlC;
            
        for(a=0; a<2; a++)
        {
          for(i=istart-1; i<iend+2; i++)
          {              
            beforeE1=nowE1;
            beforeB1=nowB1;
            beforePr=nowPr;
            beforePl=nowPl;
            beforeSr=nowSr;
            beforeSl=nowSl;
            beforeE1C=nowE1C;
            beforeB1C=nowB1C;
            beforePrC=nowPrC;
            beforePlC=nowPlC;
            beforeSrC=nowSrC;
            beforeSlC=nowSlC;
            nowE1=field[i][j].E1;
            nowE1=field[i][j].E1;
            nowPr=field[i][j].Pr;
            nowPl=field[i][j].Pl;
            nowSr=field[i][j].Sr;
            nowSl=field[i][j].Sl;
            nowE1C=field[i][j].E1C;
            nowB1C=field[i][j].B1C;
            nowPrC=field[i][j].PrC;
            nowPlC=field[i][j].PlC;
            nowSrC=field[i][j].SrC;
            nowSlC=field[i][j].SlC;
            field[i][j].E1=alpha[a]*nowE1+(1.0-alpha[a])*(beforeE1+field[i+1][j].E1)*0.5;       
            field[i][j].B1=alpha[a]*nowB1+(1.0-alpha[a])*(beforeB1+field[i+1][j].B1)*0.5;       
            field[i][j].Pr=alpha[a]*nowPr+(1.0-alpha[a])*(beforePr+field[i+1][j].Pr)*0.5;       
            field[i][j].Pl=alpha[a]*nowPl+(1.0-alpha[a])*(beforePl+field[i+1][j].Pl)*0.5;       
            field[i][j].Sr=alpha[a]*nowSr+(1.0-alpha[a])*(beforeSr+field[i+1][j].Sr)*0.5;     
            field[i][j].Sl=alpha[a]*nowSl+(1.0-alpha[a])*(beforeSl+field[i+1][j].Sl)*0.5;       
            field[i][j].E1C=alpha[a]*nowE1C+(1.0-alpha[a])*(beforeE1C+field[i+1][j].E1C)*0.5;       
            field[i][j].B1C=alpha[a]*nowB1C+(1.0-alpha[a])*(beforeB1C+field[i+1][j].B1C)*0.5;       
            field[i][j].PrC=alpha[a]*nowPrC+(1.0-alpha[a])*(beforePrC+field[i+1][j].PrC)*0.5;       
            field[i][j].PlC=alpha[a]*nowPlC+(1.0-alpha[a])*(beforePlC+field[i+1][j].PlC)*0.5;       
            field[i][j].SrC=alpha[a]*nowSrC+(1.0-alpha[a])*(beforeSrC+field[i+1][j].SrC)*0.5;     
            field[i][j].SlC=alpha[a]*nowSlC+(1.0-alpha[a])*(beforeSlC+field[i+1][j].SlC)*0.5;       
          }
        }		//end of for(i)
      }           //end of for(j)
    } 			//end of fieldType==1
    
    

}


