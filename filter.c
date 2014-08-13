#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void filterField(Domain *D)
{
    int i,j,a,istart,iend,jstart,jend;  
    float alpha[2];
    float nowE1,nowPr,nowPl,nowSr,nowSl;
    float beforeE1,beforePr,beforePl,beforeSr,beforeSl;
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

      for(j=0; j<jend+jstart; j++)
      {
        nowE1=field[istart-1][j].E1;
        nowPr=field[istart-1][j].Pr;
        nowPl=field[istart-1][j].Pl;
        nowSr=field[istart-1][j].Sr;
        nowSl=field[istart-1][j].Pl;
            
        for(a=0; a<2; a++)
        {
          for(i=istart; i<iend; i++)
          {              
            beforeE1=nowE1;
            beforePr=nowPr;
            beforePl=nowPl;
            beforeSr=nowSr;
            beforeSl=nowSl;
            nowE1=field[i][j].E1;
            nowPr=field[i][j].Pr;
            nowPl=field[i][j].Pl;
            nowSr=field[i][j].Sr;
            nowSl=field[i][j].Sl;
//          field[i].E1=alpha[a]*nowE1+(1.0-alpha[a])*(beforeE1+field[i+1].E1)*0.5;       
            field[i][j].Pr=alpha[a]*nowPr+(1.0-alpha[a])*(beforePr+field[i+1][j].Pr)*0.5;       
            field[i][j].Pl=alpha[a]*nowPl+(1.0-alpha[a])*(beforePl+field[i+1][j].Pl)*0.5;       
            field[i][j].Sr=alpha[a]*nowSr+(1.0-alpha[a])*(beforeSr+field[i+1][j].Sr)*0.5;     
            field[i][j].Sl=alpha[a]*nowSl+(1.0-alpha[a])*(beforeSl+field[i+1][j].Sl)*0.5;       
          }
        }		//end of for(i)
      }           //end of for(j)
    } 			//end of fieldType==1
    
    

}


