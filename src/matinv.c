/* Invert square matrix
 *
 * mat1 = "left" matrix
 * n1 = number of rows and columns of mat1
 *
 * Daniel Koehn
 * Kiel, the 30th of August 2013
 */

#include "fd.h"

void matinv(float ** a, int n){

    /* extern variables */
    /*extern int NPROCX, NPROCY, MYID, POS[3];*/
	
    /* local variables */
    int m,i,j,p,q;
    float tem=0,temp=0,temp1=0,temp2=0,temp4=0,temp5=0, **b;
    b = matrix(1,2,1,2);
    
    for(i=1;i<=n;i++){
       for(j=1;j<=n;j++){
          if(i==j){b[i][j]=1.0;}
          if(i!=j){b[i][j]=0.0;}
       }
    }
        

	
	for(i=1;i<n;i++)
	{
		temp=a[i][i];
		if(temp<0)
		temp=temp*(-1);
		p=i;
		for(j=i+1;j<n;j++)
		{
			if(a[j][i]<0)
				tem=a[j][i]*(-1);
			else
				tem=a[j][i];
			if(temp<0)
				temp=temp*(-1);
			if(tem>temp)
			{
				p=j;
				temp=a[j][i];
			}
		}
		//row exchange in both the matrix
		for(j=1;j<n;j++)
		{
			temp1=a[i][j];
			a[i][j]=a[p][j];
			a[p][j]=temp1;
			temp2=b[i][j];
			b[i][j]=b[p][j];
			b[p][j]=temp2;
		}
		//dividing the row by a[i][i]
		temp4=a[i][i];
		for(j=1;j<n;j++)
		{
			a[i][j]=(float)a[i][j]/temp4;
			b[i][j]=(float)b[i][j]/temp4;
		}
		//making other elements 0 in order to make the matrix a[][] an indentity matrix and obtaining a inverse b[][] matrix
		for(q=1;q<n;q++)
		{
			if(q==i)
				continue;
			temp5=a[q][i];
			for(j=1;j<n;j++)
			{
				a[q][j]=a[q][j]-(temp5*a[i][j]);
				b[q][j]=b[q][j]-(temp5*b[i][j]);
			}
		}
	}

    /* output of the inverse matrix */
    for(i=1;i<=n;i++){
       for(j=1;j<=n;j++){
          a[i][j]=b[i][j];
       }
    }

free_matrix(b,1,2,1,2);    
    	

}



