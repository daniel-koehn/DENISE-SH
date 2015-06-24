/* Calculate matrix product of two matrices
 *
 * mat1 = "left" matrix
 * mat2 = "right" matrix 
 * matr = result of the matrix product
 * n1 = number of rows of mat1
 * n2 = number of columns of mat1
 * m1 = number of rows of mat2
 * m2 = number of columns of mat2
 * sws = 1 (multiply two different matrices)
 * sws = 2 (multiply a matrix with its transpose) 
 *
 * Daniel Koehn
 * Kiel, the 30th of August 2013
 */

#include "fd.h"

void matp(float ** mat1, float ** mat2, float ** matr, int n1, int n2, int m1, int m2, int sws){

	/* extern variables */
	/*extern int NPROCX, NPROCY, MYID, POS[3];*/
	
	/* local variables */
	int i,j,h;
	float tmp,sum;
	
	/* loop over all elements of the resulting matrix*/
	/* matr = mat1*mat2 */
        if(sws==1){
          for (i=1;i<=n1;i++){
	      for (j=1;j<=m2;j++){
	    
	          matr[i][j] = 0.0;
	        
	          for (h=1;h<=n2;h++){
	              matr[i][j]+=mat1[i][h]*mat2[h][j];
                  }
	             
	       }
	   }
        }  

        /* matr = mat1t * mat1 */
	if(sws==2){
           for (i=1;i<=m2;i++){
	      for (j=1;j<=m2;j++){
	    
	          matr[i][j] = 0.0;
	       
                  for (h=1;h<=n1;h++){
                    matr[i][j]+=mat1[h][i]*mat1[h][j];
                  }  
               }   
	    
	   }
	} 	

}



