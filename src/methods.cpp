#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

double median(double *x, int n);

extern "C" void corEdgeWeights(double * X,
    int * EDGELIST,
    int * SAMEGENE,
    double * WEIGHT,
    int *NEDGES,
    int * NOBS,
    int * NCOR)
{
    int nobs = (int)(*NOBS);
    int nedges = (int)(*NEDGES);
    int ncor = (int)(*NCOR);
    int indx, i;

    // For each edge
    for (indx = 0;indx < nedges;indx = indx + 1) {
        int to_indx = EDGELIST[indx+nedges];
        int from_indx = EDGELIST[indx];

        if(to_indx == NA_INTEGER || from_indx == NA_INTEGER){
        	WEIGHT[indx] = NA_REAL;
        	continue;
        }

        WEIGHT[indx] = 0.0; // if all else fails the weight will be assigned to -1 i.e. the most inprobable edge

        // Compute the correlation
        if (SAMEGENE[indx] == 0) {
        	double corlist[ncor];

            for(int j=0; j<ncor; j++){
				double Exy = 0.0, Exx = 0.0, Ex = 0.0, Eyy = 0.0, Ey = 0.0;
				double xp = 0.0, yp = 0.0;
				double n = (double)nobs;

				for (i = 0;i < nobs;i = i + 1) {
					if(ncor >1){  //If multiple correlations, sample the columns and take the median.
						int sample = rand() % nobs;
						xp = X[from_indx*nobs + sample]; yp = X[to_indx*nobs + sample];
					}else{
					xp = X[from_indx*nobs + i]; yp = X[to_indx*nobs + i];
					}

					if (!isnan(xp) && !isnan(yp)) {
						Ex = Ex + xp; Exx = Exx + xp*xp; Ey = Ey + yp; Eyy = Eyy + yp*yp; Exy = Exy + xp*yp;
					} else n = n - 1.0; // If it is a missing value skip that observation completely and reduce the dataset size by 1.
				}
				if (n > 2) {
					if (Exy != 0.0 && Exx != 0.0 && Eyy != 0.0 && Ex != 0.0 && Ey != 0.0)
						corlist[j] = (n*Exy - Ex*Ey)/ sqrt( (n*Exx - Ex*Ex) * (n*Eyy - Ey*Ey) );
				}
            }
            WEIGHT[indx] = median(corlist, ncor);
        } else {
            // If it is the same gene set to a correlation of 1.0
            WEIGHT[indx] = 1.0;
        }
    }
}

double median(double x[], int n) {
    double temp;
    int i, j;
    if(n==1){return(x[0]);}

    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}
