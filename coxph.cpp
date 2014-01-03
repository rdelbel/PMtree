//#include <Rcpp.h>
#include <vector>
#include <math.h>
using namespace std;
using namespace Rcpp;

#define LARGE 22
#define SMALL -200

double coxsafe(double x) {
    if (x< SMALL) return(SMALL);
    if (x> LARGE) return(LARGE);
    return (x);
    }

double **dmatrix(double *array, int ncol, int nrow)
    {

    int i;
    double **pointer;

    pointer = new double * [nrow];
    for (i=0; i<nrow; i++) {
  pointer[i] = array;
	array += ncol;
	}
    return(pointer);
    }
    
void chinv2(double **matrix , int n)
     {
     register double temp;
     register int i,j,k;

     /*
     ** invert the cholesky in the lower triangle
     **   take full advantage of the cholesky's diagonal of 1's
     */
     for (i=0; i<n; i++){
    if (matrix[i][i] >0) {
	      matrix[i][i] = 1/matrix[i][i];   /*this line inverts D */
	      for (j= (i+1); j<n; j++) {
		   matrix[j][i] = -matrix[j][i];
		   for (k=0; k<i; k++)     /*sweep operator */
			matrix[j][k] += matrix[j][i]*matrix[i][k];
		   }
	      }
	  }

     /*
     ** lower triangle now contains inverse of cholesky
     ** calculate F'DF (inverse of cholesky decomp process) to get inverse
     **   of original matrix
     */
     for (i=0; i<n; i++) {
	  if (matrix[i][i]==0) {  /* singular row */
		for (j=0; j<i; j++) matrix[j][i]=0;
		for (j=i; j<n; j++) matrix[i][j]=0;
		}
	  else {
	      for (j=(i+1); j<n; j++) {
		   temp = matrix[j][i]*matrix[j][j];
		   if (j!=i) matrix[i][j] = temp;
		   for (k=i; k<j; k++)
			matrix[i][k] += temp*matrix[j][k];
		   }
	      }
	  }
     }

int cholesky2(double **matrix, int n, double toler)
    {
    double temp;
    int  i,j,k;
    double eps, pivot;
    int rank;
    int nonneg;

    nonneg=1;
    eps =0;
    for (i=0; i<n; i++) {
  if (matrix[i][i] > eps)  eps = matrix[i][i];
	for (j=(i+1); j<n; j++)  matrix[j][i] = matrix[i][j];
	}
    eps *= toler;

    rank =0;
    for (i=0; i<n; i++) {
	pivot = matrix[i][i];
	if (pivot < eps) {
	    matrix[i][i] =0;
	    if (pivot < -8*eps) nonneg= -1;
	    }
	else  {
	    rank++;
	    for (j=(i+1); j<n; j++) {
		temp = matrix[j][i]/pivot;
		matrix[j][i] = temp;
		matrix[j][j] -= temp*temp*pivot;
		for (k=(j+1); k<n; k++) matrix[k][j] -= temp*matrix[k][i];
		}
	    }
	}
    return(rank * nonneg);
    }

void chsolve2(double **matrix, int n, double *y)
     {
     register int i,j;
     register double temp;

     /*
     ** solve Fb =y
     */
     for (i=0; i<n; i++) {
    temp = y[i] ;
	  for (j=0; j<i; j++)
	       temp -= y[j] * matrix[i][j] ;
	  y[i] = temp ;
	  }
     /*
     ** solve DF'z =b
     */
     for (i=(n-1); i>=0; i--) {
	  if (matrix[i][i]==0)  y[i] =0;
	  else {
	      temp = y[i]/matrix[i][i];
	      for (j= i+1; j<n; j++)
		   temp -= y[j]*matrix[j][i];
	      y[i] = temp;
	      }
	  }
     }
     
  

/// [[Rcpp::export]]
vector<double> coxph(double * time, int * status, double ** covar,
            int nused, int nvar){

    int i,j,k, person;
    
    double **cmat, **imat;  /*ragged arrays */
    double  wtave;
    double *a, *newbeta;
    double *a2, **cmat2;
    double *scale;
    double  denom=0, zbeta, risk;
    double  temp, temp2;
    int     ndead;  /* actually, the sum of their weights */
    double  newlk=0;
    double  dtime, d2;
    double  deadwt;  /*sum of case weights for the deaths*/
    double  efronwt; /* sum of weighted risk scores for the deaths*/
    int     halving;    /*are we doing step halving at the moment? */
    int     nrisk;   /* number of subjects in the current risk set */
 

    /* vector inputs */
    double weights[nused]; //will set to 1 later
    double offset[nused]; //set to 0
    int strata[nused]; // set to 0
    
    /* returned objects */
    double beta[nvar], u[nvar], loglik[2], means[nvar];
    double *sctest;
    int *flag, *iter;

    /* hardcode scalar default paramaters */
    int method = 1; // efron
    int maxiter = 20; //default
    double eps  = 0.000000001;     /* convergence criteria */
    double toler =0.000000000001818989;  /* tolerance for cholesky */
    int doscale = 1; // default

    //initialize weights to 1
    for (int i = 0; i < nused; ++i){
        weights[i] = 1;
        offset[i]=0;
        strata[i]=0;
    }
        
    
    /*
    **  Set up the ragged arrays and scratch space
    **  Normally covar2 does not need to be duplicated, even though
    **  we are going to modify it, due to the way this routine was
    **  was called.  In this case NAMED(covar2) will =0
    */
    
    double imat2[nvar*nvar];
    imat = dmatrix(imat2,  nvar, nvar);
    //is this ok?
    // (double) 5 makes 5 a double
    a = new double[2*nvar*nvar + 4*nvar];
    newbeta = a + nvar;
    a2 = newbeta + nvar;
    scale = a2 + nvar;
    cmat = dmatrix(scale + nvar,   nvar, nvar);
    cmat2= dmatrix(scale + nvar +nvar*nvar, nvar, nvar);


    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable.
    */
    for (i=0; i<nvar; i++) {
  temp=0;
	for (person=0; person<nused; person++) temp += covar[i][person];
	temp /= nused;
	means[i] = temp;
	for (person=0; person<nused; person++) covar[i][person] -=temp;
	if (doscale==1) {  /* and also scale it */
	    temp =0;
	    for (person=0; person<nused; person++) {
		temp += fabs(covar[i][person]);
	    }
	    if (temp > 0) temp = nused/temp;   /* scaling */
	    else temp=1.0; /* rare case of a constant covariate */
	    scale[i] = temp;
	    for (person=0; person<nused; person++)  covar[i][person] *= temp;
	    }
	}
    if (doscale==1) {
	for (i=0; i<nvar; i++) beta[i] /= scale[i]; /*rescale initial betas */
	}
    else {
	for (i=0; i<nvar; i++) scale[i] = 1.0;
	}

    /*
    ** do the initial iteration step
    */
    strata[nused-1] =1;
    loglik[1] =0;
    for (i=0; i<nvar; i++) {
	u[i] =0;
	a2[i] =0;
	for (j=0; j<nvar; j++) {
	    imat[i][j] =0 ;
	    cmat2[i][j] =0;
	    }
	}

    for (person=nused-1; person>=0; ) {
	if (strata[person] == 1) {
	    nrisk =0 ;  
	    denom = 0;
	    for (i=0; i<nvar; i++) {
		a[i] = 0;
		for (j=0; j<nvar; j++) cmat[i][j] = 0;
		}
	    }

	dtime = time[person];
	ndead =0; /*number of deaths at this time point */
	deadwt =0;  /* sum of weights for the deaths */
	efronwt=0;  /* sum of weighted risks for the deaths */
	while(person >=0 &&time[person]==dtime) {
	    /* walk through the this set of tied times */
	    nrisk++;
	    zbeta = offset[person];    /* form the term beta*z (vector mult) */
	    for (i=0; i<nvar; i++)
		zbeta += beta[i]*covar[i][person];
	    zbeta = coxsafe(zbeta);
	    risk = exp(zbeta) * weights[person];
	    denom += risk;

	    /* a is the vector of weighted sums of x, cmat sums of squares */
	    for (i=0; i<nvar; i++) {
		a[i] += risk*covar[i][person];
		for (j=0; j<=i; j++)
		    cmat[i][j] += risk*covar[i][person]*covar[j][person];
	        }

	    if (status[person]==1) {
		ndead++;
		deadwt += weights[person];
		efronwt += risk;
		loglik[1] += weights[person]*zbeta;

		for (i=0; i<nvar; i++) 
		    u[i] += weights[person]*covar[i][person];
		if (method==1) { /* Efron */
		    for (i=0; i<nvar; i++) {
			a2[i] +=  risk*covar[i][person];
			for (j=0; j<=i; j++)
			    cmat2[i][j] += risk*covar[i][person]*covar[j][person];
		        }
		    }
	        }
	    
	    person--;
	    if (strata[person]==1) break;  /*ties don't cross strata */
	    }


	if (ndead >0) {  /* we need to add to the main terms */
	    if (method==0) { /* Breslow */
		loglik[1] -= deadwt* log(denom);
	   
		for (i=0; i<nvar; i++) {
		    temp2= a[i]/ denom;  /* mean */
		    u[i] -=  deadwt* temp2;
		    for (j=0; j<=i; j++)
			imat[j][i] += deadwt*(cmat[i][j] - temp2*a[j])/denom;
		    }
		}
	    else { /* Efron */
		/*
		** If there are 3 deaths we have 3 terms: in the first the
		**  three deaths are all in, in the second they are 2/3
		**  in the sums, and in the last 1/3 in the sum.  Let k go
		**  from 0 to (ndead -1), then we will sequentially use
		**     denom - (k/ndead)*efronwt as the denominator
		**     a - (k/ndead)*a2 as the "a" term
		**     cmat - (k/ndead)*cmat2 as the "cmat" term
		**  and reprise the equations just above.
		*/
		for (k=0; k<ndead; k++) {
		    temp = (double)k/ ndead;
		    wtave = deadwt/ndead;
		    d2 = denom - temp*efronwt;
		    loglik[1] -= wtave* log(d2);
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/ d2;
			u[i] -= wtave *temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] +=  (wtave/d2) *
				((cmat[i][j] - temp*cmat2[i][j]) -
					  temp2*(a[j]-temp*a2[j]));
		        }
		    }
		
		for (i=0; i<nvar; i++) {
		    a2[i]=0;
		    for (j=0; j<nvar; j++) cmat2[i][j]=0;
		    }
		}
	    }
	}   /* end  of accumulation loop */
    loglik[0] = loglik[1]; /* save the loglik for iter 0 */

    /* am I done?
    **   update the betas and test for convergence
    */
    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
	a[i] = u[i];

    *flag= cholesky2(imat, nvar, toler);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */

    temp=0;
    for (i=0; i<nvar; i++)
	temp +=  u[i]*a[i];
    *sctest = temp;  /* score test */

    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar; i++) {
	newbeta[i] = beta[i] + a[i];
	}
    if (maxiter==0) {
	chinv2(imat,nvar);
	for (i=0; i<nvar; i++) {
	    beta[i] *= scale[i];  /*return to original scale */
	    u[i] /= scale[i];
	    imat[i][i] *= scale[i]*scale[i];
	    for (j=0; j<i; j++) {
		imat[j][i] *= scale[i]*scale[j];
		imat[i][j] = imat[j][i];
		}
	    }
	goto finish;
    }

    /*
    ** here is the main loop
    */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (*iter=1; *iter<= maxiter; (*iter)++) {
	newlk =0;
	for (i=0; i<nvar; i++) {
	    u[i] =0;
	    for (j=0; j<nvar; j++)
		imat[i][j] =0;
	    }

	/*
	** The data is sorted from smallest time to largest
	** Start at the largest time, accumulating the risk set 1 by 1
	*/
	for (person=nused-1; person>=0; ) {
	    if (strata[person] == 1) { /* rezero temps for each strata */
		denom = 0;
		nrisk =0;
		for (i=0; i<nvar; i++) {
		    a[i] = 0;
		    for (j=0; j<nvar; j++) cmat[i][j] = 0;
		    }
		}

	    dtime = time[person];
	    deadwt =0;
	    ndead =0;
	    efronwt =0;
	    while(person>=0 && time[person]==dtime) {
		nrisk++;
		zbeta = offset[person];
		for (i=0; i<nvar; i++)
		    zbeta += newbeta[i]*covar[i][person];
		zbeta = coxsafe(zbeta);
		risk = exp(zbeta) * weights[person];
		denom += risk;

		for (i=0; i<nvar; i++) {
		    a[i] += risk*covar[i][person];
		    for (j=0; j<=i; j++)
		    cmat[i][j] += risk*covar[i][person]*covar[j][person];
		    }

		if (status[person]==1) {
		    ndead++;
		    deadwt += weights[person];
		    newlk += weights[person] *zbeta;
		    for (i=0; i<nvar; i++) 
			u[i] += weights[person] *covar[i][person];
		    if (method==1) { /* Efron */
			efronwt += risk;
			for (i=0; i<nvar; i++) {
			    a2[i] +=  risk*covar[i][person];
			    for (j=0; j<=i; j++)
				cmat2[i][j] += risk*covar[i][person]*covar[j][person];
			    }   
		        }
	  	    }
		
		person--;
		if (strata[person]==1) break; /*tied times don't cross strata*/
	        }

	    if (ndead >0) {  /* add up terms*/
		if (method==0) { /* Breslow */
		    newlk -= deadwt* log(denom);
		    for (i=0; i<nvar; i++) {
			temp2= a[i]/ denom;  /* mean */
			u[i] -= deadwt* temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] +=  (deadwt/denom)*
				(cmat[i][j] - temp2*a[j]);
		        }
    		    }
		else  { /* Efron */
		    for (k=0; k<ndead; k++) {
			temp = (double)k / ndead;
			wtave= deadwt/ ndead;
			d2= denom - temp* efronwt;
			newlk -= wtave* log(d2);
			for (i=0; i<nvar; i++) {
			    temp2 = (a[i] - temp*a2[i])/ d2;
			    u[i] -= wtave*temp2;
			    for (j=0; j<=i; j++)
				imat[j][i] +=  (wtave/d2)*
				    ((cmat[i][j] - temp*cmat2[i][j]) -
				    temp2*(a[j]-temp*a2[j]));
    		            }
    		        }

		    for (i=0; i<nvar; i++) { /*in anticipation */
			a2[i] =0;
			for (j=0; j<nvar; j++) cmat2[i][j] =0;
		        }
	            }
		}
	    }   /* end  of accumulation loop  */

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky2(imat, nvar, toler);

	if (fabs(1-(loglik[1]/newlk))<= eps && halving==0) { /* all done */
	    loglik[1] = newlk;
	    chinv2(imat, nvar);     /* invert the information matrix */
	    for (i=0; i<nvar; i++) {
		beta[i] = newbeta[i]*scale[i];
		u[i] /= scale[i];
		imat[i][i] *= scale[i]*scale[i];
		for (j=0; j<i; j++) {
		    imat[j][i] *= scale[i]*scale[j];
		    imat[i][j] = imat[j][i];
		    }
	    }
	    goto finish;
	}

	if (*iter== maxiter) break;  /*skip the step halving calc*/

	if (newlk < loglik[1])   {    /*it is not converging ! */
		halving =1;
		for (i=0; i<nvar; i++)
		    newbeta[i] = (newbeta[i] + beta[i]) /2; /*half of old increment */
		}
	else {
	    halving=0;
	    loglik[1] = newlk;
	    chsolve2(imat,nvar,u);
	    j=0;
	    for (i=0; i<nvar; i++) {
		beta[i] = newbeta[i];
		newbeta[i] = newbeta[i] +  u[i];
	        }
	    }
	}   /* return for another iteration */

    /*
    ** We end up here only if we ran out of iterations 
    */
    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=0; i<nvar; i++) {
	beta[i] = newbeta[i]*scale[i];
	u[i] /= scale[i];
	imat[i][i] *= scale[i]*scale[i];
	for (j=0; j<i; j++) {
	    imat[j][i] *= scale[i]*scale[j];
	    imat[i][j] = imat[j][i];
	    }
	}
    *flag = 1000;


finish:
   
    //return flag, loglik, beta. use vector as we will return to c++ land
    vector<double> result(nvar+2);
    result[0]=*flag;
    result[1]=loglik[1];
    for(int i=2; i<nvar;i++){
        result[i]=beta[i-2];
    }
    return result;
 }

