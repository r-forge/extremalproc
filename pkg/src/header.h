#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>

#define LOW -1.0e8


// 1)
/*----------------------------------------------------------------
File name: CorrelationFunction.c
Description: procedures for computation of correlation functions
Start
 ---------------------------------------------------------------*/

double CorrelationFct(int *corrmod, double lag, int *model, double *par);

void GradientCorrFct(double corr, int *corrmod, double *eps, int *flag, double *grad, 
		     double lag, int *model, double *par);


/*----------------------------------------------------------------
File name: CorrelationFunction.c
Description: procedures for computation of correlation functions
End
 ---------------------------------------------------------------*/

// 2)
/*----------------------------------------------------------------
File name: CompositeLikelihood.c
Description: functions for composite log-likelihood evaluation
Start
 ---------------------------------------------------------------*/

void CompLikelihood(int *corrmod, double *data, double *lags, int *model, 
		    int *ndata, int *nsite, double *par, double *res);

double PairLikelihood(double df, double corr, int *model, double u, double v);

double PairLikelihood_g(double corr, double u, double v);

double PairLikelihood_t(double df, double corr, double u, double v);


/*----------------------------------------------------------------
File name: CompositeLikelihood.c
Description: functions for composite log-likelihood evaluation
End
 ---------------------------------------------------------------*/

// 3)
/*----------------------------------------------------------------
File name: Gradient.c
Description: procedures for computation of gradients quantities
Start
 ---------------------------------------------------------------*/

void Gradient(int *corrmod, double *data, double *eps, int *flag, 
	      double *lags, int *model, int *ndata, int *nflag, 
	      int *npar, int *nsite, double *par, double *gradient);

void Gradient_g(double corr, double *gradcorr,  int ngrcor, 
		double u, double v, double *gradient);

void Gradient_t(double corr, int flagc, double *gradcorr, double df, 
		int flagdf, int ngrcor, double u, double v, double *gradient);

void SquaredScore(int *corrmod, double *data, double *eps, int *flag, double *lags, 
		  int *model, int *ndata, int *nflag, int *npar, int *nsite,  
		  double *par, double *varmat);


/*----------------------------------------------------------------
File name: Gradient.c
Description: procedures for computation of gradients quantities
End
 ---------------------------------------------------------------*/

// 4)
/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for computation of different quatities
Start
 ---------------------------------------------------------------*/

void ComputeMaxima(double *df, double *maxima, int *model, int *nblock, 
		   int *nsite, double *simu);

double dgev(double x, double loc, double scale, double shape);

double d1x_dt(double x, double df);

double d2x_dt(double x, double df);

double ddf_pt(double x, double df);

double ddf_dt(double x, double df);

double ddf_d1x_dt(double x, double df);

double ddf_t_d1x_dt(double x, double df, double a, double fc);

double dts(double x, double df);

double frechet2gev(double x, double alpha, double beta, double gamma,
		   double loc, double scale, double shape);

void FromDistToDist(double *data, double *alpha, double *beta, double *gamma, 
		    double *loc, int *ndata, int *nsite, double *scale, 
		    double *shape, int *type, double *res);

double gev2gumbel(double x, double loc, double scale, double shape);

double gev2unitfrechet(double x, double loc, double scale, double shape);

void GevLogLik(double *data, double *loc, int *ndata, double *scale,
	       double *shape, double *res);

double gumbel2gev(double x, double loc, double scale, double shape);

double int_pt(double x, double df);

void integr_pt(double *x, int n, void *ex);

double pgev(double x, double loc, double scale, double shape);

double qgev(double x, double loc, double scale, double shape);

double qgumbel(double x, double loc, double scale);

double ProbIntTrans(double data, double alpha, double beta, double gamma, 
		    double loc, double scale, double shape, int type);

double unitfrechet2gev(double x, double loc, double scale, double shape);

/*----------------------------------------------------------------
File name: Utility.c
Description: procedures for computation of different quatities
End
 ---------------------------------------------------------------*/
