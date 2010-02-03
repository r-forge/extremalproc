#include "header.h"

double dts(double x, double df)
{
  double res=0.0;

  res = gamma((df + 1) / 2) / (sqrt(df * M_PI) * gamma(df / 2)) *
    pow(1 + pow(x, 2) / df, - (df + 1) / 2);
  
  return res;
}

double d1x_dt(double x, double df)
{
  double res=0.0;

  res = - dts(x, df) * (df + 1) * x /
    (1 + pow(x, 2) / df) / df;

  return res;
}

double d2x_dt(double x, double df)
  {
    double df1=0.0, res=0.0, y=0.0;

    df1 = (df + 1) / df;
    y = 1 + pow(x, 2) / df;

    res = - d1x_dt(x, df) * df1 * x / y - 
      dts(x, df) * df1 * (1 / y - 2 * pow(x, 2) / (df * pow(y, 2)));

      return res;
  }

double int_pt(double x, double df)
  {
    double res=0.0, x2=0.0, y=0.0;

    x2 = pow(x, 2);
    y = 1 + x2 / df;
      
    res = dts(x, df) * ((df + 1) * x2 / pow(df, 2) / y - log(y)) / 2;

    return res;
  }

void integr_pt(double *x, int n, void *ex)
{
  int i;
  double *df, d=0.0;

  *df = *((double *) ex);
  d = *df;

  for(i = 0; i < n; i++)
    x[i] = int_pt(x[i], d);

  return;
}

double ddf_pt(double x, double df)
  {
    double epsabs=0.0, epsrel=0.0, integ=0.0, tinteg=0.0;
    double abserr=0.0, origin=0.0, q=0.0, res=0.0, *work;
    int inf=0.0, neval=0.0, ier=0.0, limit=0.0, lenw=0.0;
    int last=0.0, *iwork;

    inf = -1;
    epsabs = 1e-5;
    epsrel = 1e-5;
    limit = 100;
    lenw = 4 * limit;
    iwork = (int *) R_alloc(limit, sizeof(int));
    work = (double *) R_alloc(lenw, sizeof(double));

    if(x <= 0)
      {
	Rdqagi(integr_pt, (void*)&df, &x, &inf, &epsabs, &epsrel, 
	       &integ, &abserr, &neval, &ier, &limit, &lenw, &last, 
	       iwork, work);
      }
    else
      {
	q = - x;
	Rdqagi(integr_pt, (void*)&df, &origin, &inf, &epsabs, &epsrel, 
	       &tinteg, &abserr, &neval, &ier, &limit, &lenw, &last, 
	       iwork, work);
	Rdqagi(integr_pt, (void*)&df, &q, &inf, &epsabs, &epsrel, 
	       &integ, &abserr, &neval, &ier, &limit, &lenw, &last, 
	       iwork, work);
	integ = 2 * tinteg - integ;
	}

    res =  pt(x, df, 1, 0) * (digamma((df + 1) / 2) - 
			      digamma(df / 2)  - 1 / df) / 2 + integ;

    return res;
  }

double ddf_dt(double x, double df)
  {
    double df1=0.0, res=0.0, x2=0.0, y=0.0;

    df1 = df + 1;
    x2 = pow(x, 2);
    y = 1 + x2 / df;

    res = dts(x, df) * (digamma(df1 / 2) - digamma(df / 2) - 1 / df + 
			df1 * x2 / pow(df, 2) / y - log(y)) / 2;

    return res;
  }

double ddf_d1x_dt(double x, double df)
  {
    double df1=0.0, df2=0.0, res=0.0, x2=0.0, xdf=0.0, y=0.0;

    df1 = df + 1;
    df2 =  pow(df, 2);
    x2 = pow(x, 2);
    y = 1 + x2 / df;
    xdf = x / df;

    res = - dts(x, df) * xdf * 
      (df1 / 2 * (digamma(df1 / 2) + pow(xdf, 2) * (df + 3) / y - 
		  digamma(df / 2) - log(y) - 3 / df) + 1) / y;

      return res;
  }

double ddf_t_d1x_dt(double x, double df, double a, double fc)
  {
    double df1=0.0, df2=0.0, res=0.0, x2=0.0, y=0.0;

    df1 = df + 1;
    df2 =  pow(df - 1, 2);
    x2 = pow(x, 2);
    y = 1 + x2 / df;
    
    res = ddf_d1x_dt(x, df) - dts(x, df) * df1 * 
      (fc / a / df2 * (x2 * (df1 + 2) / df / y - 1) + x / 2 / df - 
       pow(x, 3) * (df1 + 2) / 2 / pow(df, 2) / y) / df / y;

    return res;
  }

void ComputeMaxima(double *df, double *maxima, int *model, int *nblock, 
		   int *nsite, double *simu)
{
  int i=0, k=0;
  double an=0.0, bn=0.0, chi2=1.0, n=0.0;

  // Seeting: normalizing constants for the Student t
  // and Gaussian model
  n = (double) *nblock;

  /* switch(*model)
    {
    case 1:
      an = qt(1 - 1 / n, *df, 1, 0);
      bn = 0;
      // First loop: number of blocks
      for(i = 0; i < *nblock; i++)
	{
	  chi2 = sqrt(rchisq(*df) / *df);
	  // Second loop: number of sites
	  for(k = 0; k < *nsite; k++)
	    {
	      simu[k + i * *nsite] = simu[k + i * *nsite] / chi2;
	      maxima[k] = fmax(maxima[k], (simu[k + i * *nsite] - bn) / an);
	    }
	}
      break;
    case 2:
      bn = sqrt(2 * log(n)) - 
	(0.5 * log(log(n)) + log(2 * sqrt(M_PI))) / 
	sqrt(2 * log(n));
      an = 1 / bn;
      // First loop: number of blocks
      for(i = 0; i < *nblock; i++)
	{
	  // Second loop: number of sites
	  for(k = 0; k < *nsite; k++)
	    maxima[k] = fmax(maxima[k], simu[k + i * *nsite]);
	  if(i == (*nblock - 1))
	    maxima[k] = (maxima[k] - bn) / an;
	}
	break;
	}*/

  if(*model==1)
    {
       bn = sqrt(2 * log(n)) - 
	(0.5 * log(log(n)) + log(2 * sqrt(M_PI))) / 
	 sqrt(2 * log(n));
       an = 1 / bn;
       // First loop: number of blocks  
       for(i = 0; i < *nblock; i++)
	 {
	   // Second loop: number of sites
	   for(k = 0; k < *nsite; k++)
	     { 
	       // Compute the componentwise maxima
	       maxima[k] = fmax(maxima[k], simu[k + i * *nsite]);
	       if(i == (*nblock - 1))
		 maxima[k] = (maxima[k] - bn) / an;
	     }
	 }
    }

  if(*model==2)
    {
      an = qt(1 - 1 / n, *df, 1, 0);
      bn = 0;
      // First loop: number of blocks  
      for(i = 0; i < *nblock; i++)
	{
	  chi2 = sqrt(rchisq(*df) / *df);
	  // Second loop: number of sites
	  for(k = 0; k < *nsite; k++)
	    { 
	      // Compute the componentwise maxima
	      maxima[k] = fmax(maxima[k], simu[k + i * *nsite] / chi2);
	      if(i == (*nblock - 1))
		maxima[k] = (maxima[k] - bn) / an;
	    }
	}
    }


      // First loop: number of blocks
  /*     for(i = 0; i < *nblock; i++)
	{
	  // Second loop: number of sites
	  for(k = 0; k < *nsite; k++)
	    maxima[k] = fmax(maxima[k], simu[k + i * *nsite] - bn);
	  if(i == (*nblock - 1))
	    maxima[k] = (maxima[k] - bn) / an;
	    }*/


    // First loop: number of blocks  
  /*  for(i = 0; i < *nblock; i++)
    {
      if(*model == 1)
	chi2 = sqrt(rchisq(*df) / *df);
      // Second loop: number of sites
      for(k = 0; k < *nsite; k++)
	{ // Compute the componentwise maxima
	  maxima[k] = fmax(maxima[k], simu[k + i * *nsite] / chi2);
	  if(i == (*nblock - 1))
	    maxima[k] = (maxima[k] - bn) / an;
	}
	}*/
  return;
}

double dgev(double x, double loc, double scale, double shape)
{
  double y=0.0, result=0.0;

  y = (x - loc) / scale;

  if(shape==0)
    result = exp(-exp(-y) - y) / scale;
  else
    result = exp(-pow(fmax(1 + shape * y, 0), - 1 / shape)) * 
      pow(fmax(1 + shape * y, 0), - 1 / shape - 1) / scale;

  return result;
}

double qgumbel(double x, double loc, double scale)
{
  double result=0.0;

  result = loc - scale * log(-log(x));

  return result;
}


double qgev(double x, double loc, double scale, double shape)
{
  double result=0.0;

  if(shape==0)
    result = qgumbel(x, loc, scale);
  else
    result = loc + scale * (pow(-log(x), -shape) - 1) / shape;

  return result;
}


double frechet2gev(double x, double alpha, double beta, double gamma,
		   double loc, double scale, double shape)
{
  double u=0.0, result=0.0;

  u = pgev(x, alpha, beta, gamma);

  result = qgev(u, loc, scale, shape);

  return result;

}

void GevLogLik(double *data, double *loc, int *ndata, double *scale,
	       double *shape, double *res)
{
  int n;
  
  if((*scale <= 0) || (*shape < -1))
    {
      *res = LOW;
      return;
    }

  for(n = 0; n < *ndata; n++)
    *res += log(dgev(data[n], *loc, *scale, *shape));

  return;
}


double gumbel2gev(double x, double loc, double scale, double shape)
{
  double y=0.0, result=0.0;

  y = exp(x);

  result = unitfrechet2gev(y, loc, scale, shape);

  return result;

}


void FromDistToDist(double *data, double *alpha, double *beta, double *gamma, 
		    double *loc, int *ndata, int *nsite, double *scale, 
		    double *shape, int *type, double *res)
{
  int i=0, k=0;

  for(k = 0; k < *nsite; k++)
    for(i = 0; i < *ndata; i++)
      res[i + k * *ndata] = ProbIntTrans(data[i + k * *ndata], alpha[k], 
					 beta[k], gamma[k], loc[k], 
					 scale[k], shape[k], *type);
   
  return;
}


double ProbIntTrans(double data, double alpha, double beta, double gamma, 
		    double loc, double scale, double shape, int type)
{
  double result;

  switch(type)
    {
    case 1: // Transform from GEV to UNIFORM
      result = pgev(data, loc, scale, shape);
      break;
    case 2: // Transform from GEV to unit FRECHET
      result = gev2unitfrechet(data, loc, scale, shape);
      break;
    case 3: // Transform from FRECHET to GEV
      result = frechet2gev(data, alpha, beta, gamma, 
			   loc, scale, shape);
      break;
    case 4: // Transform from GUMBEL to GEV
      result = gumbel2gev(data, loc, scale, shape);
      break;
    case 5: // Transform from GEV to unit GUMBEL
      result = gev2gumbel(data, loc, scale, shape);
      break;
    }
   
  return result;
}

double gev2gumbel(double x, double loc, double scale, double shape)
{
  double u=0.0, result=0.0;

  u = pgev(x, loc, scale, shape);
  // Transform to unit Gumbel
  result = qgumbel(u, 0, 1);

  return result;
}

double gev2unitfrechet(double x, double loc, double scale, double shape)
{
  double y=0.0, result=0.0;

  y = (x - loc) / scale;

  result = pow(fmax(1 + shape * y, 0), shape);

  return result;

}

double pgev(double x, double loc, double scale, double shape)
{
  double y=0.0, result=0.0;

  y = (x - loc) / scale;

  if(shape==0)
    result = exp(-exp(-y));
  else
    result = exp(-pow(fmax(1 + shape * y, 0), - 1 / shape));

  return result;
}

double unitfrechet2gev(double x, double loc, double scale, double shape)
{
  double result=0.0;

  result = scale * (pow(x, 1 / shape) - 1) /  shape + loc;

  return result;

}
