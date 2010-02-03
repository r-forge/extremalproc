#include "header.h"

// Composite log-likelihood for extremal models:

void CompLikelihood(int *corrmod, double *data, double *lags, int *model, 
		    int *ndata, int *nsite, double *par, double *res)
{
  int h=0, i=0, j=0, n=0; 
  double corr=0.0, df=0.0;

  df = par[0];
 
  for(i = 0; i < (*nsite - 1); i++)
    for(j = (i + 1); j < *nsite; j++)
      {
	corr = CorrelationFct(corrmod, lags[h], model, par);
	for(n = 0; n < *ndata; n++)
	  *res += PairLikelihood(df, corr, model, 
				 data[(n + i * *ndata)],
				 data[(n + j * *ndata)]);
	h++;
      }    

  if(!R_FINITE(*res))
    *res = LOW;

  return;
}

// Pairwise log-likelihood for extremal models on frechet scale:

double PairLikelihood(double df, double corr, int *model, double u, double v)
{
  double res=0.0, jac=0.0;

  jac = log(u * v * pow(log(u), 2) * pow(log(v), 2));

  switch(*model)
    {
    case 1:
      res = PairLikelihood_g(corr, u, v) + jac;
      break;
    case 2:
      res = PairLikelihood_t(df, corr, u, v) + jac;
      break;
    }
  
  return res;
}

// Pairwise log-likelihood for Extremal Gaussian model:

double PairLikelihood_g(double corr, double u, double v)
{
  double c=0.0, d2V=0.0, du=0.0, duV=0.0, dv=0.0, dvV=0.0;
  double theta2=0.0, ltheta=0.0, lnc=0.0, lnu=0.0, lnv=0.0;
  double pu=0.0, pv=0.0, res=0.0, theta=0.0, V=0.0, twu=0.0;
  double twv=0.0, w=0.0, z=0.0;
  
  theta = sqrt(corr);
  theta2 = 2 * theta;
  twu = theta2 * u;
  twv = theta2 * v;
  lnu = log(u);
  lnv = log(v);
  c = lnu / lnv;
  lnc = log(c);
  ltheta = lnc / theta2;
  z = theta + ltheta;
  w = theta - ltheta;
  pu = pnorm(z, 0, 1, 1, 0);
  pv = pnorm(w, 0, 1, 1, 0);
  du = dnorm(z, 0, 1, 0);
  dv = dnorm(w, 0, 1, 0);

  V = lnu * pu + lnv * pv;
  duV = (pu * theta2 + du - dv / c) / twu;
  dvV = (pv * theta2 + dv - du * c) / twv;
  d2V = (du * (z - theta2) + dv * (w - theta2) / c) / 
    (twu * twv * lnv);

  res = V + log(duV * dvV + d2V);  
 
  return res;
  }

// Pairwise log-likelihood for Extremal t model:

double PairLikelihood_t(double df,  double corr, double u, double v)
{
  double df1=0.0, df2=0.0;
  double lnu=0.0, lnv=0.0, c=0.0, cpdf=0.0, cmdf1=0.0, cmdf=0.0, a=0.0;
  double x=0.0, y=0.0, ptu=0.0, ptv=0.0, dtu=0.0, dtv=0.0;
  double denu=0.0, denv=0.0, denuv=0.0, denuvdf=0.0, duV=0.0;
  double dvV=0.0, d2Va=0.0, d2Vb=0.0, d2V=0.0, V=0.0, res=0.0;
  
  df1 = df + 1;
  df2 = pow(df, 2);

  lnu = log(u);
  lnv = log(v);
  c = lnu / lnv;
  cpdf = pow(c, 1 / df);
  cmdf = 1 / cpdf;
  cmdf1 = cmdf / c;
  a = sqrt((1 - pow(corr, 2)) / df1);

  x = (cpdf - corr) / a;
  y = (cmdf - corr) / a;

  ptu = pt(x, df1, 1, 0);
  ptv = pt(y, df1, 1, 0);
  dtu = dts(x, df1);
  dtv = dts(y, df1);
  denu = u * a;
  denv = v * a;
  denuv = denu * v;
  denuvdf = denuv * df2;
     
  V = lnu * ptu + lnv * ptv;

  duV = ptu / u + (dtu * cpdf - dtv * cmdf1) / (denu * df);

  dvV = ptv / v + (dtv * cmdf - dtu * cpdf * c) / (denv * df);

  d2Va = - cpdf * (dtu * df1 + d1x_dt(x, df1) * cpdf / a) /
    (denuvdf * lnv);
      
  d2Vb = - cmdf1 * (dtv * df1 + d1x_dt(y, df1) * cmdf / a) /
    (denuvdf * lnv);

  d2V = d2Va + d2Vb;

  res = V + log(duV * dvV + d2V);
 
  return res;
  }

