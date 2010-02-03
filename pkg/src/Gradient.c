#include "header.h"

// Bivariate gradient of Extremal Gaussian model:

void Gradient_g(double corr, double *gradcorr,  int ngrcor, 
		double u, double v, double *gradient)
{
  // Initialization variables:

  double c=0.0, d2V=0.0, du=0.0, duV=0.0, dv=0.0, dvV=0.0;
  double theta2=0.0, ltheta=0.0, lnc=0.0, lnu=0.0, lnv=0.0;
  double pu=0.0, pv=0.0, grad=0.0, theta=0.0, dV=0.0, twu=0.0;
  double twv=0.0, w=0.0, z=0.0, den=0.0, dtheta=0.0, dz=0.0;
  double dw=0.0, dduV=0.0, ddvV=0.0, dd2V=0.0;
  int i;
  
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

  duV = (pu * theta2 + du - dv / c) / twu;
  dvV = (pv * theta2 + dv - du * c) / twv;
  d2V = (du * (z - theta2) + dv * (w - theta2) / c) / 
    (twu * twv * lnv);

  den = duV * dvV + d2V;

  dtheta = 1 / theta2;
  dz = dtheta * (1 - lnc / theta2 / theta);
  dw = dtheta * (1 + lnc / theta2 / theta);

  dV = lnu * du * dz + lnv * dv * dw;

  dduV = (du * (dz * (theta2 - z) - dtheta / theta) + 
	  dv * (dtheta / theta + w * dw) / c) / twu;
 
  ddvV = (dv * (dw * (theta2 - w) - dtheta / theta) + 
	  du * (dtheta / theta + z * dz) * c) / twv;

  dd2V = - d2V * dtheta / theta + 
    (du * (dz - 2 * dtheta - z * (z - theta2)) + 
     dv * (dw - 2 * dtheta - w * (w - theta2)) / c) / 
    (twu * twv * lnv);

  grad = dV + (dduV * dvV + duV * ddvV + dd2V) / den;

  for(i = 0; i < ngrcor; i++)
    gradient[i] = grad * gradcorr[i];
 
  return;
  }

// Bivariate gradient of Extremal t model:

void Gradient_t(double corr, int flagc, double *gradcorr, double df, 
		int flagdf, int ngrcor, double u, double v, double *gradient)
{
  // Initialization variables:

  double c=0.0, df1=0.0, gradfcor=0.0, lu=0.0, lv=0.0, lc=0.0;
  double cpdf=0.0, cmdf=0.0, cpdf1=0.0, cmdf1=0.0, cpdfcorr1=0.0, cmdfcorr1=0.0;
  double a=0.0, udfa=0.0, vdfa=0.0, uvadflv=0.0, x=0.0, y=0.0, ptu=0.0, ptv=0.0;
  double dtu=0.0, dtv=0.0, d1tu=0.0, d1tv=0.0, duV=0.0, dvV=0.0, d2V=0.0;
  double den=0.0, a2=0.0, aomcorr2=0.0, dtucorr=0.0, dtvcorr=0.0, d1tua=0.0;
  double d1tva=0.0, d2tu=0.0, d2tv=0.0, dV=0.0, dduV=0.0, ddvV=0.0, dd2V=0.0;
  double dptu=0.0, dptv=0.0, ddtu=0.0, ddtv=0.0, df2=0.0, cpdflc=0.0;
  double cmdflc=0.0, df12=0.0, x2df=0.0, y2df=0.0, dd1tu=0.0, dd1tv=0.0, fx=0.0;
  double fy=0.0, fmlc=0.0, fmlcdf=0.0, fplc=0.0, fplcdf=0.0, ddfptu=0.0, ddfptv=0.0;
  double ddfdtu=0.0, ddfdtv=0.0, ddfd1tu=0.0, ddfd1tv=0.0, ddfV=0.0, ddfduV=0.0;
  double ddfdvV=0.0, ddfd2V=0.0;
  int i=0.0; 


  df1 = df + 1;
  lu = log(u);
  lv = log(v);
  c = lu / lv;
  lc = log(c);

  cpdf = pow(c, 1 / df);
  cmdf = 1 / cpdf;
  cpdf1 = cpdf * c;
  cmdf1 = cmdf / c;
  cpdfcorr1 = cpdf * corr - 1;
  cmdfcorr1 = cmdf * corr - 1;

  a = sqrt((1 - pow(corr, 2)) / df1);

  udfa = u * df * a;
  vdfa = v * df * a;
  uvadflv = v * df * udfa * lv;

  x = (cpdf - corr) / a;
  y = (cmdf - corr) / a;

  ptu = pt(x, df1, 1, 0);      // t distribution in x
  ptv = pt(y, df1, 1, 0);      // t distribution in y 
  dtu = dts(x, df1);           // t density in x
  dtv = dts(y, df1);           // t density in y
  d1tu = d1x_dt(x, df1);       // derivative of the t density respect u
  d1tv = d1x_dt(y, df1);       // derivative of the t density respect v

  // Start derivation density quantities
  duV = ptu / u + (dtu * cpdf - dtv * cmdf1) / udfa;
  dvV = ptv / v + (dtv * cmdf - dtu * cpdf1) / vdfa;

  d2V = - (cpdf * (dtu * df1 + d1tu * cpdf / a) + 
	   cmdf1 * (dtv * df1 + d1tv * cmdf / a)) / uvadflv;

  den = duV * dvV + d2V;
  //End derivation density quantities

  if(flagc > 0)
    {
      // variables definition for the computation
      // of the derivaties quantities respect with
      // the correlation parameter
      a2 = pow(a, 2);
      aomcorr2 = a * (1 - pow(corr, 2));
      dtucorr = dtu * corr;
      dtvcorr = dtv * corr;
      d1tua = d1tu / a;
      d1tva = d1tv / a;
      d2tu = d2x_dt(x, df1);
      d2tv = d2x_dt(y, df1);
      
      //Start derivatives respect with the correlation parameter
      dV = (lu * dtu * cpdfcorr1 + lv * dtv * cmdfcorr1) / aomcorr2;

      dduV = ((cpdf * (dtucorr + d1tua * cpdfcorr1) - 
	       cmdf1 * (dtvcorr + d1tva * cmdfcorr1)) / df + dtu * cpdfcorr1) / 
	(u * aomcorr2);

      ddvV = ((cmdf * (dtvcorr + d1tva * cmdfcorr1) - 
	       cpdf1 * (dtucorr + d1tua * cpdfcorr1)) / df + dtv * cmdfcorr1) / 
	(v * aomcorr2);

      dd2V = - (cpdf * ((dtucorr + d1tua * cpdfcorr1) / a2 + 
			cpdf * (2 * d1tu * corr + d2tu * cpdfcorr1 / a) / aomcorr2) + 
		cmdf1 * ((dtvcorr + d1tva * cmdfcorr1) / a2 + 
			 cmdf * (2 * d1tv * corr + d2tv * cmdfcorr1/ a) / aomcorr2)) / 
	uvadflv;

      gradfcor = dV + (dduV * dvV + duV * ddvV + dd2V) / den;
      for(i = 0; i < ngrcor; i++)
	gradient[i + flagdf] = gradfcor * gradcorr[i];
      // End derivatives respect with the correlation parameter
    }

  if(flagdf == 1)
    {
      // variables definition for the computation
      // of the derivaties quantities respect with
      // the degree of freedom parameter

      dptu = ddf_pt(x, df1);
      dptv = ddf_pt(y, df1);
      ddtu = ddf_dt(x, df1);
      ddtv = ddf_dt(y, df1);
      df2 = pow(df, 2);
      cpdflc = cpdf * lc;
      cmdflc = cmdf * lc;
      dd1tu = ddf_t_d1x_dt(x, df1, a, cpdflc);
      dd1tv = ddf_t_d1x_dt(y, df1, a, -cmdflc);
      df12 = 2 * df1;
      x2df = 1 + pow(x, 2) / df1;
      y2df = 1 + pow(y, 2) / df1;
      fx = x / df12 - cpdflc / (a * df2);
      fy = y / df12 + cmdflc / (a * df2);
      fmlc = 1 / df12 - lc / df2;
      fmlcdf = fmlc - 1 / df;
      fplc = 1 / df12 + lc / df2;
      fplcdf = fplc - 1 / df;
      ddfptu = dptu + dtu * fx;
      ddfptv = dptv + dtv * fy;
      ddfdtu = ddtu - dtu * (df1 + 1) * x / df1 * fx / x2df;
      ddfdtv = ddtv - dtv * (df1 + 1) * y / df1 * fy / y2df;
      ddfd1tu = ddfdtu + dtu * fmlcdf;
      ddfd1tv = ddfdtv + dtv * fplcdf;

 
      // Start derivatives respect with the degree of freedom parameter
      ddfV = lu * ddfptu  + lv * ddfptv;

      ddfduV = ddfptu / u + (cpdf * ddfd1tu - cmdf1 * ddfd1tv) / udfa;

      ddfdvV = ddfptv / v + (cmdf * ddfd1tv - cpdf1 * ddfd1tu) / vdfa;

      ddfd2V = - (cpdf * ((fmlc - 2 / df) * (dtu * df1 + d1tu * cpdf / a) + 
			  (dtu + df1 * ddfdtu + cpdf * (dd1tu + d1tu * fmlc) / a)) + 
		  cmdf1 * ((fplc - 2 / df) * (dtv * df1 + d1tv * cmdf / a) + 
			   (dtv + df1 * ddfdtv + cmdf * (dd1tv + d1tv * fplc) / a))) / 
	uvadflv;
      
      gradient[0] = ddfV + (ddfduV * dvV + duV * ddfdvV + ddfd2V) / den;
    
      // End derivatives respect with the degree of freedom parameter
    }

  return;
}

void SquaredScore(int *corrmod, double *data, double *eps, int *flag, double *lags, int *model, 
		  int *ndata, int *nflag, int *nsite, int *nsqsc, double *par, double *sqscore)
{
  int flagc=0, h=0, i=0, j=0, m=0, n=0, ngcorr=0, npair=0;
  double corr, df, *grad, *gradient, *grc;

  df = par[0];
  ngcorr = *nsqsc;
  if(flag[0] == 1) ngcorr = ngcorr - 1;

  for(i = 1; i < *nflag; i++)
    flagc = flagc + flag[i];

  npair = *nsite * (*nsite - 1) / 2;
  grc = (double *) R_alloc(ngcorr, sizeof(double));
  grad = (double *) R_alloc(*nsqsc, sizeof(double));
  gradient = (double *) R_alloc(*nsqsc, sizeof(double));

  for(i = 0; i < *nsqsc; i++)
    for(j = 0; j < *nsqsc; j++)
      sqscore[i * *nsqsc + j] = 0;

  for(n = 0; n < *ndata; n++)
    {
      h = 0;
      for(i = 0; i < *nsqsc; i++)
        gradient[i] = 0;

      for(i = 0; i < (*nsite - 1); i++)
        for(j = (i + 1); j < *nsite; j++)
	  {
	    corr = CorrelationFct(corrmod, lags[h], model, par);
	    GradientCorrFct(corr, corrmod, eps, flag, grc, lags[h], model, par);
	    switch(*model)
	      {
	      case 1:// Extremal Gaussian process
		Gradient_g(corr, grc, ngcorr, data[(n + i * *ndata)], 
			   data[(n + j * *ndata)], grad);
		break;
	      case 2:// Extremal t process
		Gradient_t(corr, flagc, grc, df, flag[0], ngcorr, 
			   data[(n + i * *ndata)], data[(n + j * *ndata)], grad);
		break;
	      }

	    for(m = 0; m < *nsqsc; m++)
	      gradient[m] = gradient[m] + grad[m];

	    h++;
	  }

      for(i = 0; i < *nsqsc; i++)
	for(j = 0; j < *nsqsc; j++)
	  sqscore[i * *nsqsc + j] = sqscore[i * *nsqsc + j] + 
	    gradient[i] * gradient[j];

    }

  return;
}
