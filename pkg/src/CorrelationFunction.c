#include "header.h"

double CorrelationFct(int *corrmod, double lag, int *model, double *par)
{
  double corr=0.0, nugget=0.0, power=0.0;
  double power1=0.0, power2=0.0, scale=0.0, smooth=0.0;

  nugget = par[1];

  if(*model == 1)
    {
      switch(*corrmod) // Extremal dependence functions are in alphabetical order
	{
          case 2:// Exponential correlation function
            scale = par[2];
            corr = (1 - nugget) * lag / scale;
            break;
          case 3:// Gaussian correlation function
            scale = par[2];
            corr = (1 - nugget) * pow(lag / scale, 2);
            break;
          case 4: // Generalised Cuachy correlation function 
            power1 = par[2];
            power2 = par[3];
            scale = par[4];
            corr = (1 - nugget) * power2 / power1 * pow(lag / scale, power1);
            break;
          case 5:// Stable correlation function
            power = par[2];
            scale = par[3];
            corr = (1 - nugget) * pow(lag / scale, power);
            break;
	}
    }

  if(*model == 2)
    {
      switch(*corrmod) // Correlation functions are in alphabetical order
	{
          case 1:// Cauchy correlation function
            power2 = par[2];
            scale = par[3];
            corr = (1 - nugget) * pow((1 + pow(lag / scale, 2)), - power2);
            break;
          case 2:// Exponential correlation function
            scale = par[2];
            corr = (1 - nugget) * exp(- lag / scale);
            break;
          case 3:// Gaussian correlation function
            scale = par[2];
            corr = (1 - nugget) * exp(-pow(lag / scale, 2));
            break;
          case 4: // Generalised Cuachy correlation function 
            power1 = par[2];
            power2 = par[3];
            scale = par[4];
            corr = (1 - nugget) * pow((1 + pow(lag / scale, power1)), - power2 / power1);
            break;
          case 5:// Stable correlation function
            power = par[2];
            scale = par[3];
            corr = (1 - nugget) * exp(-pow(lag / scale, power));
            break;
          case 6://  Whittle-Matern correlation function
            scale = par[2];
            smooth = par[3];
            corr = (1 - nugget) * pow(2, 1 - smooth) / gamma(smooth) * 
	      pow(lag / scale, smooth) * bessel_k(lag / scale, smooth, 1);
            break;
	}
    }

  return corr;
}

void GradientCorrFct(double corr, int *corrmod, double *eps, int *flag, double *grad, 
		     double lag, int *model, double *par)
{
  int i=0;
  double nugget=0.0, power=0.0, power1=0.0, power2=0.0, scale=0.0, smooth=0.0;
  double parscale=0.0, parsmooth=0.0;

  nugget = par[1];
  if(flag[1] == 1)
    {
      grad[i] = - corr;
      i++;
    }

  if(*model == 1)
    {
      switch(*corrmod)// Extremal dependence functions are in alphabetical order
	{
	case 2:// Exponential correlation function
          if(flag[2] == 1)
	    {
	      scale = par[2];
	      grad[i] = - corr / scale;
	      i++;
	    }
	  break;
	case 3:// Gaussian correlation function
          if(flag[2] == 1)
	    {
	      scale = par[2];
	      grad[i] = - 2 * pow(corr, 2) / scale;
	    }
	  break;
	case 4:// Generalised Cuachy correlation function
	  power1 = par[2];
	  power2 = par[3];
	  scale = par[4];
	  if(flag[2] == 1)
	    {
	      grad[i] = corr * (log(lag / scale) - 1 / power1);
	      i++;
	    }
	  if(flag[4] == 1)
	    {
	      grad[i] = power2 / scale * pow(lag / scale, power1 - 1);
	    }
	  break;
	case 5:// Stable correlation function
	  power = par[2];
	  scale = par[3];
	  if(flag[2] == 1)
	    {
	      grad[i] = corr * log(lag / scale);
	      i++;
	    }
	  if(flag[3] == 1)
	    {
	      grad[i] = - corr * power / scale;
	    }
	  break;
	}
    }

  if(*model == 2)
    {
      switch(*corrmod)// Correlation functions are in alphabetical order
	{
	case 1:// Cauchy correlation function
          power2 = par[2];
	  scale = par[3];
	  if(flag[2] == 1)
	    {
	      grad[i] = - corr * log(pow(corr, - 1 / power2));
	      i++;
	    }
	  if(flag[3] == 1)
	    {
	      grad[i] = 2 * power2 * corr * pow(corr, 1 / power2) * 
		pow(lag, 2) / pow(scale, 3);
	    }
	  break;
	case 2:// Exponential correlation function
          if(flag[2] == 1)
	    {
	      scale = par[2];
	      grad[i] = corr * lag / pow(scale, 2);
	      i++;
	    }
	  break;
	case 3:// Gaussian correlation function
          if(flag[2] == 1)
	    {
	      scale = par[2];
	      grad[i] = 2 * corr * pow(lag, 2) / pow(scale, 3);
	    }
	  break;
	case 4:// Generalised Cuachy correlation function
	  power1 = par[2];
	  power2 = par[3];
	  scale = par[4];
	  if(flag[2] == 1)
	    {
	      grad[i] = power2 * corr / power1 * (log(1 + pow(lag / scale, power1)) / 
						  power1 - pow(lag / scale, power1) * 
						  log(lag / scale) / (1 + pow(lag/ scale, power1)));
	      i++;
	    }
	  if(flag[3] == 1)
	    {
	      grad[i] = - corr * log(1 + pow(lag / scale, power1)) / power1;
	      i++;
	    }
	  if(flag[4] == 1)
	    {
	      grad[i] = corr / (1 + pow(lag / scale, 2)) * power2 * 
		pow(lag, power1) / pow(scale, power1 + 1);
	    }
	  break;
	case 5:// Stable correlation function
	  power = par[2];
	  scale = par[3];
	  if(flag[2] == 1)
	    {
	      grad[i] = - corr * pow(lag / scale, power) * log(lag / scale);
	      i++;
	    }
	  if(flag[3] == 1)
	    {
	      grad[i] = corr * pow(lag / scale, power - 1) * 
		power * lag / pow(scale, 2);
	    }
	  break;
	case 6:// Whittle-Matern correlation function
	  scale = par[2];
	  smooth = par[3];
	  if(flag[2] == 1)
	    {
	      parscale = (bessel_k(lag / (scale + *eps), smooth, 1) - 
			  bessel_k(lag / scale, smooth, 1)) / *eps;
	      grad[i] = (1 - nugget) * pow(2, 1 - smooth) / gamma(smooth) * 
		pow(lag / scale, smooth) * (parscale - smooth * 
					    bessel_k(lag / scale, smooth, 1) / scale);

	      i++;
	    }
	  if(flag[3] == 1)
	    {
	      parsmooth = (bessel_k(lag / scale, smooth + *eps, 1) - 
			   bessel_k(lag / scale, smooth, 1)) / *eps;
	      grad[i] = (1 - nugget) * pow(2, 1 - smooth) * 
		pow(lag / scale, smooth) / gamma(smooth) * 
		(log(lag / scale) - log(2) - digamma(smooth) * 
		 bessel_k(lag / scale, smooth, 1) + parsmooth);
	    }
	  break;
	}
    }

  return;
}


/*int expo(int *K, double *h, double *par, double *rho, double *nu)
{
  double sill, scale;
  int i, pcond;

  sill = par[0];
  scale = par[1];
  *nu = par[2];
  
  if((sill <= 0 || sill > 1) || (scale <= 0) || (*nu < 0))
    pcond = 0;
  else
    {
      pcond = 1;
      for(i = 0; i < *K; i++)
	rho[i] = sill * exp(- h[i] / scale);
    }
  return pcond;
}

int Gauss(int *K, double *h, double *par, double *rho, double *nu)
{
  double sill, scale;
  int i, pcond;

  sill = par[0];
  scale = par[1];
  *nu = par[2];
  
  if((sill <= 0 || sill > 1) || (scale <= 0) || (*nu < 0))
    pcond = 0;
  else
    {
      pcond = 1;
      for(i = 0; i < *K; i++)
	rho[i] = sill * exp(-pow(h[i] / scale, 2));
    }

  return pcond;
}

int stable(int *K, double *h, double *par, double *rho, double *nu)
{
  double sill, smooth, scale;
  int i, pcond;

  sill = par[0];
  scale = par[1];
  smooth = par[2];
  *nu = par[3];

  if((sill <= 0 || sill > 1) || 
     (smooth <= 0 || smooth > 2) || (scale <= 0) || (*nu < 0))
    pcond = 0;
  else
    {
      pcond = 1;
      for(i = 0; i < *K; i++)
	rho[i] = sill * exp(-pow(h[i] / scale, smooth));
    }

  return pcond;
}

double Cauchy(double dist, double range, double sill, double smooth)
{
  double rho;

  rho = sill * pow((1 + pow(dist / range, 2)), - smooth);

  return rho;
}

int genCauchy(int *K, double *h, double *par, double *rho, double *nu)
{
  double sill, smooth, power, scale;
  int i, pcond;

  sill = par[0];
  scale = par[1];
  smooth = par[2];
  power = par[3];
  *nu = par[4];

  if((sill <= 0 || sill > 1) || (smooth <= 0 || smooth > 2) || 
     (power < 0) || (scale <= 0) || (*nu < 0))
    pcond = 0;
  else
    {
      pcond = 1;
      for(i = 0; i < *K; i++)
	rho[i] = sill * pow((1 + pow(h[i] / scale, smooth)), -power / smooth);
    }

  return pcond;
}

int whittleMatern(int *K, double *h, double *par, double *rho, double *nu)
{
  double sill, smooth, scale;
  int i, pcond;

  sill = par[0];
  scale = par[1];
  smooth = par[2];
  *nu = par[3];

  if((sill <= 0 || sill > 1) || (smooth <= 0) || (scale <= 0) || (*nu < 0))
    pcond = 0;
  else
    {
      pcond = 1;
      for(i = 0; i < *K; i++)
	rho[i] = sill * pow(2, 1 - smooth) / gammafn(smooth) * pow(h[i] / scale, smooth) * bessel_k(h[i] / scale, smooth, 1);
    }

  return pcond;
  }*/
