// Body for solving Lane-Emden equation
#include <sstream>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <vector>
using namespace std;

/***********************************************************************//**
 * @brief Equations
 ***************************************************************************/
int func (double xi, 
          const double y[], 
          double f[],
          void *params)
{
  double alpha = *(double *)params;

  double epsilon = 1e-5;
  xi = xi + epsilon;

  double factor = 1+2*alpha*pow(y[0], 1.5);
  f[0] = y[1];
  f[1] = - pow(y[0], 1.5)/factor - (2/xi + 3*alpha*xi/factor)*y[1] 
         - 3*alpha*sqrt(y[0])/factor*y[1]*y[1];
  // return
  return GSL_SUCCESS;
}


int
jac (double xi, 
     const double y[],
     double *dfdy,
     double dfdxi[],
     void *params)
{
  double alpha = *(double *)params;

  double epsilon = 0.001;
  xi = xi + epsilon;

  double factor = 1+2*alpha*pow(y[0], 1.5);
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, 3*sqrt(y[0])/factor*(-0.5 + alpha*pow(y[0], 1.5)/factor)
    + 9*alpha*alpha*sqrt(y[0])*xi*y[1]/factor/factor + 
    3*alpha/factor*(-0.5*pow(y[0], -0.5) + 
                    3*alpha*pow(y[0], 0.25)*y[1]*y[1]/factor));
  
  gsl_matrix_set (m, 1, 1, -2/xi + 3*alpha*xi/factor
    -6*alpha*sqrt(y[0])*y[1]/factor);
  dfdxi[0] = 0.0; 
  dfdxi[1] = 2./xi/xi - 3*alpha*y[1]/factor; // Is this interpretation correct?
  return GSL_SUCCESS;
}


vector<double> theta(double alpha)
{
  gsl_odeiv2_system sys = {func, jac, 2, &alpha};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  1e-6, 1e-6, 0.0);
  int i;
  double xi = 0.0, xi1 = 10.0;
  double y[2] = { 1.0, 0.0 };

  vector<double> _theta;

  for (i = 1; i <= 100; i++)
    {
      double xii = i * xi1 / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &xi, xii, y);

      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }
      _theta.push_back(y[0]);
      //printf ("%.5e %.5e %.5e\n", xi, y[0], y[1]);
    }
  gsl_odeiv2_driver_free (d);
  // return
  return _theta;
}

int main()
{

  double alpha = 0.1;
  vector<double> _theta;
  _theta = theta(alpha);

  // output file
  stringstream filename;
  filename << "../../data/theta_alpha=" << alpha << ".dat";
  FILE * outdata = fopen (filename.str().c_str(),"w"); 
  fprintf(outdata,"# xi   theta\n");

  int i;
  double xi1 = 10.0;

  for (i = 1; i <= 100; i++)
    {
      double xii = i * xi1 / 100.0;

      fprintf(outdata,"%.5e %.5e\n", xii, _theta[i]);
    }
  fclose (outdata);

  // return
  return 0;
}
