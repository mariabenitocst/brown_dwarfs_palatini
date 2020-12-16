// ====================================================================
// Calculate evolution of degeneracy parameter with time
//
// ====================================================================
#include <cstdio>
#include <iostream>
#include <cmath>
#include "Li2.hpp"
#include "complex.hpp"
#include <gsl/gsl_odeiv2.h>
#include <limits>
#include <cfloat>
using namespace std;

double Li2(double x) noexcept;

/**
// TODO check units of mu and mu_e!!
double b (const double xi,
          double T,
          double mu,
         )
{
    // return
    return (-5./16.*xi*log(1+exp(-mu/k_B/T)) + 
            15./8.*xi*xi*(M_PI**2/3. + Li2(2, -exp(-1/xi))));
}

int func (double t, 
          const double y[], 
          double f[],
          void *params)
{
  (void)(t); // avoid unused parameter warning
  double *_params = (double *) params;
  // get parameters from params
  double M     = _params[0];
  double T     = _params[1];
  double b1    = _params[2];
  double nu    = _params[3];
  double mu_1  = _params[4];
  double mu_e  = _params[5];
  double kR    = _params[6];
  double omega = _params[7];
  double gamma = _params[8];
  double delta = _params[9];
  double alpha = _params[10];

  double b     = _b()

  f[0] = y[1];

  // return
  return GSL_SUCCESS;
}
*/


int main()
{

    double x = Li2(1.);
    //cout << Li2(1.) << endl;

    // return 
    return 0;
}
