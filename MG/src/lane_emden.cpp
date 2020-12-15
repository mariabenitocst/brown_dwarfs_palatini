// Body for solving Lane-Emden equation
#include "gsl_odeiv2.h"
#include <cmath>


/***********************************************************************//**
 * @brief Equations
 ***************************************************************************/
int func (double xi, 
          const double theta[], 
          double y[],
          void *params)
{
  double alpha = *(double *)params;
  double factor = 1+2*alpha*pow(theta[0], 1.5)
  y[0] = theta[1];
  y[1] = - pow(theta[0], 1.5)/factor - (2/xi + 3*alpha*xi/factor)*theta[1] 
         - 3*alpha*sqrt(theta[0])/factor*theta[1]*theta[1];
  // return
  return GSL_SUCCESS;
}


int main()
{


    // return
    return 0;
}
