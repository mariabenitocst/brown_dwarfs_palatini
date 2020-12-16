// ====================================================================
// Calculate evolution of degeneracy parameter with time
//
// ====================================================================
#include <sstream>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_dilog.h>
using namespace std;

/***********************************************************************//**
 * @brief Polylogarithm equation
 ***************************************************************************/
double Li2(double x) 
{
   gsl_sf_result li2_gsl{};
   gsl_sf_dilog_e(x, &li2_gsl);
   // return
   return li2_gsl.val;
}

// TODO check units of mu and mu_e!!
double _b (const double psi,
           double T,
           double mu)
{
    double k_B = 1.380649e-23; // Boltzmann constant [J/K]
    // return
    return (-5./16.*psi*log(1+exp(-mu/k_B/T)) + 
            15./8.*psi*psi*(pow(M_PI, 2/3.) + Li2(-exp(-1/psi))));
}

/***********************************************************************//**
 * @brief Function to be integrated
 ***************************************************************************/
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

  double b     = _b(y[0], T, mu_e);
  double a     = 2.5*mu_e/mu_1;

  double f_1 = 1.1634e-18*pow(b1, 2.856)*mu_1/(pow(kR, 1.1424)*pow(mu_e, 8./3.));
  double f_2 = pow(gamma, 0.7143)*pow(1-1.33*alpha/delta, 1.143)/omega;

  f[0] = f_1*f_2*pow(M, -1.094)*pow(y[0], 2.856*nu)*pow(1+b+a*y[0], 1.715);

  // return
  return GSL_SUCCESS;
}


int main()
{
    // ask user for alpha, delta and gamma values
    double alpha, gamma, delta;
    cout << "Enter alpha value: ";
    cin >> alpha;
    cout << "Enter gamma value: ";
    cin >> gamma;
    cout << "Enter delta value: ";
    cin >> delta;

    // output file
    stringstream filename;
    filename << "../../data/evolution_degeneracy_alpha=" << alpha << ".dat";
    FILE * outdata = fopen (filename.str().c_str(),"w"); 
    fprintf(outdata,"# time [year]  degeneracy");

    size_t dim = 1;
    double params[11];
    params[0]  = 0.03; // mass [Msun]
    params[1]  = pow(10, 3.94);
    params[2]  = 2.;
    params[3]  = 1.60;
    double X   = 0.75; // mass fraction H
    double Y   = 0.25; // mass fracion He
    // mean molecual weight for He and ionized H mixture
    params[4]  = pow((1+0.5*0.51)*X+Y/4., -1);
    params[5]  = X + 0.5*Y; // number of baryons per electron
    params[6]  = 0.01; // Rossland opacity [cm2/g]
    params[7]  = 1.; // Omega
    params[8]  = gamma; // gamma
    params[9]  = delta; // delta
    params[10] = alpha; //alpha
    // declare ODE system
    gsl_odeiv2_system sys = {func, NULL, dim, &params};

    double h       = 1.e-6; // starter step size for OCE solver
    double eps_abs = 1.e-8; // absolute error requested
    double eps_rel = 1.e-10; // relative error requested

    // initialize GSL driver function
    gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                   h, eps_abs, eps_rel);

    double tmin    = 0.; // starting t value
    double tmax    = 1.e8; // final t value [year]
    double delta_t = 100.; // step in t [year]

    double t       = tmin; // initialize t
    double y[dim];
    y[0] = 1.0; // initial value psi
    int status; // status of driver function

    for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {
        status = gsl_odeiv2_driver_apply(d, &t, t_next, y);
        if (status != GSL_SUCCESS) 
        {
            printf("Error: status = %d \n", status);
            break;
        }
        printf ("%.5e %.5e \n", t, y[0]); // print psi(t=t_next)
        fprintf(outdata,"%.5e  %.5e\n", t, y[0]);
    }
    gsl_odeiv2_driver_free(d);
    fclose (outdata);

    // return 
    return 0;
}
