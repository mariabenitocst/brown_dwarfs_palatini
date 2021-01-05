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
#include <vector>
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
  double b1    = _params[0];
  double nu    = _params[1];
  double mu_1  = _params[2];
  double mu_e  = _params[3];
  double kR    = _params[4];
  double omega = _params[5];
  double gamma = _params[6];
  double delta = _params[7];
  double alpha = _params[8];

  double b     = (-5./16.*y[0]*log(1+exp(-1/y[0])) + 
                  15./8.*y[0]*y[0]*(pow(M_PI, 2)/2. + Li2(-exp(-1/y[0]))));
  double a     = 2.5*mu_e/mu_1;

  double F = 0.0141*pow(b1, 0.266)*pow(kR, -0.106)*pow(y[0], 0.266*nu-0.59)*
                    pow(1+b+a*y[0], 1.506);
  double M = F*pow(gamma, 1.506)*pow(delta, -0.512)*pow(1-1.33*alpha/delta, 0.106);

  double f_1 = 1.1634e-18*pow(b1, 2.856)*mu_1/(pow(kR, 1.1424)*pow(mu_e, 8./3.));
  double f_2 = pow(gamma, 0.7143)*pow(1-1.33*alpha/delta, 1.143)/omega;

  double change_to_year = 3.1536e7; // seconds/year

  f[0] = -f_1*change_to_year*f_2*pow(M, -1.094)*
         pow(y[0], 2.856*nu)*pow(1+b+a*y[0], 1.715);

  // return
  return GSL_SUCCESS;
}


vector<double> degeneracy(double alpha,
                          double gamma,
                          double delta,
                          double tmin,
                          double tmax,
                          double delta_t)
{
    size_t dim = 1;
    double params[9];
    /** Model D   */
    params[0] = 2.; //b1
    params[1] = 1.60; //nu
    /** */
    double _X = 0.75; // mass fraction H
    double _Y = 0.25; // mass fracion He
    // mean molecual weight for He and ionized H mixture
    params[2] = pow((1+0.5*0.51)*_X+_Y/4., -1);
    params[3] = pow(_X + 0.5*_Y, -1); // number of baryons per electron
    params[4] = 0.01; // Rossland opacity [cm2/g]
    params[5] = 1.; // Omega
    params[6] = gamma; // gamma
    params[7] = delta; // delta
    params[8] = alpha; //alpha

    // declare ODE system
    gsl_odeiv2_system sys = {func, NULL, dim, &params};

    double h       = 1.e-6; // starter step size for OCE solver
    double eps_abs = 1.e-8; // absolute error requested
    double eps_rel = 1.e-10; // relative error requested

    // initialize GSL driver function
    gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                   h, eps_abs, eps_rel);

    //double tmin    = 1e6; // starting t value
    //double tmax    = 1e10; // final t value [year]
    //double delta_t = 1e4; // step in t [year]

    double t       = tmin; // initialize t
    double y[dim];
    y[0] = 1.0; // initial value psi
    int status; // status of driver function

    vector<double> psi; // degeneracy

    for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {
        status = gsl_odeiv2_driver_apply(d, &t, t_next, y);
        if (status != GSL_SUCCESS) 
        {
            printf("Error: status = %d \n", status);
            break;
        }
        psi.push_back(y[0]);
    }
    gsl_odeiv2_driver_free(d);

    // return
    return psi;
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

    double tmin    = 1e7; // starting t value
    double tmax    = 1e10; // ending t value
    double delta_t = 1e4; // step in t [year]

    vector<double> psi;

    psi = degeneracy(alpha, gamma, delta, tmin, tmax, delta_t);

    // output file
    stringstream filename;
    filename << "../../data/evolution_degeneracy_MMSM_alpha=" << alpha << ".dat";
    FILE * outdata = fopen (filename.str().c_str(),"w"); 
    fprintf(outdata,"# time [year]  degeneracy");

    int j = 0;
    for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {
        fprintf(outdata,"%.5e %.5e\n", t_next, psi[j]);
        j += 1;
    }
    fclose (outdata);

    // return 
    return 0;
}
