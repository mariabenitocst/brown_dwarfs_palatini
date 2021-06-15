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
  double M     = _params[0];
  double b1    = _params[1];
  double nu    = _params[2];
  double mu_1  = _params[3];
  double mu_e  = _params[4];
  double kR    = _params[5];
  double omega = _params[6];
  double gamma = _params[7];
  double delta = _params[8];
  double alpha = _params[9];

  double b     = (-5./16.*y[0]*log(1+exp(-1/y[0])) + 
                  15./8.*y[0]*y[0]*(pow(M_PI, 2)/2. + Li2(-exp(-1/y[0]))));
  double a     = 2.5*mu_e/mu_1;


  double f_1 = 1.1634e-18*pow(b1, 2.856)*mu_1/(pow(kR, 1.1424)*pow(mu_e, 8./3.));
  double f_2 = pow(gamma, 0.7143)*pow(1-1.33*alpha/delta, 1.143)/omega;

  double change_to_year = 3.1536e7; // seconds/year

  f[0] = -f_1*change_to_year*f_2*pow(M, -1.094)*
         pow(y[0], 2.856*nu)*pow(1+b+a*y[0], 1.715);

  // return
  return GSL_SUCCESS;
}


vector<double> degeneracy(double mass,
                          double alpha,
                          double gamma,
                          double delta,
                          double tmin,
                          double tmax,
                          double delta_t)
{
    size_t dim = 1;
    double params[11];
    params[0] = mass; // mass [Msun]
    /** Model D   */
    params[1] = 2.; //b1
    params[2] = 1.60; //nu
    /** */
    double _X = 0.75; // mass fraction H
    double _Y = 0.25; // mass fracion He
    // mean molecual weight for He and ionized H mixture
    params[3] = pow((1+0.5*0.51)*_X+_Y/4., -1);
    params[4] = pow(_X + 0.5*_Y, -1); // number of baryons per electron
    params[5] = 0.01; // Rossland opacity [cm2/g]
    params[6] = 1.; // Omega
    params[7] = gamma; // gamma
    params[8] = delta; // delta
    params[9] = alpha; //alpha

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

    // output file
    stringstream filename;
    filename << "../../data/evolution_degeneracy_alpha=" << alpha << ".dat";
    FILE * outdata = fopen (filename.str().c_str(),"w"); 
    fprintf(outdata,"# time [year]  degeneracy\n");

    double mass[5]; // mass [Msun]
    mass[0] = 0.01;
    mass[1] = 0.03;
    mass[2] = 0.05;
    mass[3] = 0.08;
    mass[4] = 0.09;

    double tmin    = 1e6; // starting t value
    double tmax    = 1e10; // final t value [year]
    double delta_t = 1e4; // step in t [year]

    vector<double> psi[5];

    for (int i=0; i<=4; i++)
    {
        psi[i] = degeneracy(mass[i], alpha, gamma, delta, tmin, tmax, delta_t);
    }

    int j = 0;
    for (double t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {
        fprintf(outdata,"%.5e %.5e %.5e %.5e %.5e %.5e\n", t_next, psi[0][j],
                        psi[1][j], psi[2][j], psi[3][j], psi[4][j]);
        j += 1;
    }
    fclose (outdata);

    // return 
    return 0;
}
