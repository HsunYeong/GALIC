#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include "allvars.h"
#include "proto.h"


/* this file contains auxiliary routines for the description of the halo,
 * here modeled as a Hernquist sphere
 */

double interpolation_M2R(double m)
{
  int i = 0, j = All.Nbin_Profile - 1;
  while(i <= j)
  {
    int mid = i + (j - i)/2;
    if(EnMassTable[mid] == m)
    {
       i = mid;
       break;
    }
    else if(EnMassTable[mid] > m) j = mid - 1;
    else i = mid + 1;
  }
  return exp(log(ProfileTable_r[i-1])+(log(ProfileTable_r[i])-log(ProfileTable_r[i-1]))*(log(m)-log(EnMassTable[i-1]))/(log(EnMassTable[i])-log(EnMassTable[i-1])));
}
double get_value_from_interpolation_log(double r, double *array_value)
{
  int i = (int)(log(r/All.r_min)/log(All.r_ratio));
  while(r > ProfileTable_r[i])
  {
    i = i + 1;
  }
  return exp(log(array_value[i-1])+(log(array_value[i])-log(array_value[i-1]))*(log(r)-log(ProfileTable_r[i-1]))/(log(ProfileTable_r[i])-log(ProfileTable_r[i-1])));
}
double get_value_from_interpolation(double r, double *array_value)
{
  int i = (int)(log(r/All.r_min)/log(All.r_ratio));
  while(r > ProfileTable_r[i])
  {
    i = i + 1;
  }
  return array_value[i-1]+(array_value[i]-array_value[i-1])*(log(r)-log(ProfileTable_r[i-1]))/(log(ProfileTable_r[i])-log(ProfileTable_r[i-1]));
}

/* this function returns a new random coordinate for the halo */
void halo_get_fresh_coordinate(double *pos)
{
  double r;
  double targetM, rnew, dr, M;
  int count = 0;
  do
    {
      double q = gsl_rng_uniform(random_generator);
      if (All.HaloUseTable == 0)
      {

        if(q > 0)
          r = All.Halo_A * (q + sqrt(q)) / (1 - q);
        else
          r = 0;
      }
      else
      {
         r = halo_mass_to_radius(q*All.Halo_Mass);
      }
      double phi = gsl_rng_uniform(random_generator) * M_PI * 2;
      double theta = acos(gsl_rng_uniform(random_generator) * 2 - 1);

      pos[0] = r * sin(theta) * cos(phi);
      pos[1] = r * sin(theta) * sin(phi);
      pos[2] = r * cos(theta) / All.HaloStretch;

      r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
      count += 1;
      if (count > 1e3)
      {
         mpi_printf("Too many loops! r = %g\n", r);
         count = 0;
      }
    }
  while(r > All.Rmax);
}

double halo_mass_to_radius(double m)
{
  double r, M, rnew, dr;
  if (m >= EnMassTable[All.Nbin_Profile-1]) return ProfileTable_r[All.Nbin_Profile-1];
  else if(m <= EnMassTable[0])
  {
    if(m > 0)
    {
      r = 0.5*ProfileTable_r[0];
      do
      {
        M = soliton_enclosed_mass(r);
        rnew = m / M * r;
        dr = rnew - r;
        if (fabs(dr) > 0.5*r) dr = 0.5* r * dr / fabs(dr);
        else dr = dr * 0.1;
        r = r + dr;
      }
      while(fabs(dr)/r > 1.0e-4);
    }
    else r = 0;
    return r;
  }
  else return interpolation_M2R(m);
}
double halo_get_density(double *pos) {

  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  if(All.HaloUseTable == 0)
  {
     return All.Halo_Mass / (2 * M_PI) * All.Halo_A / (r + 1.0e-6 * All.Halo_A) / pow(r + All.Halo_A, 3);

  }
  else
  {
     return halo_get_density_from_radius(r);
  }
}

double halo_get_density_from_radius(double r) {
  if (r < ProfileTable_r[0]) return DensTable[0];
  else if (r > ProfileTable_r[All.Nbin_Profile-1]) return 0;
  else return get_value_from_interpolation_log(r, DensTable);
}

double soliton_density(double r) {

  double rho =  0.0019 / All.m_22 / All.m_22 * pow(1.0 / All.CoreRadius/ pow(1 + 0.091 * pow(r / All.CoreRadius, 2), 2), 4);

  if ( fabs(rho) <  MIN_DENSITY) rho = 0;

  return rho;
}

double soliton_enclosed_mass(double r)
{
  if (r==0) return 0;
  else
  {
    double r_sol = All.CoreRadius;
    double A = 0.0019/All.m_22/All.m_22/pow(r_sol, 4);
    double B = 0.091;
    double soliton_mass = M_PI*pow(r_sol,3)*A/53760/pow(B,1.5)*(r_sol*sqrt(B)*r/pow(r_sol*r_sol+B*r*r,7)
           *(-3465*pow(r_sol,12)+48580*pow(r_sol,10)*B*r*r+92323*pow(r_sol,8)*B*B*pow(r,4)+101376*pow(r_sol,6)*pow(B,3)*pow(r,6)
           +65373*pow(r_sol,4)*pow(B,4)*pow(r,8)+23100*r_sol*r_sol*pow(B,5)*pow(r,10)+3465*pow(B,6)*pow(r,12))
           +3465*atan(sqrt(B)*r/r_sol));
    return soliton_mass;
  }
}
double halo_get_mass_inside_radius(double r)
{
  if (All.HaloUseTable == 0)
  {
    return All.Halo_Mass * pow(r / (r + All.Halo_A), 2);
  }
  else
  {
    if (r < ProfileTable_r[0]) return soliton_enclosed_mass(r);
    else if( r > ProfileTable_r[All.Nbin_Profile-1]) return EnMassTable[All.Nbin_Profile-1];
    else return get_value_from_interpolation_log(r, EnMassTable);
  }
}
double hyperg_z_GT1(double a, double b, double c, double z)
{
  // calculate 2F1 for z < -1
  double coef1, coef2;

  coef1 = gsl_sf_gamma(c)*gsl_sf_gamma(b-a)*pow(1-z, -a)/(gsl_sf_gamma(b)*gsl_sf_gamma(c-a));
  coef2 = gsl_sf_gamma(c)*gsl_sf_gamma(a-b)*pow(1-z, -b)/(gsl_sf_gamma(a)*gsl_sf_gamma(c-b));

  return coef1*gsl_sf_hyperg_2F1(a, c-b, a-b+1, 1/(1-z)) + coef2*gsl_sf_hyperg_2F1(b, c-a, b-a+1, 1/(1-z));
}
static double potential_int(double r, void *param) {
   return All.G*EnMassTable[All.Nbin_Profile-1]/r/r;
}

double halo_get_potential(double *pos)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  return halo_get_potential_from_radius(r);
}

double halo_get_potential_from_radius(double r)
{
  if (All.HaloUseTable == 0)
  {
    double phi = -All.G * All.Halo_Mass / (r + All.Halo_A);
    return phi;
  }
  else
  {
    if (r < ProfileTable_r[0]) return PoteTable[0];
    else if (r > ProfileTable_r[All.Nbin_Profile-1])
    {
      gsl_function F;
      gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
      F.function = &potential_int;
      double result, abserr;
      gsl_integration_qagiu(&F, r, 0, 1.0e-8, WORKSIZE, workspace, &result, &abserr);
      gsl_integration_workspace_free(workspace);
      return -1*result;
    }
    else return get_value_from_interpolation(r, PoteTable);
  }
}

/* returns the acceleration at coordinate pos[] */
void halo_get_acceleration(double *pos, double *acc)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  double fac;
  if (All.HaloUseTable == 0)
  {
     fac = All.G * All.Halo_Mass / ((r + 1.0e-6 * All.Halo_A)* (r + All.Halo_A) * (r + All.Halo_A));
  }
  else
  {
    if (r == 0) fac = 0;
    else fac = All.G * halo_get_mass_inside_radius(r)/pow(r, 3);
  }
  acc[0] = -fac * pos[0];
  acc[1] = -fac * pos[1];
  acc[2] = -fac * pos[2];
}

static double energy_int(double r, void *param) {
   return 2*M_PI*r*All.G*halo_get_mass_inside_radius(r)*(halo_get_density_from_radius(r));
}
double halo_get_energy(void)
{
  gsl_function F;
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &energy_int;
  double result, abserr;
  gsl_integration_qagiu(&F, 0, 0, 1.0e-8, WORKSIZE, workspace, &result, &abserr);
  gsl_integration_workspace_free(workspace);

  return result;
}
/*double halo_get_escape_speed(double *pos)
{
  double r = sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
  double phi = -All.G * All.Halo_Mass / (r + All.Halo_A);
  double vesc = sqrt(-2.0 * phi);

  return vesc;
}*/

/*double halo_get_sigma2(double *pos) {

	long double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);

	long double m = All.Halo_Mass;
	long double r0 = All.Halo_A;
	long double r_over_r0 = r/r0;

	long double _sigma2 =
	(long double)(All.G*m)/(12.0*r0)*
	fabs( 12*r*powl(r+r0,3)/powl(r0,4)*logl((r+r0)/r)
			-
			r/(r+r0)*(25 + r_over_r0*(52 + 42*r_over_r0 + 12*(r_over_r0*r_over_r0) ) )
		 );

	// precicion big rip so let it be like this for a while
	if (65000<r)
		_sigma2 = pow(halo_get_escape_speed(pos)/3.14, 2);
		//_sigma2 = SQR(vesc(_r)/3.14);


	return _sigma2;

}*/




/*E to q conversion*/
double halo_E_to_q(double E)
{
  return sqrt(-E * All.Halo_A / (All.G * All.Halo_Mass));
}



/*Hernquist density of states (as a function of q)*/
double halo_g_q(double q)
{
  double pre =
    2 * sqrt(2) * M_PI * M_PI * All.Halo_A * All.Halo_A * All.Halo_A * sqrt(All.G * All.Halo_Mass /
									    All.Halo_A);

  return pre * (3 * (8 * q * q * q * q - 4 * q * q + 1) * acos(q) -
		q * sqrt(1 - q * q) * (4 * q * q - 1) * (2 * q * q + 3)) / (3 * q * q * q * q * q);
}


/*Hernquist distribution function (as a function of q)*/
double halo_f_q(double q)
{
  double pre =
    (All.Halo_Mass / (All.Halo_A * All.Halo_A * All.Halo_A)) / (4 * M_PI * M_PI * M_PI *
								pow(2 * All.G * All.Halo_Mass / All.Halo_A,
								    1.5));

  return pre * (3 * asin(q) +
		q * sqrt(1 - q * q) * (1 - 2 * q * q) * (8 * q * q * q * q - 8 * q * q - 3)) / pow(1 - q * q,
												   2.5);
}


/*Hernquist distribution function (as a function of radius and velocity)*/
double halo_f(double rad, double vel)
{
  double E = 0.5 * vel * vel + halo_get_potential_from_radius(rad);
  double q = halo_E_to_q(E);

  return halo_f_q(q);
}


/*generate velocities for Hernquist distribution function with von Neumann rejection technique*/
double halo_generate_v(double rad)
{
  double pot = halo_get_potential_from_radius(rad);
  double v_max = sqrt(-2 * pot);	// escape velocity
  double v_guess, x_aux;
  double f_max = v_max * v_max * halo_f(rad, 0);

  v_guess = gsl_rng_uniform(random_generator) * v_max;
  x_aux = gsl_rng_uniform(random_generator) * f_max;

  while(x_aux > v_guess * v_guess * halo_f(rad, v_guess))
    {
      v_guess = gsl_rng_uniform(random_generator) * v_max;
      x_aux = gsl_rng_uniform(random_generator) * f_max;
    }
  return v_guess;
}
