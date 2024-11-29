#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "allvars.h"
#include "proto.h"



static double fc(double c)
{
  return c * (0.5 - 0.5 / pow(1 + c, 2) - log(1 + c) / (1 + c)) / pow(log(1 + c) - c / (1 + c), 2);
}

static double jdisk_int(double x, void *param)
{
  double vc2, Sigma0, vc, y;

  if(x > ProfileTable_r[0])
    vc2 = All.G * (halo_get_mass_inside_radius(x) + bulge_get_mass_inside_radius(x)) / x;
  else
    vc2 = 0;

  if(vc2 < 0)
    terminate("vc2 < 0");

  Sigma0 = All.Disk_Mass / (2 * M_PI * All.Disk_R * All.Disk_R);
  y = x / (2 * All.Disk_R);

  if(y > 1e-4)
    vc2 +=
      x * 2 * M_PI * All.G * Sigma0 * y * (gsl_sf_bessel_I0(y) * gsl_sf_bessel_K0(y) -
					   gsl_sf_bessel_I1(y) * gsl_sf_bessel_K1(y));

  vc = sqrt(vc2);

  return pow(x / All.Disk_R, 2) * vc * exp(-x / All.Disk_R);
}


static double gc_int(double x, void *param)
{
  return pow(log(1 + x) - x / (1 + x), 0.5) * pow(x, 1.5) / pow(1 + x, 2);
}




void structure_determination(void)
{
  if (All.HaloUseTable != 0)
  {
    load_profile_table();
    All.PeakDens = DensTable[0];
    All.CoreRadius = pow(0.0019/All.m_22/All.m_22/All.PeakDens,0.25)*3.085678e21/All.UnitLength_in_cm;
    All.CoreEnMass = soliton_enclosed_mass(0.8*All.CoreRadius);
    All.Halo_Mass = halo_get_mass_inside_radius(1.0e+4);
  }
  else
  {
    All.Halo_Mass = 7.0;
  }
  All.M200 = All.Halo_Mass + All.Disk_Mass;
  All.MD = All.Disk_Mass/All.M200;
  All.Bulge_Mass = All.MB * All.M200;
  All.R347 = pow(All.G*All.Halo_Mass/100./All.Hubble/All.Hubble*(347./200.), 1./3.);
  All.R200 = pow(All.G*All.M200/100./All.Hubble/All.Hubble, 1./3.);
  All.V200 = pow(All.G*All.M200/All.R200, 0.5);
  All.LowerDispLimit = pow(All.DispLimitRatio * All.V200, 2);
  All.Halo_C = 8;
  All.Halo_Rs = All.R347 / All.Halo_C;
  All.Halo_A = All.Halo_Rs * sqrt(2 * (log(1 + All.Halo_C) - All.Halo_C / (1 + All.Halo_C)));

  All.BH_Mass = All.MBH * All.M200;
  if(All.MBH > 0)
    All.BH_N = 1;
  else
    All.BH_N = 0;
  All.Disk_Z0 = All.DiskHeight * All.Disk_R;
  mpi_printf("\nStructural parameters:\n");
  if (All.HaloUseTable)
  {
    mpi_printf("PeakDens        = %g\n", All.PeakDens);
    mpi_printf("CoreRadius      = %g\n", All.CoreRadius);
    mpi_printf("CoreEnMass      = %g\n", All.CoreEnMass);
  }
  else
  {
    mpi_printf("A (halo)        = %g\n", All.Halo_A);
  }
  mpi_printf("R200            = %g\n", All.R200);
  mpi_printf("V200            = %g\n", All.V200);
  mpi_printf("M200            = %g  (this is the total mass)\n", All.M200);
  mpi_printf("Halo_Mass       = %g\n", All.Halo_Mass);
  mpi_printf("Disk_Mass       = %g\n", All.Disk_Mass);

  All.Bulge_A = All.BulgeSize * All.Halo_A;	// this will be used if no disk is present

  MType[1] = All.Halo_Mass;
  MType[2] = All.Disk_Mass;
  MType[3] = All.Bulge_Mass;

  NType[1] = All.Halo_N;
  NType[2] = All.Disk_N;
  NType[3] = All.Bulge_N;

  mpi_printf("R  (disk)       = %g\n", All.Disk_R);
  mpi_printf("Z0 (disk)       = %g\n", All.Disk_Z0);
  mpi_printf("MD (disk)       = %g\n", All.MD);
//  mpi_printf("A (bulge)       = %g\n", All.Bulge_A);
}


double structure_disk_angmomentum(void)
{
  gsl_function F;
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &jdisk_int;

  double result, abserr;

  gsl_integration_qag(&F, 0, dmin(30 * All.Disk_R, All.R200),
		      0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  result *= All.Disk_Mass;

  gsl_integration_workspace_free(workspace);

  return result;
}


double structure_gc(double c)
{
  gsl_function F;
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &gc_int;

  double result, abserr;

  gsl_integration_qag(&F, 0, c, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return result;
}

void load_profile_table(void)
{
  int STR_SIZE = 1024;
  int counter, k, i;
  FILE *input_table;
  char buff[STR_SIZE];
  double radius, dens, enmass, pot;
  ProfileTable_r   = mymalloc("ProfileTable_r"  , sizeof(double)*All.Nbin_Profile);
  DensTable     = mymalloc("DensTable"    , sizeof(double)*All.Nbin_Profile);
  EnMassTable   = mymalloc("EnMassTable"  , sizeof(double)*All.Nbin_Profile);
  PoteTable     = mymalloc("PoteTable"    , sizeof(double)*All.Nbin_Profile);

  //Read Table
  counter = 0;
  input_table = fopen(All.ProfileTableFile,"r");
  if (!input_table)
  {
    mpi_printf("File %s cannot be opened! Exit!\n", All.ProfileTableFile);
    exit(1);
  }
  else
    mpi_printf("Loading %s for reading profile table succeed!\n", All.ProfileTableFile);
  while (!feof(input_table))
  {
    fgets(buff, STR_SIZE, input_table);
    if (buff==NULL||buff[0]=='\0')
      break;
    else if (buff[0]!='#')
    {
      sscanf(buff, "%lf %lf %lf %lf", &radius, &dens, &enmass, &pot);
      ProfileTable_r[counter] = radius;
      DensTable[counter] = dens*pow(All.UnitLength_in_cm, 3)/All.UnitMass_in_g;
      EnMassTable[counter] = enmass/1e+10;
      PoteTable[counter] = pot/All.UnitVelocity_in_cm_per_s/All.UnitVelocity_in_cm_per_s;

      counter ++;
    }
    memset(buff, '\0', STR_SIZE);
  }
  fclose(input_table);
  All.r_min = ProfileTable_r[0];
  All.r_max = ProfileTable_r[All.Nbin_Profile-1];
  All.r_ratio = pow((All.r_max/All.r_min), 1.0/(All.Nbin_Profile-1.0));
}
void free_profile_table(void)
{
  myfree(ProfileTable_r);
  myfree(DensTable);
  myfree(EnMassTable);
  myfree(PoteTable);
}
