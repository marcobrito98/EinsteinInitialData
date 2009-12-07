#include <cstdio>
#include <cassert>
#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <bin_ns.h>

using namespace std;



extern "C"
void ID_Bin_NS_initialise (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INFO ("Setting up LORENE Bin_NS initial data");
  
  // Meudon data are distributed in SI units (MKSA).  Here are some
  // conversion factors.
  
  // Defined constants
  CCTK_REAL const c_light = 299792458.0; // speed of light [m/s]
  
  // Constants of nature (IAU, CODATA):
  CCTK_REAL const G_grav = 6.67428e-11; // gravitational constant [m^3/kg/s^2]
  CCTK_REAL const M_sun  = 1.98892e+30; // solar mass [kg]
  
  // Cactus units in terms of SI units:
  // (These are derived from M = M_sun, c = G = 1, and using 1/M_sun
  // for the magnetic field)
  CCTK_REAL const cactusM = M_sun;
  CCTK_REAL const cactusL = cactusM * G_grav / pow(c_light,2);
  CCTK_REAL const cactusT = cactusL / c_light;
  
  // Other quantities in terms of Cactus units
  CCTK_REAL const coord_unit = cactusL / 1.0e+3;         // from km
  CCTK_REAL const rho_unit   = cactusM / pow(cactusL,3); // from kg/m^3
  CCTK_REAL const ener_unit  = pow(cactusL,2);           // from c^2
  CCTK_REAL const vel_unit   = cactusL / cactusT / c_light; // from c
  
  
  
  CCTK_INFO ("Setting up coordinates");
  
  int const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];
  
  vector<double> xx(npoints), yy(npoints), zz(npoints);
  
#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    xx[i] = x[i] * coord_unit;
    yy[i] = y[i] * coord_unit;
    zz[i] = z[i] * coord_unit;
  }
  
  
  
  CCTK_VInfo (CCTK_THORNSTRING, "Reading from file \"%s\"", filename);
  
  Bin_NS bin_ns (npoints, &xx[0], &yy[0], &zz[0], filename);
  
  CCTK_VInfo (CCTK_THORNSTRING, "omega [rad/s]:       %g", bin_ns.omega);
  CCTK_VInfo (CCTK_THORNSTRING, "dist [km]:           %g", bin_ns.dist);
  CCTK_VInfo (CCTK_THORNSTRING, "dist_mass [km]:      %g", bin_ns.dist_mass);
  CCTK_VInfo (CCTK_THORNSTRING, "mass1_b [M_sun]:     %g", bin_ns.mass1_b);
  CCTK_VInfo (CCTK_THORNSTRING, "mass2_b [M_sun]:     %g", bin_ns.mass2_b);
  CCTK_VInfo (CCTK_THORNSTRING, "mass_ADM [M_sun]:    %g", bin_ns.mass_adm);
  CCTK_VInfo (CCTK_THORNSTRING, "L_tot [G M_sun^2/c]: %g", bin_ns.angu_mom);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_comp [km]:    %g", bin_ns.rad1_x_comp);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_y [km]:         %g", bin_ns.rad1_y);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_z [km]:         %g", bin_ns.rad1_z);
  CCTK_VInfo (CCTK_THORNSTRING, "rad1_x_opp [km]:     %g", bin_ns.rad1_x_opp);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_comp [km]:    %g", bin_ns.rad2_x_comp);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_y [km]:         %g", bin_ns.rad2_y);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_z [km]:         %g", bin_ns.rad2_z);
  CCTK_VInfo (CCTK_THORNSTRING, "rad2_x_opp [km]:     %g", bin_ns.rad2_x_opp);
  assert (bin_ns.np == npoints);
  
  
  
  CCTK_INFO ("Filling in Cactus grid points");
  
#pragma omp parallel for
  for (int i=0; i<npoints; ++i) {
    
    alp[i] = bin_ns.nnn[i];
    
    betax[i] = bin_ns.beta_x[i];
    betay[i] = bin_ns.beta_y[i];
    betaz[i] = bin_ns.beta_z[i];
    
    // These initial data assume a helical Killing vector field
    // TODO: calculate spatial derivatives to set this correctly
    dtalp[i] = 0.0;
    
    dtbetax[i] = 0.0;
    dtbetay[i] = 0.0;
    dtbetaz[i] = 0.0;
    
    CCTK_REAL g[3][3];
    g[0][0] = bin_ns.g_xx[i];
    g[0][1] = bin_ns.g_xy[i];
    g[0][2] = bin_ns.g_xz[i];
    g[1][1] = bin_ns.g_yy[i];
    g[1][2] = bin_ns.g_yz[i];
    g[2][2] = bin_ns.g_zz[i];
    g[1][0] = g[0][1];
    g[2][0] = g[0][2];
    g[2][1] = g[1][2];
    
    CCTK_REAL ku[3][3];
    ku[0][0] = bin_ns.k_xx[i];
    ku[0][1] = bin_ns.k_xy[i];
    ku[0][2] = bin_ns.k_xz[i];
    ku[1][1] = bin_ns.k_yy[i];
    ku[1][2] = bin_ns.k_yz[i];
    ku[2][2] = bin_ns.k_zz[i];
    ku[1][0] = ku[0][1];
    ku[2][0] = ku[0][2];
    ku[2][1] = ku[1][2];
    
    CCTK_REAL k[3][3];
    for (int a=0; a<3; ++a) {
      for (int b=0; b<3; ++b) {
        k[a][b] = 0.0;
        for (int c=0; c<3; ++c) {
          for (int d=0; d<3; ++d) {
            k[a][b] += g[a][c] * g[b][d] * ku[c][d];
          }
        }
      }
    }
    
    gxx[i] = g[0][0];
    gxy[i] = g[0][1];
    gxz[i] = g[0][2];
    gyy[i] = g[1][1];
    gyz[i] = g[1][2];
    gzz[i] = g[2][2];
    
    kxx[i] = k[0][0];
    kxy[i] = k[0][1];
    kxz[i] = k[0][2];
    kyy[i] = k[1][1];
    kyz[i] = k[1][2];
    kzz[i] = k[2][2];
    
    rho[i] = bin_ns.nbar[i] / rho_unit;
    
    eps[i] = rho[i] * bin_ns.ener_spec[i] / ener_unit;
    
    vel[i          ] = bin_ns.u_euler_x[i] / vel_unit;
    vel[i+  npoints] = bin_ns.u_euler_y[i] / vel_unit;
    vel[i+2*npoints] = bin_ns.u_euler_z[i] / vel_unit;
    
  } // for i
  
  
  
  CCTK_INFO ("Done.");
}
