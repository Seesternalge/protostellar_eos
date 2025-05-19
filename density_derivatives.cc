#include <cmath>
#include <gsl/gsl_linalg.h>

using namespace std;

void dln_n_dln_rho(double X, double Y, double rho, double abundances[7], double dln_abundances_dln_rho[7]);
void dln_n_dln_T(double T, double abundances[7], double dln_abundances_dln_T[7]);

void dln_n_dln_rho(double X, double Y, double rho, double abundances[7], double dln_abundances_dln_rho[7])
{
  double nH_tot = rho * X / (1.67262178e-24 + 9.1093829e-28), nHe_tot = rho * Y / (4.002602 / 1.007276466621 * 1.67262178e-24);
  
  double a_data[] = {1.0 , - 2.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, // Some typos in Tomida+2013; they seem to have mixed up the positions of e- and He in the matrix
                     2.0 * abundances[0] , abundances[1] , abundances[2] , 0.0 , 0.0 , 0.0 , 0.0,
                     0.0 , - 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0, 
                     0.0 , 0.0 , abundances[2] , 0.0, abundances[4] , 2.0 * abundances[5] , - abundances[6],
                     0.0 , 0.0 , 0.0 , - 1.0 , 1.0 , 0.0 , 1.0,
                     0.0 , 0.0 , 0.0 , 0.0 , - 1.0 , 1.0 , 1.0,
                     0.0 , 0.0 , 0.0 , abundances[3] , abundances[4] , abundances[5] , 0.0};
                     
  double b_data[] = {0.0 , nH_tot , 0.0 , 0.0 , 0.0 , 0.0 , nHe_tot};
  
  
  gsl_matrix_view m = gsl_matrix_view_array(a_data, 7, 7);
  gsl_vector_view b = gsl_vector_view_array(b_data, 7);
  gsl_vector *x = gsl_vector_alloc(7);
  int s;
  gsl_permutation *p = gsl_permutation_alloc(7);
  gsl_linalg_LU_decomp(&m.matrix, p, &s);
  gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);

  for(int k = 0 ; k < 7 ; k++) dln_abundances_dln_rho[k] = gsl_vector_get(x,k);

  gsl_permutation_free(p);
  gsl_vector_free(x);
}

void dln_n_dln_T(double T, double abundances[7], double dln_abundances_dln_T[7])
{ 
  double a_data[] = {1.0 , - 2.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0,
                     2.0 * abundances[0] , abundances[1] , abundances[2] , 0.0 , 0.0 , 0.0 , 0.0,
                     0.0 , - 1.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0, 
                     0.0 , 0.0 , abundances[2] , 0.0, abundances[4] , 2.0 * abundances[5] , - abundances[6],
                     0.0 , 0.0 , 0.0 , - 1.0 , 1.0 , 0.0 , 1.0,
                     0.0 , 0.0 , 0.0 , 0.0 , - 1.0 , 1.0 , 1.0,
                     0.0 , 0.0 , 0.0 , abundances[3] , abundances[4] , abundances[5] , 0.0};
                     
  double b_data[] = {dln_z_H2(T) - 2.0 * dln_z_H(T) , 0.0 , dln_z_Hp(T) + dln_z_e(T) - dln_z_H(T) , 0.0 , dln_z_Hep(T) + dln_z_e(T) - dln_z_He(T) , dln_z_Hep2(T) + dln_z_e(T) - dln_z_Hep(T) , 0.0};
  
  
  gsl_matrix_view m = gsl_matrix_view_array(a_data, 7, 7);
  gsl_vector_view b = gsl_vector_view_array(b_data, 7);
  gsl_vector *x = gsl_vector_alloc(7);
  int s;
  gsl_permutation *p = gsl_permutation_alloc(7);
  gsl_linalg_LU_decomp(&m.matrix, p, &s);
  gsl_linalg_LU_solve(&m.matrix, p, &b.vector, x);

  for(int k = 0 ; k < 7 ; k++) dln_abundances_dln_T[k] = gsl_vector_get(x,k);

  gsl_permutation_free(p);
  gsl_vector_free(x);
}
