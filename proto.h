#ifndef PROTO_H
#define PROTO_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_linalg.h>

extern class GlobalVariables
{
  public:
    double X; // Hydrogen mass-fraction
    double Y; // Helium mass-fraction
    double op_ratio; // Ratio of ortho- to parahydrogen
    
    double k = 1.38065e-16;
    double h = 6.6260695e-27;
    
    double theta_rot = 170.64;
    double theta_vib = 5984.48;
    double chi_diss = 7.17e-12;
    double chi_ion = 2.18e-11;
    double chi_He1 = 3.94e-11; 
    double chi_He2 = 8.72e-11;
    
    double mP = 1.67262178e-24;
    double me = 9.1093829e-28;
    double mu_H = 1.007276466621;
    double mu_He = 4.002602;
} All;

double K_dis(double T);
double K_ion(double T);
double K_He1(double T);
double K_He2(double T);
void calculate_abundances(double rho, double T, double abundances[7]);
double solve_for_electron_abundance(double nH_tot, double nHe_tot, double K_dis, double K_ion, double K_He1, double K_He2);
double f_ne(double ne, double nH_tot, double nHe_tot, double K_dis, double K_ion, double K_He1, double K_He2);

void dln_n_dln_rho(double rho, double abundances[7], double dln_abundances_dln_rho[7]);
void dln_n_dln_T(double T, double abundances[7], double dln_abundances_dln_T[7]);
void solve_matrix_equation(double *a_data, double *b_data, double x_data[7]);

double get_partition_derivative(int i, double T);
double dln_z_H2(double T);
double dln_z_H2_rot(double T);
double dln_z_H2_vib(double T);
double dln_z_H(double T);
double dln_z_H_elec(double T);
double dln_z_Hp(double T);
double dln_z_Hp_elec(double T);
double dln_z_He(double T);
double dln_z_Hep(double T);
double dln_z_Hep_elec(double T);
double dln_z_Hep2(double T);
double dln_z_Hep2_elec(double T);
double dln_z_e(double T);

double get_partition_function(int i, double T);
double z_H2(double T);
double z_H2_tr(double T);
double z_H2_rot(double T);
double z_H2_vib(double T);
double z_H(double T);
double z_H_tr(double T);
double z_H_elec(double T);
double z_Hp(double T);
double z_Hp_tr(double T);
double z_Hp_elec(double T);
double z_He(double T);
double z_He_tr(double T);
double z_Hep(double T);
double z_Hep_tr(double T);
double z_Hep_elec(double T);
double z_Hep2(double T);
double z_Hep2_tr(double T);
double z_Hep2_elec(double T);
double z_e(double T);
double z_e_tr(double T);

double get_partition_second_derivative(int i, double T);

double d2ln_z_H2(double T);
double d2ln_z_H2_rot(double T);
double d2ln_z_H2_vib(double T);
double d2ln_z_H(double T);
double d2ln_z_H_elec(double T);
double d2ln_z_Hp(double T);
double d2ln_z_Hp_elec(double T);
double d2ln_z_He(double T);
double d2ln_z_Hep(double T);
double d2ln_z_Hep_elec(double T);
double d2ln_z_Hep2(double T);
double d2ln_z_Hep2_elec(double T);
double d2ln_z_e(double T);

#endif
