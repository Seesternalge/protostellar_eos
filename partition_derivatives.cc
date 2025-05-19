#include "proto.h"

double get_partition_derivative(int i, double T)
{
  if(i == 0) return dln_z_H2(T);
  else if(i == 1) return dln_z_H(T);
  else if(i == 2) return dln_z_Hp(T);
  else if(i == 3) return dln_z_He(T);
  else if(i == 4) return dln_z_Hep(T);
  else if(i == 5) return dln_z_Hep2(T);
  else if(i == 6) return dln_z_e(T);
  else return 0.0;
}

// Molecular hydrogen

double dln_z_H2(double T)
{
  return 3.0 / 2.0 + dln_z_H2_rot(T) + dln_z_H2_vib(T);
}

double dln_z_H2_rot(double T)
{
  double theta_rot = All.theta_rot;
  double z_even = 0.0, z_odd = 0.0;
  
  for(int j = 0 ; j < 10000 ; j++)
  {
    if(j % 2 == 0) z_even += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T);
    else           z_odd  += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T);
  }
  
  double dln_z_even = 1.0 / z_even, dln_z_odd = 1.0 / z_odd;
  double numerator_even = 0.0, numerator_odd = 0.0;
  
  for(int j = 0 ; j < 10000 ; j++)
  {
    if(j % 2 == 0) numerator_even += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T) * j * (j + 1.0) / 2.0;
    else           numerator_odd  += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T) * j * (j + 1.0) / 2.0;
  }  
  
  dln_z_even *= numerator_even;
  dln_z_odd *= numerator_odd;
  
  double result = (1.0 / (All.op_ratio + 1.0) * dln_z_even + All.op_ratio / (All.op_ratio + 1.0) * dln_z_odd - All.op_ratio / (All.op_ratio + 1.0)) * theta_rot / T / T;
  result *= T; 
  
  return result;
}

double dln_z_H2_vib(double T) /* Appears to be a mistake in Tomida 2013; this is the standard vibrational partition function*/
{
  double theta_vib = All.theta_vib;
  double result = theta_vib / T / T / (exp(theta_vib / T) - 1.0);
  result *= T;
  return result;
}

// Atomic hydrogen

double dln_z_H(double T)
{
  return 3.0 / 2.0 + T * dln_z_H_elec(T);
}

double dln_z_H_elec(double T)
{
  double k = All.k;
  double chi_diss = All.chi_diss;
  return chi_diss / 2.0 / k / T / T;
}

// Ionized hydrogen

double dln_z_Hp(double T)
{
  return 3.0 / 2.0 + T * dln_z_Hp_elec(T);
}

double dln_z_Hp_elec(double T)
{
  double k = All.k;
  double chi_diss = All.chi_diss, chi_ion = All.chi_ion;
  return (chi_diss + 2.0 * chi_ion) / 2.0 / k / T / T;
}

// Atomic helium

double dln_z_He(double T)
{
  return 3.0 / 2.0;
}

// Ionized helium

double dln_z_Hep(double T)
{
  return 3.0 / 2.0 + T * dln_z_Hep_elec(T);
}

double dln_z_Hep_elec(double T)
{
  double k = All.k;
  double chi_He1 = All.chi_He1;
  return chi_He1 / k / T / T;
}

// Doubly ionized helium

double dln_z_Hep2(double T)
{
  return 3.0 / 2.0 + T * dln_z_Hep2_elec(T);
}

double dln_z_Hep2_elec(double T)
{
  double k = All.k;
  double chi_He1 = All.chi_He1, chi_He2 = All.chi_He2;
  return (chi_He1 + chi_He2) / k / T / T;
}

// Electron

double dln_z_e(double T)
{
  return 3.0 / 2.0;
}

