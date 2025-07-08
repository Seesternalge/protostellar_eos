#include "proto.h"

double get_partition_second_derivative(int i, double T)
{
  if(i == 0) return d2ln_z_H2(T) + dln_z_H2(T);
  else if(i == 1) return d2ln_z_H(T) + dln_z_H(T);
  else if(i == 2) return d2ln_z_Hp(T) + dln_z_Hp(T);
  else if(i == 3) return d2ln_z_He(T) + dln_z_He(T);
  else if(i == 4) return d2ln_z_Hep(T) + dln_z_Hep(T);
  else if(i == 5) return d2ln_z_Hep2(T) + dln_z_Hep2(T);
  else if(i == 6) return d2ln_z_e(T) + dln_z_e(T);
  else return 0.0;
}

/* These now only include only the part of the second logarithmic derivative that has the second linear derivative! */

// Molecular hydrogen

double d2ln_z_H2(double T)
{
  return - 3.0 / 2.0 + T * T * (d2ln_z_H2_rot(T) + d2ln_z_H2_vib(T));
}

double d2ln_z_H2_rot(double T) 
{
  double theta_rot = All.theta_rot;
  double z_even = 0.0, z_odd = 0.0;
  
  for(int j = 0 ; j < 10000 ; j++)
  {
    if(j % 2 == 0) z_even += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T);
    else           z_odd  += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T);
  }
  
  double dln_z_even = 1.0 / z_even / z_even, dln_z_odd = 1.0 / z_odd / z_odd;
  double numerator_even = 0.0, numerator_odd = 0.0;
  
  for(int j = 0 ; j < 10000 ; j++)
  {
    if(j % 2 == 0) numerator_even += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T) * j * (j + 1.0) / 2.0 * theta_rot / T / T;
    else           numerator_odd  += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T) * j * (j + 1.0) / 2.0 * theta_rot / T / T;
  }  
  
  dln_z_even *= numerator_even * numerator_even;
  dln_z_odd *= numerator_odd * numerator_odd;
  
  double result = 0.0;
  
  result += - (1.0 / (All.op_ratio + 1.0) * dln_z_even + All.op_ratio / (All.op_ratio + 1.0) * dln_z_odd - All.op_ratio / (All.op_ratio + 1.0) * 2.0 * theta_rot / pow(T,3.0));
  
  double d2ln_z_even = 1.0 / z_even, d2ln_z_odd = 1.0 / z_odd;
  double new_numerator_even = 0.0, new_numerator_odd = 0.0;
  
  for(int j = 0 ; j < 10000 ; j++)
  {
    if(j % 2 == 0) new_numerator_even += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T) * (pow(j * (j + 1.0) / 2.0 * theta_rot / T / T , 2.0) - (j * (j + 1.0) * theta_rot / pow(T,3.0)));
    else           new_numerator_odd  += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T) * (pow(j * (j + 1.0) / 2.0 * theta_rot / T / T , 2.0) - (j * (j + 1.0) * theta_rot / pow(T,3.0)));
  }  
  
  d2ln_z_even *= new_numerator_even;
  d2ln_z_odd *= new_numerator_odd;
  
  result += 1.0 / (All.op_ratio + 1.0) * d2ln_z_even + All.op_ratio / (All.op_ratio + 1.0) * d2ln_z_odd;
  
  return result;
}

double d2ln_z_H2_vib(double T)
{
  double theta_vib = All.theta_vib;
  if(T < theta_vib / 10.0) T = theta_vib / 10.0; // To prevent overflow
  //return theta_vib * (theta_vib * exp(theta_vib / T) - 2.0 * T * (exp(theta_vib / T) - 1.0)) / (pow(T,4.0) * pow(exp(theta_vib / T) - 1.0 , 2.0));
  return pow(theta_vib , 2.0) / (pow(T , 4.0) * pow(exp(theta_vib / T) - 1.0 , 2.0)) + (pow(theta_vib , 2.0) - 2.0 * theta_vib * T)/ (pow(T , 4.0) * (exp(theta_vib / T) - 1.0));
}

// Atomic hydrogen

double d2ln_z_H(double T)
{
  return - 3.0 / 2.0 + T * T * d2ln_z_H_elec(T);
}

double d2ln_z_H_elec(double T)
{
  double k = All.k;
  double chi_diss = All.chi_diss;
  return - chi_diss / k / T / T / T;
}

// Ionized hydrogen

double d2ln_z_Hp(double T)
{
  return - 3.0 / 2.0 + T * T * d2ln_z_Hp_elec(T);
}

double d2ln_z_Hp_elec(double T)
{
  double k = All.k;
  double chi_diss = All.chi_diss, chi_ion = All.chi_ion;
  return - (chi_diss + 2.0 * chi_ion) / k / T / T / T;
}

// Atomic helium

double d2ln_z_He(double T)
{
  return - 3.0 / 2.0;
}

// Ionized helium

double d2ln_z_Hep(double T)
{
  return - 3.0 / 2.0 + T * T * d2ln_z_Hep_elec(T);
}

double d2ln_z_Hep_elec(double T)
{
  double k = All.k;
  double chi_He1 = All.chi_He1;
  return - 2.0 * chi_He1 / k / T / T / T;
}

// Doubly ionized helium

double d2ln_z_Hep2(double T)
{
  return - 3.0 / 2.0 + T * T * d2ln_z_Hep2_elec(T);
}

double d2ln_z_Hep2_elec(double T)
{
  double k = All.k;
  double chi_He1 = All.chi_He1, chi_He2 = All.chi_He2;
  return - 2.0 * (chi_He1 + chi_He2) / k / T / T / T;
}

// Electron

double d2ln_z_e(double T)
{
  return - 3.0 / 2.0;
}

