#include "proto.h"

double get_partition_function(int i, double T)
{
  if(i == 0) return z_H2(T);
  else if(i == 1) return z_H(T);
  else if(i == 2) return z_Hp(T);
  else if(i == 3) return z_He(T);
  else if(i == 4) return z_Hep(T);
  else if(i == 5) return z_Hep2(T);
  else if(i == 6) return z_e(T);
  else return 0.0;
}
 
// Molecular hydrogen
 
double z_H2(double T)
{
  return z_H2_tr(T) * z_H2_rot(T) * z_H2_vib(T) * 4.0 * 2.0;
}

double z_H2_tr(double T)
{
  double k = All.k, h = All.h;
  double m_H2 = 2.0 * All.mu_H * All.mP;
  return pow(2.0 * M_PI * m_H2 * k * T , 3.0 / 2.0) / pow(h , 3.0);
}

double z_H2_rot(double T)
{
  double theta_rot = All.theta_rot;
  double z_even = 0.0, z_odd = 0.0;
  
  for(int j = 0 ; j < 10000 ; j++)
  {
    if(j % 2 == 0) z_even += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T);
    else           z_odd  += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T);
  }
  
  return pow(z_even , 1.0 / (All.op_ratio + 1.0)) * pow(3.0 * z_odd * exp(theta_rot / T) , All.op_ratio / (All.op_ratio + 1.0));
}

double z_H2_vib(double T) /* Appears to be a mistake in Tomida 2013; this is the standard vibrational partition function. */
{
  double theta_vib = All.theta_vib;
  return 1.0 / (1.0 - exp(- theta_vib / T));
}

// Atomic hydrogen

double z_H(double T)
{
  return z_H_tr(T) * 2.0 * z_H_elec(T);
}

double z_H_tr(double T)
{
  double k = All.k, h = All.h;
  double m_H = All.mP * All.mu_H;
  return pow(2.0 * M_PI * m_H * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

double z_H_elec(double T)
{
  double k = All.k;
  double chi_diss = All.chi_diss;
  return 2.0 * exp( - chi_diss / 2.0 / k / T);
}

// Ionized hydrogen

double z_Hp(double T)
{
  return z_Hp_tr(T) * 2.0 * z_Hp_elec(T);
}

double z_Hp_tr(double T)
{
  double k = All.k, h = All.h;
  double m_Hp = All.mP;
  return pow(2.0 * M_PI * m_Hp * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

double z_Hp_elec(double T)
{
  double k = All.k;
  double chi_diss = All.chi_diss, chi_ion = All.chi_ion;
  return 2.0 * exp( - (chi_diss + 2.0 * chi_ion) / 2.0 / k / T);
}

// Atomic helium

double z_He(double T)
{
  return z_He_tr(T);
}

double z_He_tr(double T)
{
  double k = All.k, h = All.h;
  double m_He = All.mu_He / All.mu_H * All.mP;
  return pow(2.0 * M_PI * m_He * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

// Ionized helium

double z_Hep(double T)
{
  return z_Hep_tr(T) * z_Hep_elec(T);
}

double z_Hep_tr(double T)
{
  double k = All.k, h = All.h;
  double m_Hep =  All.mu_He / All.mu_H * All.mP - All.me;
  return pow(2.0 * M_PI * m_Hep * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

double z_Hep_elec(double T)
{
  double k = All.k;
  double chi_He1 = All.chi_He1;
  return 2.0 * exp( - chi_He1 / k / T);
}

// Doubly ionized helium

double z_Hep2(double T)
{
  return z_Hep2_tr(T) * z_Hep2_elec(T);
}

double z_Hep2_tr(double T)
{
  double k = All.k, h = All.h;
  double m_Hep2 =  All.mu_He / All.mu_H * All.mP - 2.0 * All.me;;
  return pow(2.0 * M_PI * m_Hep2 * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

double z_Hep2_elec(double T)
{
  double k = All.k;
  double chi_He1 = All.chi_He1, chi_He2 = All.chi_He2;
  return 2.0 * exp( - (chi_He1 + chi_He2) / k / T);
}

// Electron

double z_e(double T)
{
  return z_e_tr(T) * 2.0;
}

double z_e_tr(double T)
{
  double k = All.k, h = All.h;
  double m_e = All.me;
  return pow(2.0 * M_PI * m_e * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

