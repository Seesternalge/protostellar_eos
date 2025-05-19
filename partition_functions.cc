#include <cmath>

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


// Molecular hydrogen

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
 
double z_H2(double T)
{
  return z_H2_tr(T) * z_H2_rot(T) * z_H2_vib(T) * 4.0 * 2.0;
}

double z_H2_tr(double T)
{
  double k = 1.38065e-16, h = 6.6260695e-27;
  double m_H2 = 2.0 * (1.67262178e-24 + 9.1093829e-28);
  return pow(2.0 * M_PI * m_H2 * k * T , 3.0 / 2.0) / pow(h , 3.0);
}

double z_H2_rot(double T)
{
  double theta_rot = 170.64;
  double z_even = 0.0, z_odd = 0.0;
  
  for(int j = 0 ; j < 10000 ; j++)
  {
    if(j % 2 == 0) z_even += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T);
    else           z_odd  += (2.0 * j + 1.0) * exp( - j * (j + 1.0) * theta_rot / 2.0 / T);
  }
  
  return pow(z_even , 1.0 / 4.0) * pow(3.0 * z_odd * exp(theta_rot / T) , 3.0 / 4.0);
}

double z_H2_vib(double T) /* Appears to be a mistake in Tomida 2013; this is the standard vibrational partition function. However, might need to change theta_vib by a factor of 2... */
{
  double theta_vib = 5984.48;
  return 1.0 / (1.0 - exp(- theta_vib / T));
}

// Atomic hydrogen

double z_H(double T)
{
  return z_H_tr(T) * 2.0 * z_H_elec(T);
}

double z_H_tr(double T)
{
  double k = 1.38065e-16, h = 6.6260695e-27;
  double m_H = 1.67262178e-24 + 9.1093829e-28;
  return pow(2.0 * M_PI * m_H * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

double z_H_elec(double T)
{
  double k = 1.38065e-16;
  double chi_diss = 7.17e-12;
  return 2.0 * exp( - chi_diss / 2.0 / k / T);
}

// Ionized hydrogen

double z_Hp(double T)
{
  return z_Hp_tr(T) * 2.0 * z_Hp_elec(T);
}

double z_Hp_tr(double T)
{
  double k = 1.38065e-16, h = 6.6260695e-27;
  double m_Hp = 1.67262178e-24;
  return pow(2.0 * M_PI * m_Hp * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

double z_Hp_elec(double T)
{
  double k = 1.38065e-16;
  double chi_diss = 7.17e-12, chi_ion = 2.18e-11;
  return 2.0 * exp( - (chi_diss + 2.0 * chi_ion) / 2.0 / k / T);
}

// Atomic helium

double z_He(double T)
{
  return z_He_tr(T);
}

double z_He_tr(double T)
{
  double k = 1.38065e-16, h = 6.6260695e-27;
  double m_He = 4.002602 / 1.007276466621 * 1.67262178e-24;
  return pow(2.0 * M_PI * m_He * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

// Ionized helium

double z_Hep(double T)
{
  return z_Hep_tr(T) * z_Hep_elec(T);
}

double z_Hep_tr(double T)
{
  double k = 1.38065e-16, h = 6.6260695e-27;
  double m_Hep =  4.002602 / 1.007276466621 * 1.67262178e-24 - 9.1093829e-28;
  return pow(2.0 * M_PI * m_Hep * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

double z_Hep_elec(double T)
{
  double k = 1.38065e-16;
  double chi_He1 = 3.94e-11;
  return 2.0 * exp( - chi_He1 / k / T);
}

// Doubly ionized helium

double z_Hep2(double T)
{
  return z_Hep2_tr(T) * z_Hep2_elec(T);
}

double z_Hep2_tr(double T)
{
  double k = 1.38065e-16, h = 6.6260695e-27;
  double m_Hep2 =  4.002602 / 1.007276466621 * 1.67262178e-24 - 2.0 * 9.1093829e-28;
  return pow(2.0 * M_PI * m_Hep2 * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

double z_Hep2_elec(double T)
{
  double k = 1.38065e-16;
  double chi_He1 = 3.94e-11, chi_He2 = 8.72e-11;
  return 2.0 * exp( - (chi_He1 + chi_He2) / k / T);
}

// Electron

double z_e(double T)
{
  return z_e_tr(T) * 2.0;
}

double z_e_tr(double T)
{
  double k = 1.38065e-16, h = 6.6260695e-27;
  double m_e = 9.1093829e-28;
  return pow(2.0 * M_PI * m_e * k * T , 3.0 / 2.0) / pow(h , 3.0); 
}

