#include <cmath>

double K_dis(double T);
double K_ion(double T);
double K_He1(double T);
double K_He2(double T);

void calculate_abundances(double X, double Y, double rho, double T, double abundances[7]);
double solve_for_electron_abundance(double nH_tot, double nHe_tot, double K_dis, double K_ion, double K_He1, double K_He2);
double f_ne(double ne, double nH_tot, double nHe_tot, double K_dis, double K_ion, double K_He1, double K_He2);

double K_dis(double T)
{
  return pow(z_H(T) , 2.0) / z_H2(T);
}

double K_ion(double T)
{
  return z_Hp(T) * z_e(T) / z_H(T);
}

double K_He1(double T)
{
  return z_Hep(T) * z_e(T) / z_He(T);
}

double K_He2(double T)
{
  return z_Hep2(T) * z_e(T) / z_Hep(T);
}

void calculate_abundances(double X, double Y, double rho, double T, double abundances[7])
{
  double ne, nH, nHe, nH2, nHp, nHep, nHep2;
  
  double Kdis = fmax(1e-200 , K_dis(T)), Kion = fmax(1e-200 , K_ion(T)), KHe1 = fmax(1e-200 , K_He1(T)), KHe2 = fmax(1e-200 , K_He2(T)); 
  
  double nH_tot = rho * X / (1.67262178e-24 + 9.1093829e-28), nHe_tot = rho * Y / (4.002602 / 1.007276466621 * 1.67262178e-24);
  
  ne = solve_for_electron_abundance(nH_tot, nHe_tot, Kdis, Kion, KHe1, KHe2);
  
  nH = 2.0 * nH_tot / ((1.0 + Kion / ne) + sqrt((1.0 + Kion / ne) * (1.0 + Kion / ne) + 4.0 * 2.0 / Kdis * nH_tot));
  nHp = Kion * nH / ne;

  nHe = nHe_tot / (1.0 + KHe1 / ne + KHe1 * KHe2 / ne / ne);
  nH2 = nH * nH / Kdis;

  nHep = KHe1 * nHe / ne;
  nHep2 = KHe2 * nHep / ne;
  
  abundances[0] = fmax(nH2,1e-200);
  abundances[1] = fmax(nH,1e-200);
  abundances[2] = fmax(nHp,1e-200);
  abundances[3] = fmax(nHe,1e-200);
  abundances[4] = fmax(nHep,1e-200);
  abundances[5] = fmax(nHep2,1e-200);
  abundances[6] = fmax(ne,1e-200);
}

double solve_for_electron_abundance(double nH_tot, double nHe_tot, double Kdis, double Kion, double KHe1, double KHe2)
{
  int iter;
  double err, dne;
  double ne_old = nH_tot;
  
  double ne = ne_old, ne_lower = ne / 10.0, ne_upper = ne * 10.0;
  
  while(f_ne(ne_lower, nH_tot, nHe_tot, Kdis, Kion, KHe1, KHe2) * f_ne(ne_upper, nH_tot, nHe_tot, Kdis, Kion, KHe1, KHe2) > 0.0)
  {
    ne_lower /= 10.0;
    ne_upper *= 10.0;
  }

  do
    {
      ne = 0.5 * (ne_lower + ne_upper);

      err = f_ne(ne, nH_tot, nHe_tot, Kdis, Kion, KHe1, KHe2);

      if(err * f_ne(ne_lower, nH_tot, nHe_tot, Kdis, Kion, KHe1, KHe2) > 0.0)
        {
          ne_lower = ne;
        }
      else
        {
          ne_upper = ne;
        }

      dne = ne_upper - ne_lower;
      iter++;
    }
  while(fabs(dne / ne) > 1.0e-10);  
  
  return ne;
}

double f_ne(double ne, double nH_tot, double nHe_tot, double Kdis, double Kion, double KHe1, double KHe2)
{
  return (2.0 * ne * ne * nH_tot * Kion) / (sqrt((ne + Kion) * (ne + Kion) + 8.0 / Kdis * nH_tot * ne * ne) + ne + Kion) + (KHe1 * ne + 2.0 * KHe1 * KHe2) / (ne * ne + KHe1 * ne + KHe1 * KHe2) * nHe_tot * ne * ne - ne * ne * ne;
}
