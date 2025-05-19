#include "proto.h"

using namespace std;

GlobalVariables All;

int main()
{ 
  All.X = 0.72; // Hydrogen mass-fraction
  All.Y = 0.28; // Helium mass-fraction
  All.op_ratio = 3.0; // Ratio of ortho- to parahydrogen
  
  double rho, T;
  double n, mu, U, gammaE, C_V, S_T, S_rho, P, P_T, P_rho, c_s, n_T, n_rho, gammaC;
  
  double abundances[7], dln_abundances_dln_rho[7], dln_abundances_dln_T[7];
  
  ofstream file_P ("./data/P.txt");
  ofstream file_U ("./data/U.txt");
  ofstream file_Cv ("./data/Cv.txt");
  ofstream file_gammaC ("./data/GammaC.txt");
  
  for(double logT = 0.5 ; logT < 6.01 ; logT += 0.02) 
  {
    T = pow(10.0 , logT);
    cout << "\n\n\n";
    cout << T << "\n\n\n";   
    for(double logrho = -22.0 ; logrho <= 1.01 ; logrho += 0.05)
    {
    	rho = pow(10.0 , logrho);
    	
	P = 0.0, U = 0.0, C_V = 0.0, gammaC = 0.0;
	n = 0.0, S_T = 0.0, S_rho = 0.0, n_T = 0.0, n_rho = 0.0, P_T = 0.0, P_rho = 0.0;

	calculate_abundances(rho, T, abundances);
	
	for(int i = 0 ; i < 7 ; i++) n += abundances[i];
	for(int i = 0 ; i < 7 ; i++) U += abundances[i] * get_partition_derivative(i,T);
	
	U *= All.k * T / rho; 
	P = n * All.k * T;  
      
	dln_n_dln_rho(rho, abundances, dln_abundances_dln_rho);
	dln_n_dln_T(T, abundances, dln_abundances_dln_T);

	for(int i = 0 ; i < 7 ; i++)
	{
	  S_T += All.k / rho / T * abundances[i] * (get_partition_second_derivative(i,T) + (1.0 + dln_abundances_dln_T[i]) * get_partition_derivative(i,T)); 
	  S_rho += All.k / rho / rho * abundances[i] * ((dln_abundances_dln_rho[i] - 1.0) * get_partition_derivative(i,T) - 1.0);
	  n_T += abundances[i] / n * dln_abundances_dln_T[i];
	  n_rho += abundances[i] / n * dln_abundances_dln_rho[i];
	}

	P_T = 1.0 + n_T;
	P_rho = n_rho;

	C_V = rho * T * S_T;
	c_s = sqrt(P / rho * P_rho - P / T * P_T * S_rho / S_T);
	gammaC = rho / P * c_s * c_s;

	if((gammaC < 1.0) || (gammaC > 1.67)) terminate(); // Should not happen anymore
	
	file_P << setprecision(20) << log10(P) << "   ";
        file_U << setprecision(20) << log10(U) << "   ";
	file_Cv << setprecision(20) << log10(C_V) << "   ";
        file_gammaC << setprecision(20) << gammaC << "   ";
    }
    file_P << "\n";
    file_U << "\n";
    file_Cv << "\n";
    file_gammaC << "\n";
  }
  file_P.close();
  file_U.close();
  file_Cv.close();
  file_gammaC.close();
    
  return 0;
}
