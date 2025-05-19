#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "partition_functions.cc"
#include "partition_derivatives.cc"
#include "second_derivatives.cc"
#include "calculate_abundances.cc"
#include "density_derivatives.cc"

using namespace std;

int main()
{ 
  double X = 0.72, Y = 0.28;
  double k = 1.38065e-16;
  
  double rho, T;
  double n, mu, U, gammaE, C_V, S_T, S_rho, P, P_T, P_rho, c_s, n_T, n_rho, gammaC;
  
  double abundances[7];    
  double dln_abundances_dln_rho[7], dln_abundances_dln_T[7];
  
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
	P = 0.0;
	U = 0.0;

	rho = pow(10.0 , logrho);

	//cout << "\n\n\n";
	//cout << rho << "\n\n\n";
    
	n = 0.0;

	calculate_abundances(X, Y, rho, T, abundances);
	for(int i = 0 ; i < 7 ; i++) n += abundances[i];

	for(int i = 0 ; i < 7 ; i++) U += abundances[i] * get_partition_derivative(i,T);
	U *= k * T / rho; 

	P = n * k * T;  
      
	dln_n_dln_rho(X, Y, rho, abundances, dln_abundances_dln_rho);
	dln_n_dln_T(T, abundances, dln_abundances_dln_T);

	//for(int i = 0 ; i < 7 ; i++) cout << dln_abundances_dln_T[i] << "\n";

	S_T = 0.0, S_rho = 0.0;

	for(int i = 0 ; i < 7 ; i++) S_T += k / rho / T * abundances[i] * (get_partition_second_derivative(i,T) + (1.0 + dln_abundances_dln_T[i]) * get_partition_derivative(i,T)); 
	for(int i = 0 ; i < 7 ; i++) S_rho += k / rho / rho * abundances[i] * ((dln_abundances_dln_rho[i] - 1.0) * get_partition_derivative(i,T) - 1.0);
	
	C_V = 0.0;

	C_V = rho * T * S_T;

	//cout << C_V / n / k << "\n";
	//cout << U / T * rho / n / k << "\n";

	n_T = 0.0, n_rho = 0.0;

	for(int i = 0 ; i < 7 ; i++)
	{
	  n_T += abundances[i] / n * dln_abundances_dln_T[i];
	  n_rho += abundances[i] / n * dln_abundances_dln_rho[i];
	}

	P_T = 1.0 + n_T;
	P_rho = n_rho;

	c_s = sqrt(P / rho * P_rho - P / T * P_T * S_rho / S_T);

	gammaC = rho / P * c_s * c_s;

	//cout << c_s << "\n";

	//cout << gammaC << "\n";

	if((gammaC < 1.0) || (gammaC > 1.67)) terminate();
	
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
