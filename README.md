# protostellar_eos
Custom implementation of the equation of state for protostar formation used and described by Tomida et al. 2013 (https://arxiv.org/abs/1206.3567) with the update from Tomida et al. 2015 (https://arxiv.org/abs/1501.04102). This reimplementation was used in Mayer et al. 2025 (https://arxiv.org/abs/2510.12620). 

The code can be compiled to produce an executable called "eos" by running 'make eos' in the folder containing the code files; the Makefile may need to be modified to change the C++ compiler (here g++) and the linking of the GNU Scientific Library (here via -lgsl). Running 'make clean' removes the executable. 

The default output, which are the equation of state tables used in Mayer et al. 2025, is found in data. "U.txt" contains the internal energy per unit mass, "P.txt" the pressure, "CV.txt" the derivative of the internal energy density w.r.t. temperature (at constant density) and "GammaC.txt" the adiabatic index.
These four quantities are all tabulated as a function of temperature (log10(T_min [K]) = 0.5 to log10(T_max [K]) in steps of log(DeltaT [K]) = 0.01 )and density
