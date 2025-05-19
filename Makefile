eos: main.cc calculate_abundances.cc partition_functions.cc partition_derivatives.cc second_derivatives.cc density_derivatives.cc proto.h
	g++ main.cc calculate_abundances.cc partition_functions.cc partition_derivatives.cc second_derivatives.cc density_derivatives.cc -lgsl -o eos
clean:
	rm eos
