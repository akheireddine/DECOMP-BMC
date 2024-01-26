
#pragma once

#include <vector>
#include <string>

#include "../../nusmv_wrapper.hh"
#include "utils/EnvBMC.hh"

class NuSMV
{
public:
	NuSMV() = default;

	// create new SMV problem in inputFile
	NuSMV(const char *inputFile, EnvBMC *envi);

	~NuSMV();

	// Getter for k_ variable
	inline int get_k() { return k_; };

	// Setter for k_ variable
	inline void set_k(int k) { k_ = k; };

	// NuSMV read model
	void read_model();

	/* generate DIMACS file with specific ltl type :
		prop_value : 0 no property (add a true property)
		prop_value : 1 generate a random property (TODO)
		prop_value : 2 use current property defined in smv file */
	void generate_cnf_formula(int k_val, int prop_value);

	// Executes NuSMV commands from file source (batch mode)
	void run_NuSMV(std::string source);

	// Print properties in SMV program
	void print_properties();

protected:
	// NuSMV Envrionment
	NuSMVEnv_ptr env_nusmv;

	// Filename of input model
	std::string modelFile_;

	// Bound of the transition relation
	int k_;

	EnvBMC *env;
};