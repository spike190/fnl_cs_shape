#pragma once
#ifndef FUNCTION_DECLARATIONS	
#define FUNCTION_DECLARATIONS

#include "includes_and_definitions.h"

void get_accelerators();

int Accumulate(array_view<int> sum);	//returns the sum of a list of integers

bool check_completed(Extent_1 Trajectories, array_view<int> sum, array_view<real> dN, int clo);

void Update_functions(real* New, real* Old, real* ak, real &N, real dN, int i, int L)restrict(amp);

void ODE(real* Parameters, real* ak, int i, int L) restrict(amp);
real dEps_dN(real eps, real eta) restrict(amp);
real dLambda_dN(real eps, real eta, real lambda, real lambda1, real l)restrict(amp);
real dLambda_dN(real eps, real eta, real lambda, real l)restrict(amp);

real dH_dN(real eps, real H) restrict(amp);

#endif