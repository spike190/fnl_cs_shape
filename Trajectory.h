#pragma once
#ifndef TRAJECTORY	
#define TRAJECTORY




#include "includes_and_definitions.h"

struct Wrapper {
	real b[5][5];
	real c[6];
	real dc[6];
};

template<int L> void RKCK(Extent_1 &Trajectories, array_view<real> N_, array_view<real> dN_, array_view<real, 2> Parameters_, Wrapper &constants);	//performs an RKCK integration step for background integration

template<int L> void try_a_step(real* Param_old, real &N, real Nend, real &dN, const Wrapper &constants) restrict(amp);		//attempts an integration step

template<int L>	real Evaluate_ODEs_and_error(real* Param_new, real* Param_old, real N, real dN, const Wrapper &constants) restrict(amp);

template<int L>	
class Trajectory {
	private:

		ofstream NSR;	//file for printing main results
		ofstream time;	//file for printing time dependence
		int clo;	//stores initial time

		int No_Modes;	//number of modes to integrate
		int No_Traj;	//number of trajectories to integrate
		
		//const int L = 4;	//number of SR parameters

		Extent_2 Traj_Modes;	//2D loop space for trajectories and modes

		//underlying data source for array_views
		vector<real> Parameters_v;		
		vector<real> Parameters_0_v;
		vector<real> Initial_N_v;

		vector<real> Parameters_k_v;	
		vector<real> N_k_v;
		vector<real> K_v;


		array_view<real, 2> Parameters;		//background parameters to be integrated over
		array_view<real, 2> Parameters_0;	//stores initial conditions on background parameters
		array_view<real> Initial_N;		//stores initial N_total

		array_view<real, 3> Parameters_k;	//parameters for mode integration
		array_view<real, 2> N_k;		//N to be integrated over for each trajectory mode
		array_view<real, 2> K;			//K depends on each trajectory as H will vary
	
		array_view<const real> Shape;		//stores the shape parameter (beta) of each trajectory
		array_view<const real> cs;		//stores 1/cs^2 for each trajectory

		real A, B, cutoff;		//variables for controlling the integration limits

		Wrapper constants;

	public:

		bool broke = false;		//flags broken trajectories (those with eps > 1?)

		//constructors
		Trajectory(int Total_traj, int Total_modes, vector<real> Shape_v, vector<real> cs_v);

		//destructors
		~Trajectory() { time.close(); NSR.close(); }

		void Generate_Trajectories();	//generates random trajectories using the mersenne twistor

		void integrate(real N_end);		//background integration function, integrates backwards for Initial_N -> N_end
		
};	

//trajectory class constructor
template<int L> Trajectory<L>::Trajectory(int Total_traj, int Total_modes, vector<real> Shape_v, vector<real> cs_v) :

		Parameters_v((L + 1)*Total_traj, 0.0f), Parameters(Total_traj, L + 1, Parameters_v),		//initialises arrays for background variables
		Parameters_0_v((L + 1)*Total_traj, 0.0f), Parameters_0(Total_traj, L + 1, Parameters_v),	//arrays for storing initial conditions

		//array for storing initial N
		Initial_N_v(Total_traj, 0.0f), Initial_N(Total_traj, Initial_N_v),

		//arrays for mode integration
		Parameters_k_v((L + 1)*Total_traj*Total_modes, 0.0f), Parameters_k(Total_traj, Total_modes, L + 1, Parameters_k_v),

		//2D loop space for double loop (loops over trajectories and modes)
		Traj_Modes(Total_traj, Total_modes),

		N_k_v(Total_traj*Total_modes, 0.0f), N_k(Traj_Modes, N_k_v),		//stores N variable for each trajectory mode, no need for underlying data source for N_k?
		K_v(Total_traj*Total_modes, 0.0f), K(Traj_Modes, K_v),		//stores k variable for each trajectory mode

		Shape(Total_traj, Shape_v),		//stores beta variable for each trajectory
		cs(Total_traj, cs_v)			//stores 1/cs^2 variable for each trajectory

{

		time.open("time.txt");	//opens the output files
		NSR.open("NSR.txt");

		clo = clock();		//starts the timer

		//assigns integration constants
		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++)constants.b[i][j] = 0.0f;
			constants.c[i] = 0.0f;
			constants.dc[i] = 0.0f;
		}

		constants.b[0][0] = b21;
		constants.b[1][0] = b31;
		constants.b[1][1] = b32;
		constants.b[2][0] = b41;
		constants.b[2][1] = b42;
		constants.b[2][2] = b43;
		constants.b[3][0] = b51;
		constants.b[3][1] = b52;
		constants.b[3][2] = b53;
		constants.b[3][3] = b54;
		constants.b[4][0] = b61;
		constants.b[4][1] = b62;
		constants.b[4][2] = b63;
		constants.b[4][3] = b64;
		constants.b[4][4] = b65;

		constants.c[0] = c1;
		constants.c[2] = c3;
		constants.c[3] = c4;
		constants.c[5] = c6;

		constants.dc[0] = dc1;
		constants.dc[2] = dc3;
		constants.dc[3] = dc4;
		constants.dc[4] = dc5;
		constants.dc[5] = dc6;


		Parameters.discard_data();	//do not copy over data yet
		Parameters_0.discard_data();
		Parameters_k.discard_data();
		Initial_N.discard_data();
		N_k.discard_data();
		K.discard_data();
//		Shape.discard_data();
//		cs.discard_data();

		No_Modes = Total_modes;		//sets number of Fourier modes to consider
		No_Traj = Total_traj;		//sets number of trajectories to consider

		A = (real)exp(6.0), cutoff = (real)exp(5.0), B = 0.05f;	//sets mode integration limits

		Generate_Trajectories();	//assigns initial conditions for the background

}

//generates random trajectories using the mersenne twistor
template<int L>void Trajectory<L>::Generate_Trajectories() {

	//randomly assign initial condition		
	mt19937_64 generator;					//mersenne twistor generator
	uniform_real_distribution<real> distribution;		//converts random generator to a uniform real distribution

	auto rng = bind(distribution, generator);	//rng() = distribution(generator) = [0,1]

	real suppress = 1.5;

	for(int i = 0; i < No_Traj; i++) {

		//int traj = i/256;

		//planck_rng rng2(traj);
		//rng = rng2;

		//initialises eps
		Parameters_v[i*(L + 1)] = 1.0f;
		Parameters_0_v[i*(L + 1)] = 1.0f;


		Initial_N_v[i] = 20.0f*rng() + 60.0f + logf(A);		//sets total e-foldings

		//initialises lambda_l
		if(L > 1) {
			for(int l = 1; l < L; l++) {
				Parameters_v[i*(L + 1) + l] = (2.0f*rng() - 1.0f)*exp(-suppress*(real)l);
				Parameters_0_v[i*(L + 1) + l] = Parameters_v[i*(L + 1) + l];
			}
		}

		//initialises H
		Parameters_0_v[i*(L + 1) + L] = (0.5f + rng());
		Parameters_v[i*(L + 1) + L] = Parameters_0_v[i*(L + 1) + L];		//initialise H //*4.0f*pi*sqrt(2.0f*pi)*1.0e-5f;

	}

	//copy data to accelerator
	Parameters.refresh();
	Initial_N.refresh();

	//CPU -> GPU: use refresh
	//GPU -> CPU: use synchronise

}

//background integration function
template<int L>void Trajectory<L>::integrate(real N_end) {

	Extent_1 Trajectories(Initial_N.extent);

	array_view<real> dN(Trajectories);		//creates an array_view for dN, no need for an underlying data source for dN

	array_view<int> sum(Trajectories);		//array for counting completed trajectories

	//assigns initial dN
	parallel_for_each(Trajectories, [dN, sum](index<1> Traj)restrict(amp) {dN[Traj] = -0.01f; sum[Traj] = 0; });

	int check_start = 0;		//determines how often, integrator should check for completeness
	int check_freq = 1;			

	if(N_end > 5.0f) {		//checks less often for first stage as the integration occurs over many more e-foldings
		check_start = 30;
		check_freq = 10;
	}
	
	//perform integration
	for(int n = 0; n < 20000; n++) {

		//RKCK step
		RKCK<L>(N_end, Trajectories, Initial_N, dN, Parameters, constants);

		//controls often to check for completion
		if(n >= check_start && n % check_freq == 0) {
			//cout << n << "\t";
			if (check_completed(Trajectories, sum, dN, clo)) break; //checks to see if trajectories are finished, breaks if all trajectories finished
		
		}
	}
}

//performs an RKCK integration step for background integration
template<int L> void RKCK(real N_end, Extent_1 &Trajectories, array_view<real> N_, array_view<real> dN_, array_view<real, 2> Parameters_, Wrapper &constants) {

	//loop over all trajectories
	parallel_for_each(Trajectories, [=](index<1> Traj) restrict(amp) {
		
		real N = N_[Traj];	//reads in N
		
		//reads in parameter values
		real Param_old[L + 1];
		for(int l = 0; l < L + 1; l++) Param_old[l] = Parameters_[Traj[0]][l];
		
		real dN = dN_[Traj];	//reads in dN
		
		//checks to see if trajectory finished and guards against floating point errors
		try_a_step<L>(Param_old, N, N_end, dN, constants);
		//if(N > N_end && Param_old[0] > 1.0e-10f)	try_a_step<L>(Param_old, N, N_end, dN, constants);
		//if(N > N_end && Param_old[0] > 1.0e-10f)	try_a_step<L>(Param_old, N, N_end, dN, constants);
		//if(N > N_end && Param_old[0] > 1.0e-10f)	try_a_step<L>(Param_old, N, N_end, dN, constants);
		
		//else dN = 0.0f;		//mark trajectory as finished	

		
		N_[Traj] = N;	//updates N
		dN_[Traj] = dN;	//update step size
		for(int l = 0; l < L + 1; l++)Parameters_[Traj[0]][l] = Param_old[l];		//updates Parameters
	});

	//cout << N_[0] << "\t" << Parameters_[0][L] << "\t" << Parameters_[0][0] << "\t" << Parameters_[0][1] << "\t" << Parameters_[0][L - 1] * Parameters_[0][L] * pow(Parameters_[0][0], 0.5*(2-L)) << endl;
}

//attempts an integration step
template<int L> void try_a_step(real* Param_old, real &N, real Nend, real &dN, const Wrapper &constants) restrict(amp) {

	real Param_new[L + 1];

	real err_max = Evaluate_ODEs_and_error<L>(Param_new, Param_old, N, dN, constants);		//N is passed by value here so it does not change

	//adaptive step size******************************

	if(err_max < 1.0f) {
		//if step successful...
		N = N + dN;		//new value of N		

		//calculate new step size
		if(err_max > ERRCON) dN = SAFETY*dN*powf(err_max, PGROW);	//increase step size by no more than a factor of 5
		else dN = 5.0f*dN;

		//checks to see if integrator shoots past final integration point on next step
		if(N + dN < Nend) dN = Nend - N;

		if(dN < -dN_cap) dN = -dN_cap;

		for(int l = 0; l < L + 1; l++)Param_old[l] = Param_new[l];		//update N and Parameters

	}
	else {
		//error too large, reduce step size but by no more than a factor of 10
		real dNtemp = SAFETY*dN*powf(err_max, PSHRNK);
		if(dNtemp*dNtemp < dN*dN*0.01f)	dN = dN*0.1f;
		else dN = dNtemp;
	};


	//constant step size***************************
	/*N = N + dN;
	
	//checks to see if integrator shoots past final integration point on next step
	if (N + dN < Nend) dN = Nend - N;

	for (int l = 0; l < L + 1; l++)Param_old[l] = Param_new[l];		//update N and Parameters*/

}


//obtains new updated functions and the worst offending error
template<int L>	real Evaluate_ODEs_and_error(real* Param_new, real* Param_old, real N, real dN, const Wrapper &coeff) restrict(amp) {

	real ak[6 * (L + 1)];

	//performs the 6 function evaluations necessary for integration

	ODE(Param_old, ak, 0, L);
	for(int l = 0; l < (L + 1); l++)Param_new[l] = Param_old[l] + (b21 * ak[l])*dN;
	N = N + 0.2f*dN;	//a2

	ODE(Param_new, ak, 1, L);
	for(int l = 0; l < (L + 1); l++)Param_new[l] = Param_old[l] + (b31 * ak[l] + b32 * ak[(L + 1) + l])*dN;
	N = N + (0.1f)*dN;	//a3 - a2

	ODE(Param_new, ak, 2, L);
	for(int l = 0; l < (L + 1); l++)Param_new[l] = Param_old[l] + (b41 * ak[l] + b42 * ak[(L + 1) + l] + b43 * ak[2 * (L + 1) + l])*dN;
	N = N + (0.3f)*dN;	//a4 - a3

	ODE(Param_new, ak, 3, L);
	for(int l = 0; l < (L + 1); l++)Param_new[l] = Param_old[l] + (b51 * ak[l] + b52 * ak[(L + 1) + l] + b53 * ak[2 * (L + 1) + l] + b54 * ak[3 * (L + 1) + l])*dN;
	N = N + (0.4f)*dN;	//a5 - a4

	ODE(Param_new, ak, 4, L);
	for(int l = 0; l < (L + 1); l++)Param_new[l] = Param_old[l] + (b61 * ak[l] + b62 * ak[(L + 1) + l] + b63 * ak[2 * (L + 1) + l] + b64 * ak[3 * (L + 1) + l] + b65 * ak[4 * (L + 1) + l])*dN;
	N = N + (-0.125f)*dN;	//a6 - a5

	ODE(Param_new, ak, 5, L);

	
	real Param_error[L + 1];

	//calculates errors and finds worst offender
	real err_max = 0.0f;
	
	for(int l = 0; l < L + 1; l++) {

		Param_new[l] = Param_old[l] + (ak[l] * c1 + ak[2 * (L + 1) + l] *c3 + ak[3 * (L + 1) + l] * c4 + ak[5 * (L + 1) + l] * c6)*dN;

		Param_error[l] = (ak[l] * dc1 + ak[2 * (L + 1) + l] *dc3 + ak[3 * (L + 1) + l] * dc4 + ak[4 * (L + 1) + l] * dc5 + ak[5 * (L + 1) + l] * dc6)*dN;

		if(Param_old[l] != 0.0f)Param_error[l] = fabsf(Param_error[l] / Param_old[l]);
		else Param_error[l] = 0.0f;
		if(Param_error[l] > err_max)err_max = Param_error[l];
	}


	err_max = err_max*tol_1;		//scales error relative to tolerance

	return err_max;

}



#endif

