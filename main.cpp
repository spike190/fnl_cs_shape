#include "includes_and_definitions.h"

int main() {


	//get_accelerators();

	//accelerator::set_default(accelerator::direct3d_warp);

	const int L = 4;

	int Total_Traj = 65536*8;
	int Total_modes = 64;

	vector<real> Shape(Total_Traj, (real)1.0);	//defines shape and sound speed for trajectories
	vector<real> cs(Total_Traj, (real)1.0);		//will probably never be randomly sampled but may be looped over to analyse functional dependance


	Trajectory<L> Trajectories(Total_Traj, Total_modes, Shape, cs);		//generates the trajectories class with L HSR parameters for the desired number of trajectories and Fourier modes with the desired shape and sound speeds

	Trajectories.integrate((real)0.0);	//integrates backwards to N = 0

	//assign initial conditions on mode functions

}