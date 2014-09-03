#include "includes_and_definitions.h"

void get_accelerators(){

	vector<accelerator> accs = accelerator::get_all();
	std::for_each(accs.begin(), accs.end(), [](accelerator acc)
	{
		std::wcout << "New accelerator: " << acc.description << std::endl;
		std::wcout << "device_path = " << acc.device_path << std::endl;
		std::wcout << "version = " << (acc.version >> 16) << '.' << (acc.version & 0xFFFF) << std::endl;
		std::wcout << "dedicated_memory = " << acc.dedicated_memory << " KB" << std::endl;
		std::wcout << "doubles = " << ((acc.supports_double_precision) ? "true" : "false") << std::endl;
		std::wcout << "limited_doubles = " << ((acc.supports_limited_double_precision) ? "true" : "false") << std::endl;
		std::wcout << "has_display = " << ((acc.has_display) ? "true" : "false") << std::endl;
		std::wcout << "is_emulated = " << ((acc.is_emulated) ? "true" : "false") << std::endl;
		std::wcout << "is_debug = " << ((acc.is_debug) ? "true" : "false") << std::endl;
		std::cout << std::endl;
	});

	//accelerator::set_default(accelerator::direct3d_warp);
}


//counts completed trajectories
//could download and use the c++ amp algorithm library but probably not worth it. 
int Accumulate(array_view<int> sum) {
	
	vector<int> sum_total(1, 0);

	int Total_elements = sum.extent[0];

	//reduction algorithm to add up a list of numbers
	for(int stride = (Total_elements / 2); stride > 0; stride /= 2) {
		parallel_for_each(Extent_1(stride), [stride, sum](index<1> idx) restrict(amp) {
			int Traj = idx[0];
			sum[Traj] += sum[Traj + stride];		//counts completed trajectories
		});
	}

	copy(sum.section(0, 1), sum_total.begin());		//section(origin_x, length_x)
	sum.discard_data();
	return sum_total[0];
}

bool check_completed(Extent_1 Trajectories, array_view<int> sum, array_view<real> dN, int clo){

	//copies labels over to array for counting completed traj and reordering
	parallel_for_each(Trajectories, [dN, sum](index<1> Traj) restrict(amp) {

		if (dN[Traj] == 0.0f) sum[Traj] = 1;		//labels completed trajectories
		else sum[Traj] = 0;

	});

	int traj_completed = Accumulate(sum);		//adds up completed trajectories

	//cout << traj_completed << "\t" << Trajectories[0] << "\t" << (double)(clock() - clo)*0.001 << endl;

	if (traj_completed == Trajectories[0]) return true;		//all trajectories finished
	else return false;
}


void Update_functions(real* New, real* Old, real* ak, real &N, real dN, int i, int L)restrict(amp) {

	if(i == 0) {
		N = N + 0.2f*dN;
		for(int l = 0; l < L; l++)New[l] = Old[l] + b21*ak[l] * dN;
	}
	else if(i == 1) {
		N = N + (0.1f)*dN;
		for(int l = 0; l < L; l++)New[l] = Old[l] + (b31*ak[l] + b32*ak[L + l])*dN;
	}
	else if(i == 2) {
		N = N + (0.3f)*dN;
		for(int l = 0; l < L; l++)New[l] = Old[l] + (b41*ak[l] + b42*ak[L + l] + b43*ak[2 * L + l])*dN;
	}
	else if(i == 3) {
		N = N + (0.4f)*dN;
		for(int l = 0; l < L; l++)New[l] = Old[l] + (b51*ak[l] + b52*ak[L + l] + b53*ak[2 * L + l] + b54*ak[3 * L + l])*dN;
	}
	else if(i == 4) {
		N = N + (-0.125f)*dN;
		for(int l = 0; l < L; l++)New[l] = Old[l] + (b61*ak[l] + b62*ak[L + l] + b63*ak[2 * L + l] + b64*ak[3 * L + l] + b65*ak[4 * L + l])*dN;
	}

}

//evalutes ODE's and gives results to ak
void ODE(real* Parameters, real* ak, int i, int L) restrict(amp) {

	//Parameters[0] = eps
	//Parameters[L] = H

	real eta = 0.0f;

	//evolution of lambda taking care of truncated hierarchy
	if(L > 1) {
		eta = Parameters[1];
		ak[i*(L + 1) + L - 1] = dLambda_dN(Parameters[0], eta, Parameters[L - 1], (real)(L - 1));		//final lambda value

		if(L > 2) {

			for(int l = 1; l < L - 1; l++) ak[i*(L + 1) + l] = dLambda_dN(Parameters[0], eta, Parameters[l], Parameters[l + 1], (real)l);

		}
	}

	ak[i*(L + 1)] = dEps_dN(Parameters[0], eta);					//evolution for eps
	ak[i*(L + 1) + L] = dH_dN(Parameters[0], Parameters[L]);		//evolution for H	
}

real dEps_dN(real eps, real eta) restrict(amp) {
	return 2.0f*eps*(eps - eta);
}

real dLambda_dN(real eps, real eta, real lambda, real lambda1, real l)restrict(amp) {
	return ((l*eps + (1.0f - l)*eta)*lambda - lambda1);
}

real dLambda_dN(real eps, real eta, real lambda, real l)restrict(amp) {
	return ((l*eps + (1.0f - l)*eta)*lambda);
}

real dH_dN(real eps, real H) restrict(amp) {
	return -eps*H;
}