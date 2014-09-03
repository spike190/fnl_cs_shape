#pragma once
#ifndef INCLUDE_DECLARATIONS

	#define INCLUDE_DECLARATIONS

	#include <vector>
	#include <stdlib.h>
	#include <cmath>
	#include <fstream>
	#include <iostream>
	#include <time.h>
	#include <random>
	#include <functional>
	#include <amp.h>
	#include <amp_math.h>
	#include <amp_short_vectors.h>

	using namespace std;
	using namespace concurrency;
	using namespace concurrency::fast_math;

	typedef float real;
	typedef concurrency::extent<1> Extent_1;
	typedef concurrency::extent<2> Extent_2;

	#include "integration constants.h"
	#include "Trajectory.h"
	#include "declarations.h"

	const real A = exp(6.0f);			//0.05	mode integration always starts at times such that k = A*aH and ends at times such that k = BaH
	const real cutoff = exp(5.0f);			//fnl cutoff time
	const real B = 0.05f;

	const real h = 0.7f;				//usual little h in cosmology

	//numbers needed in ns r formulae
	#define EulerConst 0.57721566490153286060651209008240243104215933593992f	//euler-mascheroni constant
	#define C (4.0f*(log(2.0f) + EulerConst) - 5.0f)				//= 0.08..
	#define pi 3.14159265358979323846264338327950288419716939937510582f		


	

		
#endif // !INCLUDE_DECLARATIONS


