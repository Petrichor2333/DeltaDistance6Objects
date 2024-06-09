#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <chrono>

#include "f2c.h"        
#include "clapack.h"

#pragma once

using namespace std;

constexpr auto PI = 3.141592653589793;
constexpr auto ERROR_GRAD = 0.00000001;
constexpr auto ERROR_COLLIDE = 0.00000001;
constexpr auto ERROR_COLLIDE_V = 0.5*ERROR_COLLIDE*ERROR_COLLIDE;

double epsilon = std::numeric_limits<double>::epsilon();

template<typename T>
T RandT(T _min, T _max)
{
	T temp;
	if (_min > _max)
	{
		temp = _max;
		_max = _min;
		_min = temp;
	}
	return rand() / (double)RAND_MAX *(_max - _min) + _min;
}


#define NDEBUG_SOL