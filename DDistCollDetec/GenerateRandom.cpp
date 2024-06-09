#include "pch.h"
#include "config.h"
#include "GenerateRandom.h"

void GenerateRandomRotTra(double R[][3], double *t, double interval_boundary) {
	t[0] = RandT<double>(-interval_boundary, interval_boundary);
	t[1] = RandT<double>(-interval_boundary, interval_boundary);
	t[2] = RandT<double>(-interval_boundary, interval_boundary);

	double theta_Z = RandT<double>(-PI, PI);
	double theta_Y = RandT<double>(-PI, PI);
	double theta_X = RandT<double>(-PI, PI);

	double cz = cos(theta_Z), sz = sin(theta_Z);
	double cy = cos(theta_Y), sy = sin(theta_Y);
	double cx = cos(theta_X), sx = sin(theta_X);

	R[0][0] = cy * cz, R[0][1] = -cy * sz, R[0][2] = sy;
	R[1][0] = cx * sz - sx * sy*cz, R[1][1] = cx * cz + sx * sy*sz, R[1][2] = sx * cy;
	R[2][0] = -sx * sz - cx * sy*cz, R[2][1] = -sx * cz + cx * sy*sz, R[2][2] = cx * cy;
}

void GenerateRandomVel(double *v, double velocity_boundary) {
	v[0] = RandT<double>(-velocity_boundary, velocity_boundary);
	v[1] = RandT<double>(-velocity_boundary, velocity_boundary);
	v[2] = RandT<double>(-velocity_boundary, velocity_boundary);
}
