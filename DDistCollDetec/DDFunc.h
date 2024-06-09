#pragma once

// Value of V_{\Delta}(x)
// Gradient of V_{\Delta}(x)
// Hessian matrix of V_{\Delta}(x)

void DDFuncSqGrad_Obj_EL(const double* P, double* Grad, void* PObjData);
void DDFuncSqHess_Obj_EL(double* Hess, void* PObjData);
void DDFuncSqValu_Obj_EL(double* val, void* PObjData);
void DDFuncSqValuIndep_Obj_EL(const double* P, double* val, void* PObjData);

void DDFuncSqGrad_Obj1_MS(const double* P, double* Grad, void* PObjData);
void DDFuncSqHess_Obj1_MS(double* Hess, void* PObjData);
void DDFuncSqValu_Obj1_MS(double *val, void* PObjData);
void DDFuncSqValuIndep_Obj1_MS(const double* P, double *val, void* PObjData);

void DDFuncSqGrad_Obj2_PC(const double* P, double* Grad, void* PObjData);
void DDFuncSqHess_Obj2_PC(double* Hess, void* PObjData);
void DDFuncSqValu_Obj2_PC(double *val, void* PObjData);
void DDFuncSqValuIndep_Obj2_PC(const double* P, double *val, void* PObjData);

void DDFuncSqGrad_Obj3_DM(const double* P, double* Grad, void* PObjData);
void DDFuncSqHess_Obj3_DM(double* Hess, void* PObjData);
void DDFuncSqValu_Obj3_DM(double *val, void* PObjData);
void DDFuncSqValuIndep_Obj3_DM(const double* P, double *val, void* PObjData);

void DDFuncSqGrad_Obj4_CK(const double* P, double* Grad, void* PObjData);
void DDFuncSqHess_Obj4_CK(double* Hess, void* PObjData);
void DDFuncSqValu_Obj4_CK(double *val, void* PObjData);
void DDFuncSqValuIndep_Obj4_CK(const double* P, double *val, void* PObjData);

void DDFuncSqGrad_Obj5_PL(const double* P, double* Grad, void* PObjData);
void DDFuncSqHess_Obj5_PL(double* Hess, void* PObjData);
void DDFuncSqValu_Obj5_PL(double *val, void* PObjData);
void DDFuncSqValuIndep_Obj5_PL(const double* P, double *val, void* PObjData);

void DDFuncSqGrad_Obj6_HN(const double* P, double* Grad, void* PObjData);
void DDFuncSqHess_Obj6_HN(double* Hess, void* PObjData);
void DDFuncSqValu_Obj6_HN(double *val, void* PObjData);
void DDFuncSqValuIndep_Obj6_HN(const double* P, double *val, void* PObjData);

struct ObjectDDFuncInfo {	
	void* PObjData;
	void(*DDFuncSqGrad)(const double*, double*, void*);
	void(*DDFuncSqHess)(double*, void*);
	void(*DDFuncSqValu)(double*, void*);
	void(*DDFuncSqValuIndep)(const double*, double*, void*);
};

inline double Sign(double y);




struct ObjElData {
	double a = 1.5, b = 1.2, c = 0.7;
	double x1, y1, z1;
	double F0, F1, F2, F3, F;
	bool F_activated;
};

struct Obj1MSData {
	double a[3] = { 0.7,1,1.5 }, b[3] = { 1.5,1.1,1.2 }, c[3] = { 1.2,1.2,0.6 };
	double x[3], y[3], z[3];
	double F0[3], F1[3], F2[3], F3[3];
	double F[4];
	bool F_activated[4];
};

struct Obj2PCData {
	double h = 2, r0 = 0.6;
	double l = sqrt(h * h + r0 * r0);
	double t0 = r0 / h;
	double c0 = h / l;
	double s0 = r0 / l;
	double a1 = 0.2, b1 = 0.4;
	double h2 = 1;

	double F[5];
	double r_sq, r, r_cu;
	double b_c2;
	double c0P0, c0P1, p2mh;

	bool Fco1, Fco2;
	double x1[3], y1[3];
	double xp[2], yp[2];

	double F10, F11, F12, F13;
	double F20, F21, F22, F23;
	double F30, F31, F32, F33;
	bool F_activated[5];
};

struct Obj3DMData {
	double a = 0.6, b = 1.3, c = 1.8;
	double n1 = 1.5, n2 = 1.6, n3 = 1.5;
	double x1, y1, z1;
	double sx, sy, sz;
	double x11, y11, z11;
	double x13, y13, z13;
	double F;
	bool F_activated;
};

struct Obj4CKData {
	double a = 2.1, b = 2.1, c = 0.5;
	double s = 4;
	double theta = PI / 3;
	double s0 = sin(theta), c0 = cos(theta), t0 = tan(theta);
	double x1, y1, z1, z11, z12, z13;
	double Fxy0, Fxy1, Fxy2, Fxy3;
	double Ftr0, Ftr1, Ftr2, Ftr3, Ftr;
	bool Ftr_activated;
	double F[3];
	bool F_activated[3];
};

struct Obj5PLData {
	double a = 0.7, b = 1.2, c = 0.4;
	double s = 6;
	double x1, y1, z1;
	double x11, x12, x13, y11, y12, y13;
	double Fxy0, Fxy1, Fxy2, Fxy3, Fpl0, Fpl1, Fpl2, Fpl;
	bool Fpl_activated;
};

struct Obj6HNData {
	double c = 0.4;
	double s = 16;
	double p1, q1, y1, z1;
	double p11, q11, y11, z11;
	double p13, q13, y13, z13;
	double sgp, sgq, sgy, sgz;
	double w1, w2, w3;
	double Fhn0, Fhn1, Fhn2, Fhn3, Fhn;
	bool Fhn_activated;
};

