#pragma once

// Value of V_{\Delta}(x)
// Gradient of V_{\Delta}(x)
// Hessian matrix of V_{\Delta}(x)

void DDFuncSqValu_Obj1_MS(const double* P, double *val);
void DDFuncSqGrad_Obj1_MS(const double* P, double* Grad);
void DDFuncSqHess_Obj1_MS(const double* P, double* Hess);

void DDFuncSqValu_Obj2_PC(const double* P, double *val);
void DDFuncSqGrad_Obj2_PC(const double* P, double* Grad);
void DDFuncSqHess_Obj2_PC(const double* P, double* Hess);

void DDFuncSqValu_Obj3_DM(const double* P, double *val);
void DDFuncSqGrad_Obj3_DM(const double* P, double* Grad);
void DDFuncSqHess_Obj3_DM(const double* P, double* Hess);

void DDFuncSqValu_Obj4_CK(const double* P, double *val);
void DDFuncSqGrad_Obj4_CK(const double* P, double* Grad);
void DDFuncSqHess_Obj4_CK(const double* P, double* Hess);

void DDFuncSqValu_Obj5_PL(const double* P, double *val);
void DDFuncSqGrad_Obj5_PL(const double* P, double* Grad);
void DDFuncSqHess_Obj5_PL(const double* P, double* Hess);

void DDFuncSqValu_Obj6_HN(const double* P, double *val);
void DDFuncSqGrad_Obj6_HN(const double* P, double* Grad);
void DDFuncSqHess_Obj6_HN(const double* P, double* Hess);

struct ObjectDDFuncInfo {
	void(*DDFuncSqValu)(const double*, double*);
	void(*DDFuncSqGrad)(const double*, double*);
	void(*DDFuncSqHess)(const double*, double*);
};

ObjectDDFuncInfo Object1DDFuncInfo = { DDFuncSqValu_Obj1_MS ,DDFuncSqGrad_Obj1_MS ,DDFuncSqHess_Obj1_MS };
ObjectDDFuncInfo Object2DDFuncInfo = { DDFuncSqValu_Obj2_PC ,DDFuncSqGrad_Obj2_PC ,DDFuncSqHess_Obj2_PC };
ObjectDDFuncInfo Object3DDFuncInfo = { DDFuncSqValu_Obj3_DM ,DDFuncSqGrad_Obj3_DM ,DDFuncSqHess_Obj3_DM };
ObjectDDFuncInfo Object4DDFuncInfo = { DDFuncSqValu_Obj4_CK ,DDFuncSqGrad_Obj4_CK ,DDFuncSqHess_Obj4_CK };
ObjectDDFuncInfo Object5DDFuncInfo = { DDFuncSqValu_Obj5_PL ,DDFuncSqGrad_Obj5_PL ,DDFuncSqHess_Obj5_PL };
ObjectDDFuncInfo Object6DDFuncInfo = { DDFuncSqValu_Obj6_HN ,DDFuncSqGrad_Obj6_HN ,DDFuncSqHess_Obj6_HN };

inline double Sign(double y);