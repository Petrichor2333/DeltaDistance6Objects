#include "pch.h"
#include "config.h"
#include "DDFunc.h"

/*-------------------------------   Object ellipse  ------------------------------------------*/

void DDFuncSqGrad_Obj_EL(const double* P, double* Grad, void* PObjData) {
	ObjElData* PObjElData = (ObjElData*)PObjData;
	PObjElData->x1 = P[0] / PObjElData->a, PObjElData->y1 = P[1] / PObjElData->b, PObjElData->z1 = P[2] / PObjElData->c;
	PObjElData->F0 = PObjElData->x1 * PObjElData->x1 + PObjElData->y1 * PObjElData->y1 + PObjElData->z1 * PObjElData->z1;
	PObjElData->F1 = sqrt(PObjElData->F0);
	PObjElData->F = PObjElData->F1 - 1;
	PObjElData->F_activated = PObjElData->F > 0 ? 1 : 0;
	if (PObjElData->F_activated) {
		Grad[0] = PObjElData->x1 * PObjElData->F * 2 / (PObjElData->a*PObjElData->F1);
		Grad[1] = PObjElData->y1 * PObjElData->F * 2 / (PObjElData->b*PObjElData->F1);
		Grad[2] = PObjElData->z1 * PObjElData->F * 2 / (PObjElData->c*PObjElData->F1);
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}
}

void DDFuncSqHess_Obj_EL(double* Hess, void* PObjData) {
	ObjElData* PObjElData = (ObjElData*)PObjData;
	PObjElData->F2 = PObjElData->F0 * PObjElData->F1;
	PObjElData->F3 = 1 / PObjElData->F0 - PObjElData->F / PObjElData->F2;
	if (PObjElData->F_activated) {
		Hess[0] = (PObjElData->x1*PObjElData->x1*PObjElData->F3 + PObjElData->F / PObjElData->F1) * 2 / (PObjElData->a*PObjElData->a);
		Hess[1] = 2 * PObjElData->x1*PObjElData->y1*PObjElData->F3 / (PObjElData->a*PObjElData->b);
		Hess[2] = 2 * PObjElData->x1*PObjElData->z1*PObjElData->F3 / (PObjElData->a*PObjElData->c);
		Hess[3] = (PObjElData->y1*PObjElData->y1*PObjElData->F3 + PObjElData->F / PObjElData->F1) * 2 / (PObjElData->b*PObjElData->b);
		Hess[4] = 2 * PObjElData->y1*PObjElData->z1*PObjElData->F3 / (PObjElData->b*PObjElData->c);
		Hess[5] = (PObjElData->z1*PObjElData->z1*PObjElData->F3 + PObjElData->F / PObjElData->F1) * 2 / (PObjElData->c*PObjElData->c);
	}
	else {
		Hess[0] = 0;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 0;
		Hess[4] = 0;
		Hess[5] = 0;
	}
}

void DDFuncSqValu_Obj_EL(double *val, void* PObjData) {
	ObjElData* PObjElData = (ObjElData*)PObjData;
	if (PObjElData->F_activated) {
		*val = PObjElData->F*PObjElData->F;
	}
	else {
		*val = 0;
	}

}

void DDFuncSqValuIndep_Obj_EL(const double* P, double *val, void* PObjData) {
	ObjElData* PObjElData = (ObjElData*)PObjData;
	double x1, y1, z1;
	x1 = P[0] / PObjElData->a, y1 = P[1] / PObjElData->b, z1 = P[2] / PObjElData->c;
	double F = sqrt(x1 * x1 + y1 * y1 + z1 * z1) - 1;
	if (F > 0) {
		*val = F * F;
	}
	else {
		*val = 0;
	}
}

/*-------------------------------   Object ellipse  ------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 1 Mouse  ------------------------------------------*/

void DDFuncSqGrad_Obj1_MS(const double* P, double* Grad, void* PObjData) {
	Obj1MSData* PObj1MSData = (Obj1MSData*)PObjData;
	for (int i = 0;i < 3;i++) {
		Grad[i] = 0;
	}

	for (int i = 0;i < 3;i++) {
		PObj1MSData->x[i] = P[0] / PObj1MSData->a[i];
		PObj1MSData->y[i] = P[1] / PObj1MSData->b[i];
		PObj1MSData->z[i] = P[2] / PObj1MSData->c[i];
		PObj1MSData->F0[i] = PObj1MSData->x[i] * PObj1MSData->x[i] + PObj1MSData->y[i] * PObj1MSData->y[i] + PObj1MSData->z[i] * PObj1MSData->z[i];
		PObj1MSData->F1[i] = sqrt(PObj1MSData->F0[i]);
		PObj1MSData->F[i] = PObj1MSData->F1[i] - 1;
		PObj1MSData->F_activated[i] = PObj1MSData->F[i] > 0 ? 1 : 0;
		if (PObj1MSData->F_activated[i]) {
			Grad[0] += PObj1MSData->x[i] * PObj1MSData->F[i] * 2 / (PObj1MSData->a[i] * PObj1MSData->F1[i]);
			Grad[1] += PObj1MSData->y[i] * PObj1MSData->F[i] * 2 / (PObj1MSData->b[i] * PObj1MSData->F1[i]);
			Grad[2] += PObj1MSData->z[i] * PObj1MSData->F[i] * 2 / (PObj1MSData->c[i] * PObj1MSData->F1[i]);
		}
	}

	PObj1MSData->F[3] = -P[2];
	PObj1MSData->F_activated[3] = PObj1MSData->F[3] > 0 ? 1 : 0;
	if (PObj1MSData->F_activated[3]) {
		Grad[2] -= 2 * PObj1MSData->F[3];
	}
}

void DDFuncSqHess_Obj1_MS(double* Hess, void* PObjData) {
	Obj1MSData* PObj1MSData = (Obj1MSData*)PObjData;
	for (int i = 0;i < 6;i++) {
		Hess[i] = 0;
	}
	
	for (int i = 0;i < 3;i++) {
		PObj1MSData->F2[i] = PObj1MSData->F0[i] * PObj1MSData->F1[i];
		PObj1MSData->F3[i] = 1 / PObj1MSData->F0[i] - PObj1MSData->F[i] / PObj1MSData->F2[i];
		if (PObj1MSData->F_activated[i]) {
			Hess[0] += (PObj1MSData->x[i] * PObj1MSData->x[i] * PObj1MSData->F3[i] + PObj1MSData->F[i] / PObj1MSData->F1[i]) * 2 / (PObj1MSData->a[i] * PObj1MSData->a[i]);
			Hess[1] += 2 * PObj1MSData->x[i] * PObj1MSData->y[i] * PObj1MSData->F3[i] / (PObj1MSData->a[i] * PObj1MSData->b[i]);
			Hess[2] += 2 * PObj1MSData->x[i] * PObj1MSData->z[i] * PObj1MSData->F3[i] / (PObj1MSData->a[i] * PObj1MSData->c[i]);
			Hess[3] += (PObj1MSData->y[i] * PObj1MSData->y[i] * PObj1MSData->F3[i] + PObj1MSData->F[i] / PObj1MSData->F1[i]) * 2 / (PObj1MSData->b[i] * PObj1MSData->b[i]);
			Hess[4] += 2 * PObj1MSData->y[i] * PObj1MSData->z[i] * PObj1MSData->F3[i] / (PObj1MSData->b[i] * PObj1MSData->c[i]);
			Hess[5] += (PObj1MSData->z[i] * PObj1MSData->z[i] * PObj1MSData->F3[i] + PObj1MSData->F[i] / PObj1MSData->F1[i]) * 2 / (PObj1MSData->c[i] * PObj1MSData->c[i]);
		}
	}

	if (PObj1MSData->F_activated[3]) {
		Hess[5] += 2;
	}
}

void DDFuncSqValu_Obj1_MS(double *val, void* PObjData) {
	Obj1MSData* PObj1MSData = (Obj1MSData*)PObjData;

	*val = 0;
	for (int i = 0;i < 4;i++) {
		if (PObj1MSData->F_activated[i]) {
			*val += PObj1MSData->F[i] * PObj1MSData->F[i];
		}
	}
}

void DDFuncSqValuIndep_Obj1_MS(const double* P, double *val, void* PObjData) {
	Obj1MSData* PObj1MSData = (Obj1MSData*)PObjData;
	
	*val = 0;
	double x[3], y[3], z[3];
	double F[4];

	for (int i = 0;i < 3;i++) {
		x[i] = P[0] / PObj1MSData->a[i];
		y[i] = P[1] / PObj1MSData->b[i];
		z[i] = P[2] / PObj1MSData->c[i];
		F[i] = sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]) - 1;
		if (F[i] > 0) {
			*val += F[i] * F[i];
		}
	}

	F[3] = -P[2];
	if (F[3] > 0) {
		*val += F[3] * F[3];
	}

}

/*-------------------------------   Object 1 Mouse  --------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 2 Pencil  ------------------------------------------*/



void DDFuncSqGrad_Obj2_PC(const double* P, double* Grad, void* PObjData) {
	Obj2PCData* PObj2PCData = (Obj2PCData*)PObjData;
	double sqrt3 = sqrt(3);

	PObj2PCData->r_sq = P[0] * P[0] + P[1] * P[1];
	PObj2PCData->r = sqrt(PObj2PCData->r_sq);
	
	PObj2PCData->p2mh = P[2] - PObj2PCData->h;
	PObj2PCData->Fco1 = PObj2PCData->p2mh > PObj2PCData->t0 * PObj2PCData->r ? 1 : 0;
	PObj2PCData->Fco2 = PObj2PCData->p2mh > -PObj2PCData->r / PObj2PCData->t0 ? 1 : 0;

	if (PObj2PCData->Fco1) {
		Grad[0] = 2 * P[0];
		Grad[1] = 2 * P[1];
		Grad[2] = 2 * PObj2PCData->p2mh;
	}
	else if (PObj2PCData->Fco2) {
		PObj2PCData->b_c2 = PObj2PCData->c0 * PObj2PCData->r + PObj2PCData->s0 * P[2] - PObj2PCData->c0 * PObj2PCData->r0;
		PObj2PCData->c0P0 = PObj2PCData->c0 * P[0], PObj2PCData->c0P1 = PObj2PCData->c0 * P[1];
		Grad[0] = 2 * PObj2PCData->b_c2*PObj2PCData->c0P0 / PObj2PCData->r;
		Grad[1] = 2 * PObj2PCData->b_c2*PObj2PCData->c0P1 / PObj2PCData->r;
		Grad[2] = 2 * PObj2PCData->b_c2*PObj2PCData->s0;
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}

	PObj2PCData->x1[0] = P[0] / PObj2PCData->a1, PObj2PCData->y1[0] = P[1] / PObj2PCData->b1;
	PObj2PCData->F10 = PObj2PCData->x1[0] * PObj2PCData->x1[0] + PObj2PCData->y1[0] * PObj2PCData->y1[0];
	PObj2PCData->F11 = sqrt(PObj2PCData->F10);
	PObj2PCData->F[1] = PObj2PCData->F11 - 1;
	PObj2PCData->F_activated[1] = PObj2PCData->F[1] > 0 ? 1 : 0;
	if (PObj2PCData->F_activated[1]) {
		Grad[0] += PObj2PCData->x1[0] * PObj2PCData->F[1] * 2 / (PObj2PCData->a1*PObj2PCData->F11);
		Grad[1] += PObj2PCData->y1[0] * PObj2PCData->F[1] * 2 / (PObj2PCData->b1*PObj2PCData->F11);
	}
	
	PObj2PCData->xp[0] = P[0] / 2 + sqrt3 * P[1] / 2, PObj2PCData->yp[0] = -sqrt3 * P[0] / 2 + P[1] / 2;
	PObj2PCData->x1[1] = PObj2PCData->xp[0] / PObj2PCData->a1, PObj2PCData->y1[1] = PObj2PCData->yp[0] / PObj2PCData->b1;
	PObj2PCData->F20 = PObj2PCData->x1[1] * PObj2PCData->x1[1] + PObj2PCData->y1[1] * PObj2PCData->y1[1];
	PObj2PCData->F21 = sqrt(PObj2PCData->F20);
	PObj2PCData->F[2] = PObj2PCData->F21 - 1;
	PObj2PCData->F_activated[2] = PObj2PCData->F[2] > 0 ? 1 : 0;
	if (PObj2PCData->F_activated[2]) {
		double Gradp[2];
		Gradp[0] = PObj2PCData->x1[1] * PObj2PCData->F[2] * 2 / (PObj2PCData->a1*PObj2PCData->F21);
		Gradp[1] = PObj2PCData->y1[1] * PObj2PCData->F[2] * 2 / (PObj2PCData->b1*PObj2PCData->F21);
		Grad[0] += Gradp[0] / 2 - sqrt3 * Gradp[1] / 2;
		Grad[1] += sqrt3 * Gradp[0] / 2 + Gradp[1] / 2;
	}

	PObj2PCData->xp[1] = P[0] / 2 - sqrt3 * P[1] / 2, PObj2PCData->yp[1] = sqrt3 * P[0] / 2 + P[1] / 2;
	PObj2PCData->x1[2] = PObj2PCData->xp[1] / PObj2PCData->a1, PObj2PCData->y1[2] = PObj2PCData->yp[1] / PObj2PCData->b1;
	PObj2PCData->F30 = PObj2PCData->x1[2] * PObj2PCData->x1[2] + PObj2PCData->y1[2] * PObj2PCData->y1[2];
	PObj2PCData->F31 = sqrt(PObj2PCData->F30);
	PObj2PCData->F[3] = PObj2PCData->F31 - 1;
	PObj2PCData->F_activated[3] = PObj2PCData->F[3] > 0 ? 1 : 0;
	if (PObj2PCData->F_activated[3]) {
		double Gradp[2];
		Gradp[0] = PObj2PCData->x1[2] * PObj2PCData->F[3] * 2 / (PObj2PCData->a1*PObj2PCData->F31);
		Gradp[1] = PObj2PCData->y1[2] * PObj2PCData->F[3] * 2 / (PObj2PCData->b1*PObj2PCData->F31);
		Grad[0] += Gradp[0] / 2 + sqrt3 * Gradp[1] / 2;
		Grad[1] += -sqrt3 * Gradp[0] / 2 + Gradp[1] / 2;
	}
	
	PObj2PCData->F[4] = -P[2] - PObj2PCData->h2;
	PObj2PCData->F_activated[4] = PObj2PCData->F[4] > 0 ? 1 : 0;
	if (PObj2PCData->F_activated[4]) {
		Grad[2] -= 2 * PObj2PCData->F[4];
	}

}

void DDFuncSqHess_Obj2_PC(double* Hess, void* PObjData) {
	Obj2PCData* PObj2PCData = (Obj2PCData*)PObjData;
	double sqrt3 = sqrt(3);

	if (PObj2PCData->Fco1) {
		Hess[0] = 2;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 2;
		Hess[4] = 0;
		Hess[5] = 2;
	}
	else if (PObj2PCData->Fco2) {
		double temp = (1 - PObj2PCData->b_c2 / (PObj2PCData->c0 * PObj2PCData->r)) / PObj2PCData->r_sq;
		Hess[0] = 2 * (PObj2PCData->c0P0*PObj2PCData->c0P0 * temp + PObj2PCData->c0 * PObj2PCData->b_c2 / PObj2PCData->r);
		Hess[1] = 2 * PObj2PCData->c0P0*PObj2PCData->c0P1 * temp;
		Hess[2] = 2 * PObj2PCData->c0P0*PObj2PCData->s0 / PObj2PCData->r;
		Hess[3] = 2 * (PObj2PCData->c0P1*PObj2PCData->c0P1 * temp + PObj2PCData->c0 * PObj2PCData->b_c2 / PObj2PCData->r);
		Hess[4] = 2 * PObj2PCData->c0P1*PObj2PCData->s0 / PObj2PCData->r;
		Hess[5] = 2 * PObj2PCData->s0*PObj2PCData->s0;
	}
	else {
		Hess[0] = 0;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 0;
		Hess[4] = 0;
		Hess[5] = 0;
	}

	double a1sq = PObj2PCData->a1*PObj2PCData->a1;
	double a1b1 = PObj2PCData->a1*PObj2PCData->b1;
	double b1sq = PObj2PCData->b1*PObj2PCData->b1;

	PObj2PCData->F12 = PObj2PCData->F10 * PObj2PCData->F11;
	PObj2PCData->F13 = 1 / PObj2PCData->F10 - PObj2PCData->F[1] / PObj2PCData->F12;
	if (PObj2PCData->F_activated[1]) {
		Hess[0] += (PObj2PCData->x1[0] * PObj2PCData->x1[0] * PObj2PCData->F13 + PObj2PCData->F[1] / PObj2PCData->F11) * 2 / a1sq;
		Hess[1] += 2 * PObj2PCData->x1[0] * PObj2PCData->y1[0] * PObj2PCData->F13 / a1b1;
		Hess[3] += (PObj2PCData->y1[0] * PObj2PCData->y1[0] * PObj2PCData->F13 + PObj2PCData->F[1] / PObj2PCData->F11) * 2 / b1sq;
	}
	
	PObj2PCData->F22 = PObj2PCData->F20 * PObj2PCData->F21;
	PObj2PCData->F23 = 1 / PObj2PCData->F20 - PObj2PCData->F[2] / PObj2PCData->F22;
	if (PObj2PCData->F_activated[2]) {
		double Hessp[3];
		Hessp[0] = (PObj2PCData->x1[1]* PObj2PCData->x1[1]* PObj2PCData->F23 + PObj2PCData->F[2] / PObj2PCData->F21) * 2 / a1sq;
		Hessp[1] = 2 * PObj2PCData->x1[1]* PObj2PCData->y1[1]* PObj2PCData->F23 / a1b1;
		Hessp[2] = (PObj2PCData->y1[1]* PObj2PCData->y1[1]* PObj2PCData->F23 + PObj2PCData->F[2] / PObj2PCData->F21) * 2 / b1sq;

		Hess[0] += Hessp[0] / 4 - sqrt3 * Hessp[1] / 2 + 3 * Hessp[2] / 4;
		Hess[1] += sqrt3 * Hessp[0] / 4 - Hessp[1] / 2 - sqrt3 * Hessp[2] / 4;
		Hess[3] += 3 * Hessp[0] / 4 + sqrt3 * Hessp[1] / 2 + Hessp[2] / 4;
	}

	PObj2PCData->F32 = PObj2PCData->F30 * PObj2PCData->F31;
	PObj2PCData->F33 = 1 / PObj2PCData->F30 - PObj2PCData->F[3] / PObj2PCData->F32;
	if (PObj2PCData->F_activated[3]) {
		double Hessp[3];
		Hessp[0] = (PObj2PCData->x1[2]* PObj2PCData->x1[2]* PObj2PCData->F33 + PObj2PCData->F[3] / PObj2PCData->F31) * 2 / a1sq;
		Hessp[1] = 2 * PObj2PCData->x1[2]* PObj2PCData->y1[2]* PObj2PCData->F33 / a1b1;
		Hessp[2] = (PObj2PCData->y1[2]* PObj2PCData->y1[2]* PObj2PCData->F33 + PObj2PCData->F[3] / PObj2PCData->F31) * 2 / b1sq;

		Hess[0] += Hessp[0] / 4 + sqrt3 * Hessp[1] / 2 + 3 * Hessp[2] / 4;
		Hess[1] += -sqrt3 * Hessp[0] / 4 - Hessp[1] / 2 + sqrt3 * Hessp[2] / 4;
		Hess[3] += 3 * Hessp[0] / 4 - sqrt3 * Hessp[1] / 2 + Hessp[2] / 4;
	}
	
	if (PObj2PCData->F_activated[4]) {
		Hess[5] += 2;
	}
}

void DDFuncSqValu_Obj2_PC(double *val, void* PObjData) {
	Obj2PCData* PObj2PCData = (Obj2PCData*)PObjData;
	double sqrt3 = sqrt(3);

	if (PObj2PCData->Fco1) {
		*val = PObj2PCData->r_sq + PObj2PCData->p2mh*PObj2PCData->p2mh;
	}
	else if (PObj2PCData->Fco2) {
		*val = PObj2PCData->b_c2 * PObj2PCData->b_c2;
	}
	else {
		*val = 0;
	}

	for (int i = 1;i <= 4;i++) {
		if (PObj2PCData->F_activated[i]) {
			*val += PObj2PCData->F[i] * PObj2PCData->F[i];
		}
	}
}

void DDFuncSqValuIndep_Obj2_PC(const double* P, double *val, void* PObjData) {
	Obj2PCData* PObj2PCData = (Obj2PCData*)PObjData;
	double sqrt3 = sqrt(3);

	double r_sq = P[0] * P[0] + P[1] * P[1];
	double r = sqrt(r_sq);
	double F[5];


	if (P[2] - PObj2PCData->t0 * r > PObj2PCData->h) {
		*val = r_sq + (P[2] - PObj2PCData->h)*(P[2] - PObj2PCData->h);
	}
	else if (P[2] + r / PObj2PCData->t0 > PObj2PCData->h) {
		double b_c2 = PObj2PCData->c0 * r + PObj2PCData->s0 * P[2] - PObj2PCData->c0 * PObj2PCData->r0;
		*val = b_c2 * b_c2;
	}
	else {
		*val = 0;
	}

	double x1, y1;
	double xp, yp;

	x1 = P[0] / PObj2PCData->a1, y1 = P[1] / PObj2PCData->b1;
	F[1] = sqrt(x1 * x1 + y1 * y1) - 1;

	xp = P[0] / 2 + sqrt3 * P[1] / 2, yp = -sqrt3 * P[0] / 2 + P[1] / 2;
	x1 = xp / PObj2PCData->a1, y1 = yp / PObj2PCData->b1;
	F[2] = sqrt(x1 * x1 + y1 * y1) - 1;

	xp = P[0] / 2 - sqrt3 * P[1] / 2, yp = sqrt3 * P[0] / 2 + P[1] / 2;
	x1 = xp / PObj2PCData->a1, y1 = yp / PObj2PCData->b1;
	F[3] = sqrt(x1 * x1 + y1 * y1) - 1;


	F[4] = -P[2] - PObj2PCData->h2;

	for (int i = 1;i <= 4;i++) {
		if (F[i] > 0) {
			*val += F[i] * F[i];
		}
	}
}

/*-------------------------------  Object 2 Pencil --------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 3 Diamond  ------------------------------------------*/



void DDFuncSqGrad_Obj3_DM(const double* P, double* Grad, void* PObjData) {
	Obj3DMData* pObj3DMData = (Obj3DMData*)PObjData;
	pObj3DMData->x1 = fabs(P[0] / pObj3DMData->a), pObj3DMData->y1 = fabs(P[1] / pObj3DMData->b), pObj3DMData->z1 = fabs(P[2] / pObj3DMData->c);
	pObj3DMData->sx = Sign(P[0]), pObj3DMData->sy = Sign(P[1]), pObj3DMData->sz = Sign(P[2]);

	pObj3DMData->x13 = pow(pObj3DMData->x1, pObj3DMData->n1 - 2), pObj3DMData->y13 = pow(pObj3DMData->y1, pObj3DMData->n2 - 2), pObj3DMData->z13 = pow(pObj3DMData->z1, pObj3DMData->n3 - 2);
	pObj3DMData->x11 = pow(pObj3DMData->x1, pObj3DMData->n1 - 1), pObj3DMData->y11 = pow(pObj3DMData->y1, pObj3DMData->n2 - 1), pObj3DMData->z11 = pow(pObj3DMData->z1, pObj3DMData->n3 - 1);

	pObj3DMData->F = pObj3DMData->x11 * pObj3DMData->x1 + pObj3DMData->y11 * pObj3DMData->y1 + pObj3DMData->z11 * pObj3DMData->z1 - 1;
	pObj3DMData->F_activated = pObj3DMData->F > 0 ? 1 : 0;

	if (pObj3DMData->F_activated) {
		Grad[0] = 2 * pObj3DMData->sx * pObj3DMData->n1*pObj3DMData->F*pObj3DMData->x11 / pObj3DMData->a;
		Grad[1] = 2 * pObj3DMData->sy * pObj3DMData->n2*pObj3DMData->F*pObj3DMData->y11 / pObj3DMData->b;
		Grad[2] = 2 * pObj3DMData->sz * pObj3DMData->n3*pObj3DMData->F*pObj3DMData->z11 / pObj3DMData->c;
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}
}

void DDFuncSqHess_Obj3_DM(double* Hess, void* PObjData) {
	Obj3DMData* pObj3DMData = (Obj3DMData*)PObjData;

	if (pObj3DMData->F_activated) {
		Hess[0] = 2 * pObj3DMData->n1*(pObj3DMData->n1*pObj3DMData->x11*pObj3DMData->x11 + (pObj3DMData->n1 - 1)*pObj3DMData->F*pObj3DMData->x13) / (pObj3DMData->a*pObj3DMData->a);
		Hess[1] = 2 * pObj3DMData->sx*pObj3DMData->sy * pObj3DMData->n1*pObj3DMData->n2*pObj3DMData->x11*pObj3DMData->y11 / (pObj3DMData->a*pObj3DMData->b);
		Hess[2] = 2 * pObj3DMData->sx*pObj3DMData->sz * pObj3DMData->n1*pObj3DMData->n3*pObj3DMData->x11*pObj3DMData->z11 / (pObj3DMData->a*pObj3DMData->c);
		Hess[3] = 2 * pObj3DMData->n2*(pObj3DMData->n2*pObj3DMData->y11*pObj3DMData->y11 + (pObj3DMData->n2 - 1)*pObj3DMData->F*pObj3DMData->y13) / (pObj3DMData->b*pObj3DMData->b);
		Hess[4] = 2 * pObj3DMData->sy*pObj3DMData->sz * pObj3DMData->n2*pObj3DMData->n3*pObj3DMData->y11*pObj3DMData->z11 / (pObj3DMData->b*pObj3DMData->c);
		Hess[5] = 2 * pObj3DMData->n3*(pObj3DMData->n3*pObj3DMData->z11*pObj3DMData->z11 + (pObj3DMData->n3 - 1)*pObj3DMData->F*pObj3DMData->z13) / (pObj3DMData->c*pObj3DMData->c);
	}
	else {
		Hess[0] = 0;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 0;
		Hess[4] = 0;
		Hess[5] = 0;
	}
}

void DDFuncSqValu_Obj3_DM(double *val, void* PObjData) {
	Obj3DMData* pObj3DMData = (Obj3DMData*)PObjData;
	if (pObj3DMData->F_activated) {
		*val = pObj3DMData->F * pObj3DMData->F;
	}
	else {
		*val = 0;
	}
}

void DDFuncSqValuIndep_Obj3_DM(const double* P, double *val, void* PObjData) {
	Obj3DMData* pObj3DMData = (Obj3DMData*)PObjData;

	double x1 = fabs(P[0] / pObj3DMData->a), y1 = fabs(P[1] / pObj3DMData->b), z1 = fabs(P[2] / pObj3DMData->c);
	double F = pow(x1, pObj3DMData->n1) + pow(y1, pObj3DMData->n2) + pow(z1, pObj3DMData->n3) - 1;
	if (F > 0) {
		*val = F * F;
	}
	else {
		*val = 0;
	}
}
/*-------------------------------   Object 3 Diamond  ------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 4 Piece of cake  ------------------------------------------*/



void DDFuncSqGrad_Obj4_CK(const double* P, double* Grad, void* PObjData) {
	Obj4CKData* PObj4CKData = (Obj4CKData*)PObjData;
	
	PObj4CKData->x1 = P[0] / PObj4CKData->a, PObj4CKData->y1 = P[1] / PObj4CKData->b, PObj4CKData->z1 = P[2] / PObj4CKData->c;
	
	PObj4CKData->Fxy0 = PObj4CKData->x1 * PObj4CKData->x1 + PObj4CKData->y1 * PObj4CKData->y1;
	PObj4CKData->Fxy3 = pow(PObj4CKData->Fxy0, PObj4CKData->s - 2);
	PObj4CKData->Fxy1 = PObj4CKData->Fxy3 * PObj4CKData->Fxy0;
	
	PObj4CKData->z13 = pow(PObj4CKData->z1*PObj4CKData->z1, PObj4CKData->s - 1);
	PObj4CKData->z11 = PObj4CKData->z13 * PObj4CKData->z1;
	PObj4CKData->Ftr0 = PObj4CKData->Fxy1 * PObj4CKData->Fxy0 + PObj4CKData->z11 * PObj4CKData->z1;
	PObj4CKData->Ftr3 = pow(PObj4CKData->Ftr0, 0.5 / PObj4CKData->s - 2);
	PObj4CKData->Ftr1 = PObj4CKData->Ftr3 * PObj4CKData->Ftr0;
	PObj4CKData->Ftr = PObj4CKData->Ftr0 * PObj4CKData->Ftr1 - 1;
	PObj4CKData->Ftr_activated = PObj4CKData->Ftr > 0 ? 1 : 0;

	if (PObj4CKData->Ftr_activated) {
		Grad[0] = 2 * PObj4CKData->x1 * PObj4CKData->Fxy1*PObj4CKData->Ftr1*PObj4CKData->Ftr / PObj4CKData->a;
		Grad[1] = 2 * PObj4CKData->y1 * PObj4CKData->Fxy1*PObj4CKData->Ftr1*PObj4CKData->Ftr / PObj4CKData->b;
		Grad[2] = 2 * PObj4CKData->z11 * PObj4CKData->Ftr1*PObj4CKData->Ftr / PObj4CKData->c;
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}

	PObj4CKData->F[0] = -P[2];
	PObj4CKData->F[1] = -P[0];
	PObj4CKData->F[2] = PObj4CKData->c0 * P[0] - PObj4CKData->s0 * P[1];
	PObj4CKData->F_activated[0] = PObj4CKData->F[0] > 0 ? 1 : 0;
	PObj4CKData->F_activated[1] = PObj4CKData->F[1] > 0 ? 1 : 0;
	PObj4CKData->F_activated[2] = PObj4CKData->F[2] > 0 ? 1 : 0;

	if (PObj4CKData->F_activated[0]) {
		Grad[2] += 2 * P[2];
	}
	if (PObj4CKData->F_activated[1]) {
		Grad[0] += 2 * P[0];
	}
	if (PObj4CKData->F_activated[2]) {
		Grad[0] += 2 * PObj4CKData->F[2] * PObj4CKData->c0;
		Grad[1] -= 2 * PObj4CKData->F[2] * PObj4CKData->s0;
	}

}

void DDFuncSqHess_Obj4_CK(double* Hess, void* PObjData) {
	Obj4CKData* PObj4CKData = (Obj4CKData*)PObjData;

	PObj4CKData->Fxy2 = PObj4CKData->Fxy1 * PObj4CKData->Fxy1;
	
	PObj4CKData->z12 = PObj4CKData->z11 * PObj4CKData->z11;
	
	PObj4CKData->Ftr2 = PObj4CKData->Ftr1 * PObj4CKData->Ftr1;
	

	double s1 = -2 * PObj4CKData->s + 1;
	double s2 = 2 * PObj4CKData->s - 2;

	double T1 = PObj4CKData->Ftr2 + s1 * PObj4CKData->Ftr3*PObj4CKData->Ftr;
	double T2 = s2 * PObj4CKData->Fxy3*PObj4CKData->Ftr1*PObj4CKData->Ftr;

	if (PObj4CKData->Ftr_activated) {
		Hess[0] = (PObj4CKData->x1*PObj4CKData->x1*(PObj4CKData->Fxy2*T1 + T2) + PObj4CKData->Fxy1 * PObj4CKData->Ftr1*PObj4CKData->Ftr) * 2 / (PObj4CKData->a*PObj4CKData->a);
		Hess[1] = PObj4CKData->x1 * PObj4CKData->y1*(PObj4CKData->Fxy2*T1 + T2) * 2 / (PObj4CKData->a*PObj4CKData->b);
		Hess[2] = PObj4CKData->x1 * PObj4CKData->Fxy1*PObj4CKData->z11*T1 * 2 / (PObj4CKData->a*PObj4CKData->c);
		Hess[3] = (PObj4CKData->y1*PObj4CKData->y1*(PObj4CKData->Fxy2*T1 + T2) + PObj4CKData->Fxy1 * PObj4CKData->Ftr1*PObj4CKData->Ftr) * 2 / (PObj4CKData->b*PObj4CKData->b);
		Hess[4] = PObj4CKData->y1 * PObj4CKData->Fxy1*PObj4CKData->z11*T1 * 2 / (PObj4CKData->b*PObj4CKData->c);
		Hess[5] = (PObj4CKData->z12*T1 - s1 * PObj4CKData->z13*PObj4CKData->Ftr1*PObj4CKData->Ftr) * 2 / (PObj4CKData->c*PObj4CKData->c);
	}
	else {
		Hess[0] = 0;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 0;
		Hess[4] = 0;
		Hess[5] = 0;
	}

	if (PObj4CKData->F_activated[0]) {
		Hess[5] += 2;
	}
	if (PObj4CKData->F_activated[1]) {
		Hess[0] += 2;
	}
	if (PObj4CKData->F_activated[2]) {
		Hess[0] += 2 * PObj4CKData->c0 * PObj4CKData->c0;
		Hess[1] -= 2 * PObj4CKData->c0 * PObj4CKData->s0;
		Hess[3] += 2 * PObj4CKData->s0 * PObj4CKData->s0;
	}
}

void DDFuncSqValu_Obj4_CK(double *val, void* PObjData) {
	Obj4CKData* PObj4CKData = (Obj4CKData*)PObjData;
	if (PObj4CKData->Ftr_activated) {
		*val = PObj4CKData->Ftr * PObj4CKData->Ftr;
	}
	else {
		*val = 0;
	}

	for (int i = 0;i < 3;i++) {
		if (PObj4CKData->F_activated[i]) {
			*val += PObj4CKData->F[i] * PObj4CKData->F[i];
		}
	}

}

void DDFuncSqValuIndep_Obj4_CK(const double* P, double *val, void* PObjData) {
	Obj4CKData* PObj4CKData = (Obj4CKData*)PObjData;

	double x1 = P[0] / PObj4CKData->a, y1 = P[1] / PObj4CKData->b, z1 = P[2] / PObj4CKData->c;
	double Fxy0 = x1 * x1 + y1 * y1;

	double Fse = pow(pow(Fxy0, PObj4CKData->s) + pow(z1*z1, PObj4CKData->s), 0.5 / PObj4CKData->s) - 1;
	if (Fse > 0) {
		*val = Fse * Fse;
	}
	else {
		*val = 0;
	}

	double F[3];
	int k = 3;
	F[0] = -P[2];
	F[1] = -P[0];
	F[2] = PObj4CKData->c0 * P[0] - PObj4CKData->s0 * P[1];

	for (int i = 0;i < k;i++) {
		if (F[i] > 0) {
			*val += F[i] * F[i];
		}
	}

}

/*-------------------------------   Object 4 Piece of cake  ------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 5 Pillow  ------------------------------------------*/

void DDFuncSqGrad_Obj5_PL(const double* P, double* Grad, void* PObjData) {
	Obj5PLData* PObj5PLData = (Obj5PLData*)PObjData;
	
	PObj5PLData->x1 = P[0] / PObj5PLData->a, PObj5PLData->y1 = P[1] / PObj5PLData->b, PObj5PLData->z1 = P[2] / PObj5PLData->c;
	PObj5PLData->x13 = pow(PObj5PLData->x1*PObj5PLData->x1, PObj5PLData->s - 2);
	PObj5PLData->y13 = pow(PObj5PLData->y1*PObj5PLData->y1, PObj5PLData->s - 2);
	PObj5PLData->x11 = PObj5PLData->x13 * PObj5PLData->x1*PObj5PLData->x1;
	PObj5PLData->y11 = PObj5PLData->y13 * PObj5PLData->y1*PObj5PLData->y1;

	PObj5PLData->Fxy0 = pow(PObj5PLData->x1*PObj5PLData->x1, PObj5PLData->s) + pow(PObj5PLData->y1*PObj5PLData->y1, PObj5PLData->s);
	PObj5PLData->Fxy3 = pow(PObj5PLData->Fxy0, 1 / PObj5PLData->s - 2);
	PObj5PLData->Fxy1 = PObj5PLData->Fxy3 * PObj5PLData->Fxy0;

	PObj5PLData->Fpl0 = PObj5PLData->Fxy0 * PObj5PLData->Fxy1 + PObj5PLData->z1 * PObj5PLData->z1;
	PObj5PLData->Fpl1 = sqrt(PObj5PLData->Fpl0);
	PObj5PLData->Fpl = PObj5PLData->Fpl1 - 1;

	PObj5PLData->Fpl_activated = PObj5PLData->Fpl > 0 ? 1 : 0;

	if (PObj5PLData->Fpl_activated) {
		Grad[0] = 2 * PObj5PLData->x1*PObj5PLData->x11*PObj5PLData->Fxy1*PObj5PLData->Fpl / (PObj5PLData->a*PObj5PLData->Fpl1);
		Grad[1] = 2 * PObj5PLData->y1*PObj5PLData->y11*PObj5PLData->Fxy1*PObj5PLData->Fpl / (PObj5PLData->b*PObj5PLData->Fpl1);
		Grad[2] = 2 * PObj5PLData->z1*PObj5PLData->Fpl / (PObj5PLData->c*PObj5PLData->Fpl1);
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}

}

void DDFuncSqHess_Obj5_PL(double* Hess, void* PObjData) {
	Obj5PLData* PObj5PLData = (Obj5PLData*)PObjData;

	PObj5PLData->x12 = PObj5PLData->x11 * PObj5PLData->x11, PObj5PLData->y12 = PObj5PLData->y11 * PObj5PLData->y11;

	PObj5PLData->Fxy2 = PObj5PLData->Fxy1 * PObj5PLData->Fxy1;

	PObj5PLData->Fpl2 = PObj5PLData->Fpl0 * PObj5PLData->Fpl1;

	double s1 = -2 * PObj5PLData->s + 2;
	double T0 = 1 / PObj5PLData->Fpl0 - PObj5PLData->Fpl / PObj5PLData->Fpl2;
	double T1 = PObj5PLData->Fxy2 * T0 + s1 * PObj5PLData->Fxy3 * PObj5PLData->Fpl / PObj5PLData->Fpl1;

	if (PObj5PLData->Fpl_activated) {
		Hess[0] = (PObj5PLData->x1*PObj5PLData->x1*(PObj5PLData->x12*T1 - s1 * PObj5PLData->x13 * PObj5PLData->Fxy1*PObj5PLData->Fpl / PObj5PLData->Fpl1) + PObj5PLData->x11 * PObj5PLData->Fxy1*PObj5PLData->Fpl / PObj5PLData->Fpl1) * 2 / (PObj5PLData->a*PObj5PLData->a);
		Hess[1] = PObj5PLData->x1 * PObj5PLData->y1*PObj5PLData->x11*PObj5PLData->y11*T1 * 2 / (PObj5PLData->a*PObj5PLData->b);
		Hess[2] = PObj5PLData->x1 * PObj5PLData->z1*PObj5PLData->x11*PObj5PLData->Fxy1*T0 * 2 / (PObj5PLData->a*PObj5PLData->c);
		Hess[3] = (PObj5PLData->y1*PObj5PLData->y1*(PObj5PLData->y12*T1 - s1 * PObj5PLData->y13 * PObj5PLData->Fxy1*PObj5PLData->Fpl / PObj5PLData->Fpl1) + PObj5PLData->y11 * PObj5PLData->Fxy1*PObj5PLData->Fpl / PObj5PLData->Fpl1) * 2 / (PObj5PLData->b*PObj5PLData->b);
		Hess[4] = PObj5PLData->y1 * PObj5PLData->z1*PObj5PLData->y11*PObj5PLData->Fxy1*T0 * 2 / (PObj5PLData->b*PObj5PLData->c);
		Hess[5] = (PObj5PLData->z1*PObj5PLData->z1*T0 + PObj5PLData->Fpl / PObj5PLData->Fpl1) * 2 / (PObj5PLData->c*PObj5PLData->c);
	}
	else {
		Hess[0] = 0;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 0;
		Hess[4] = 0;
		Hess[5] = 0;
	}
}

void DDFuncSqValu_Obj5_PL(double *val, void* PObjData) {
	Obj5PLData* PObj5PLData = (Obj5PLData*)PObjData;

	if (PObj5PLData->Fpl_activated) {
		*val = PObj5PLData->Fpl * PObj5PLData->Fpl;
	}
	else {
		*val = 0;
	}
}

void DDFuncSqValuIndep_Obj5_PL(const double* P, double *val, void* PObjData) {
	Obj5PLData* PObj5PLData = (Obj5PLData*)PObjData;

	double x1 = P[0] / PObj5PLData->a, y1 = P[1] / PObj5PLData->b, z1 = P[2] / PObj5PLData->c;
	double Fxy0 = pow(x1*x1, PObj5PLData->s) + pow(y1*y1, PObj5PLData->s);
	double Fpl = sqrt(pow(Fxy0, 1 / PObj5PLData->s) + z1 * z1) - 1;

	if (Fpl > 0) {
		*val = Fpl * Fpl;
	}
	else {
		*val = 0;
	}
}
/*-------------------------------   Object 5 Pillow  ------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 6 Hexagonal nut  ------------------------------------------*/



void DDFuncSqGrad_Obj6_HN(const double* P, double* Grad, void* PObjData) {
	Obj6HNData* PObj6HNData = (Obj6HNData*)PObjData;

	PObj6HNData->p1 = fabs(P[0] - 0.5*P[1]), PObj6HNData->q1 = fabs(P[0] + 0.5*P[1]);
	PObj6HNData->y1 = fabs(P[1]), PObj6HNData->z1 = fabs(P[2] / PObj6HNData->c);
	
	PObj6HNData->p13 = pow(PObj6HNData->p1, PObj6HNData->s - 2), PObj6HNData->q13 = pow(PObj6HNData->q1, PObj6HNData->s - 2);
	PObj6HNData->y13 = pow(PObj6HNData->y1, PObj6HNData->s - 2), PObj6HNData->z13 = pow(PObj6HNData->z1, PObj6HNData->s - 2);
	PObj6HNData->p11 = PObj6HNData->p13 * PObj6HNData->p1, PObj6HNData->q11 = PObj6HNData->q13 * PObj6HNData->q1;
	PObj6HNData->y11 = PObj6HNData->y13 * PObj6HNData->y1, PObj6HNData->z11 = PObj6HNData->z13 * PObj6HNData->z1;

	PObj6HNData->sgp = Sign(P[0] - 0.5*P[1]), PObj6HNData->sgq = Sign(P[0] + 0.5*P[1]);
	PObj6HNData->sgy = Sign(P[1]), PObj6HNData->sgz = Sign(P[2]);
	
	PObj6HNData->Fhn0 = PObj6HNData->p11 * PObj6HNData->p1 + PObj6HNData->q11 * PObj6HNData->q1 + PObj6HNData->y11 * PObj6HNData->y1 + PObj6HNData->z11 * PObj6HNData->z1;
	PObj6HNData->Fhn3 = pow(PObj6HNData->Fhn0, -2 + 1 / PObj6HNData->s);
	PObj6HNData->Fhn1 = PObj6HNData->Fhn3 * PObj6HNData->Fhn0;
	PObj6HNData->Fhn = PObj6HNData->Fhn0 * PObj6HNData->Fhn1 - 1;
	PObj6HNData->Fhn_activated = PObj6HNData->Fhn > 0 ? 1 : 0;

	PObj6HNData->w1 = PObj6HNData->p11 * PObj6HNData->sgp + PObj6HNData->q11 * PObj6HNData->sgq;
	PObj6HNData->w2 = -0.5*PObj6HNData->p11*PObj6HNData->sgp + 0.5*PObj6HNData->q11*PObj6HNData->sgq + PObj6HNData->y11 * PObj6HNData->sgy;
	PObj6HNData->w3 = PObj6HNData->z11 * PObj6HNData->sgz / PObj6HNData->c;

	if (PObj6HNData->Fhn_activated) {
		Grad[0] = 2 * PObj6HNData->Fhn1*PObj6HNData->Fhn*PObj6HNData->w1;
		Grad[1] = 2 * PObj6HNData->Fhn1*PObj6HNData->Fhn*PObj6HNData->w2;
		Grad[2] = 2 * PObj6HNData->Fhn1*PObj6HNData->Fhn*PObj6HNData->w3;
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}
}

void DDFuncSqHess_Obj6_HN(double* Hess, void* PObjData) {
	Obj6HNData* PObj6HNData = (Obj6HNData*)PObjData;

	PObj6HNData->Fhn2 = PObj6HNData->Fhn1 * PObj6HNData->Fhn1;

	double s1 = -PObj6HNData->s + 1;
	double T0 = PObj6HNData->Fhn2 + s1 * PObj6HNData->Fhn3*PObj6HNData->Fhn;

	if (PObj6HNData->Fhn_activated) {
		Hess[0] = (PObj6HNData->w1*PObj6HNData->w1*T0 - s1 * PObj6HNData->Fhn1*PObj6HNData->Fhn*(PObj6HNData->p13 + PObj6HNData->q13)) * 2;
		Hess[1] = (PObj6HNData->w1*PObj6HNData->w2*T0 - s1 * PObj6HNData->Fhn1*PObj6HNData->Fhn*0.5*(-PObj6HNData->p13 + PObj6HNData->q13)) * 2;
		Hess[2] = PObj6HNData->w1 * PObj6HNData->w3*T0 * 2;
		Hess[3] = (PObj6HNData->w2*PObj6HNData->w2*T0 - s1 * PObj6HNData->Fhn1*PObj6HNData->Fhn*(0.25*(PObj6HNData->p13 + PObj6HNData->q13) + PObj6HNData->y13)) * 2;
		Hess[4] = PObj6HNData->w2 * PObj6HNData->w3*T0 * 2;
		Hess[5] = (PObj6HNData->w3*PObj6HNData->w3*T0 - s1 * PObj6HNData->z13*PObj6HNData->Fhn1*PObj6HNData->Fhn / (PObj6HNData->c*PObj6HNData->c)) * 2;
	}
	else {
		Hess[0] = 0;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 0;
		Hess[4] = 0;
		Hess[5] = 0;
	}
}

void DDFuncSqValu_Obj6_HN(double *val, void* PObjData) {
	Obj6HNData* PObj6HNData = (Obj6HNData*)PObjData;
	if (PObj6HNData->Fhn_activated) {
		*val = PObj6HNData->Fhn * PObj6HNData->Fhn;
	}
	else {
		*val = 0;
	}
}

void DDFuncSqValuIndep_Obj6_HN(const double* P, double *val, void* PObjData) {
	Obj6HNData* PObj6HNData = (Obj6HNData*)PObjData;

	double p1 = fabs(P[0] - 0.5*P[1]), q1 = fabs(P[0] + 0.5*P[1]), y1 = fabs(P[1]), z1 = fabs(P[2] / PObj6HNData->c);

	double Fhn0 = pow(p1, PObj6HNData->s) + pow(q1, PObj6HNData->s) + pow(y1, PObj6HNData->s) + pow(z1, PObj6HNData->s);
	double Fhn = pow(Fhn0, 1 / PObj6HNData->s) - 1;

	if (Fhn > 0) {
		*val = Fhn * Fhn;
	}
	else {
		*val = 0;
	}
}

/*-------------------------------   Object 6 Hexagonal nut  ------------------------------------------*/

/**************************************************************************************************/




/*-----------------------------------------------------------*/
inline double Sign(double y)
{
	double s;
	if (y > 0) {
		s = 1;
	}
	else if (y < 0) {
		s = -1;
	}
	else {
		s = 0;
	}
	return s;
}