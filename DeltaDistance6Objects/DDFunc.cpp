#include "pch.h"
#include "config.h"
#include "DDFunc.h"


/**************************************************************************************************/

/*-------------------------------   Object 1 Mouse  ------------------------------------------*/
void DDFuncSqValu_Obj1_MS(const double* P, double *val) {
	int k = 4;
	double F[4];

	double x1, y1, z1;

	double a0 = 0.7, b0 = 1.5, c0 = 1.2;
	double a1 = 1, b1 = 1.1, c1 = 1.2;
	double a2 = 1.5, b2 = 1.2, c2 = 0.6;

	x1 = P[0] / a0, y1 = P[1] / b0, z1 = P[2] / c0;
	F[0] = sqrt(x1 * x1 + y1 * y1 + z1 * z1) - 1;
	
	x1 = P[0] / a1, y1 = P[1] / b1, z1 = P[2] / c1;
	F[1] = sqrt(x1 * x1 + y1 * y1 + z1 * z1) - 1;

	x1 = P[0] / a2, y1 = P[1] / b2, z1 = P[2] / c2;
	F[2] = sqrt(x1 * x1 + y1 * y1 + z1 * z1) - 1;

	F[3] = -P[2];

	*val = 0;
	for (int i = 0;i < k;i++) {
		if (F[i] > 0) {
			*val += F[i] * F[i];
		}
	}
}

void DDFuncSqGrad_Obj1_MS(const double* P, double* Grad) {
	int k = 4;
	double F[4];
	for (int i = 0;i < 3;i++) {
		Grad[i] = 0;
	}
	
	double x1, y1, z1;

	double a0 = 0.7, b0 = 1.5, c0 = 1.2;
	double a1 = 1, b1 = 1.1, c1 = 1.2;
	double a2 = 1.5, b2 = 1.2, c2 = 0.6;

	x1 = P[0] / a0, y1 = P[1] / b0, z1 = P[2] / c0;
	double F01 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
	F[0] = F01 - 1;
	if (F[0] > 0) {
		Grad[0] += x1 * F[0] * 2 / (a0*F01);
		Grad[1] += y1 * F[0] * 2 / (b0*F01);
		Grad[2] += z1 * F[0] * 2 / (c0*F01);
	}

	x1 = P[0] / a1, y1 = P[1] / b1, z1 = P[2] / c1;
	double F11 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
	F[1] = F11 - 1;
	if (F[1] > 0) {
		Grad[0] += x1 * F[1] * 2 / (a1*F11);
		Grad[1] += y1 * F[1] * 2 / (b1*F11);
		Grad[2] += z1 * F[1] * 2 / (c1*F11);
	}

	x1 = P[0] / a2, y1 = P[1] / b2, z1 = P[2] / c2;
	double F21 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
	F[2] = F21 - 1;
	if (F[2] > 0) {
		Grad[0] += x1 * F[2] * 2 / (a2*F21);
		Grad[1] += y1 * F[2] * 2 / (b2*F21);
		Grad[2] += z1 * F[2] * 2 / (c2*F21);
	}

	F[3] = -P[2];
	if (F[3] > 0) {
		Grad[2] -= 2 * F[3];
	}
}

void DDFuncSqHess_Obj1_MS(const double* P, double* Hess) {
	int k = 4;
	double F[4];
	for (int i = 0;i < 6;i++) {
		Hess[i] = 0;
	}
	
	double x1, y1, z1;

	double a0 = 0.7, b0 = 1.5, c0 = 1.2;
	double a1 = 1, b1 = 1.1, c1 = 1.2;
	double a2 = 1.5, b2 = 1.2, c2 = 0.6;

	x1 = P[0] / a0, y1 = P[1] / b0, z1 = P[2] / c0;
	double F00 = x1 * x1 + y1 * y1 + z1 * z1;
	double F01 = sqrt(F00);
	double F02 = F00 * F01;
	F[0] = F01 - 1;
	double F03 = 1 / F00 - F[0] / F02;
	if (F[0] > 0) {
		Hess[0] += (x1*x1*F03 + F[0] / F01) * 2 / (a0*a0);
		Hess[1] += 2 * x1*y1*F03 / (a0*b0);
		Hess[2] += 2 * x1*z1*F03 / (a0*c0);
		Hess[3] += (y1*y1*F03 + F[0] / F01) * 2 / (b0*b0);
		Hess[4] += 2 * y1*z1*F03 / (b0*c0);
		Hess[5] += (z1*z1*F03 + F[0] / F01) * 2 / (c0*c0);
	}

	x1 = P[0] / a1, y1 = P[1] / b1, z1 = P[2] / c1;
	double F10 = x1 * x1 + y1 * y1 + z1 * z1;
	double F11 = sqrt(F10);
	double F12 = F10 * F11;
	F[1] = F11 - 1;
	double F13 = 1 / F10 - F[1] / F12;
	if (F[1] > 0) {
		Hess[0] += (x1*x1*F13 + F[1] / F11) * 2 / (a1*a1);
		Hess[1] += 2 * x1*y1*F13 / (a1*b1);
		Hess[2] += 2 * x1*z1*F13 / (a1*c1);
		Hess[3] += (y1*y1*F13 + F[1] / F11) * 2 / (b1*b1);
		Hess[4] += 2 * y1*z1*F13 / (b1*c1);
		Hess[5] += (z1*z1*F13 + F[1] / F11) * 2 / (c1*c1);
	}

	x1 = P[0] / a2, y1 = P[1] / b2, z1 = P[2] / c2;
	double F20 = x1 * x1 + y1 * y1 + z1 * z1;
	double F21 = sqrt(F20);
	double F22 = F20 * F21;
	F[2] = F21 - 1;
	double F23 = 1 / F20 - F[2] / F22;
	if (F[2] > 0) {
		Hess[0] += (x1*x1*F23 + F[2] / F21) * 2 / (a2*a2);
		Hess[1] += 2 * x1*y1*F23 / (a2*b2);
		Hess[2] += 2 * x1*z1*F23 / (a2*c2);
		Hess[3] += (y1*y1*F23 + F[2] / F21) * 2 / (b2*b2);
		Hess[4] += 2 * y1*z1*F23 / (b2*c2);
		Hess[5] += (z1*z1*F23 + F[2] / F21) * 2 / (c2*c2);
	}

	F[3] = -P[2];
	if (F[3] > 0) {
		Hess[5] += 2;
	}

}
/*-------------------------------   Object 1 Mouse  --------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 2 Pencil  ------------------------------------------*/
void DDFuncSqValu_Obj2_PC(const double* P, double *val) {
	int k = 5;
	double F[5];

	double h = 2, r0 = 0.6;
	double l = sqrt(h * h + r0 * r0);
	double t0 = r0 / h;
	double c0 = h / l;
	double s0 = r0 / l;
	double r_sq = P[0] * P[0] + P[1] * P[1];
	double r = sqrt(r_sq);
	double sqrt3 = sqrt(3);
	double a1 = 0.2, b1 = 0.4;
	double h2 = 1;

	if (P[2] - t0 * r > h) {
		*val = r_sq + (P[2] - h)*(P[2] - h);
	}
	else if ((P[2] - t0 * r <= h) && (P[2] + r / t0 > h)) {
		double b_c2 = c0 * r + s0 * P[2] - c0 * r0;
		*val = b_c2 * b_c2;
	}
	else {
		*val = 0;
	}

	double x1, y1;
	double xp, yp;
	
	x1 = P[0] / a1, y1 = P[1] / b1;
	F[1] = sqrt(x1 * x1 + y1 * y1) - 1;
	
	xp = P[0] / 2 + sqrt3 * P[1] / 2, yp = -sqrt3 * P[0] / 2 + P[1] / 2;
	x1 = xp / a1, y1 = yp / b1;
	F[2] = sqrt(x1 * x1 + y1 * y1) - 1;

	xp = P[0] / 2 - sqrt3 * P[1] / 2, yp = sqrt3 * P[0] / 2 + P[1] / 2;
	x1 = xp / a1, y1 = yp / b1;
	F[3] = sqrt(x1 * x1 + y1 * y1) - 1;
	

	F[4] = -P[2] - h2;

	for (int i = 1;i < k;i++) {
		if (F[i] > 0) {
			*val += F[i] * F[i];
		}
	}
}

void DDFuncSqGrad_Obj2_PC(const double* P, double* Grad) {
	int k = 5;
	double F[5];

	double h = 2, r0 = 0.6;
	double l = sqrt(h * h + r0 * r0);
	double t0 = r0 / h;
	double c0 = h / l;
	double s0 = r0 / l;
	double r_sq = P[0] * P[0] + P[1] * P[1];
	double r = sqrt(r_sq);
	double sqrt3 = sqrt(3);
	double a1 = 0.2, b1 = 0.4;
	double h2 = 1;

	if (P[2] - t0 * r > h) {
		Grad[0] = 2 * P[0];
		Grad[1] = 2 * P[1];
		Grad[2] = 2 * (P[2] - h);
	}
	else if ((P[2] - t0 * r <= h) && (P[2] + r / t0 > h)) {
		double b_c2 = c0 * r + s0 * P[2] - c0 * r0;
		Grad[0] = 2 * b_c2*c0*P[0] / r;
		Grad[1] = 2 * b_c2*c0*P[1] / r;
		Grad[2] = 2 * b_c2*s0;
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}

	double x1, y1;
	double xp, yp;
	
	x1 = P[0] / a1, y1 = P[1] / b1;
	double F11 = sqrt(x1 * x1 + y1 * y1);
	F[1] = F11 - 1;
	if (F[1] > 0) {
		Grad[0] += x1 * F[1] * 2 / (a1*F11);
		Grad[1] += y1 * F[1] * 2 / (b1*F11);
	}
	
	xp = P[0] / 2 + sqrt3 * P[1] / 2, yp = -sqrt3 * P[0] / 2 + P[1] / 2;
	x1 = xp / a1, y1 = yp / b1;
	double F21 = sqrt(x1 * x1 + y1 * y1);
	F[2] = F21 - 1;
	if (F[2] > 0) {
		double Gradp[2];
		Gradp[0] = x1 * F[2] * 2 / (a1*F21);
		Gradp[1] = y1 * F[2] * 2 / (b1*F21);
		Grad[0] += Gradp[0] / 2 - sqrt3 * Gradp[1] / 2;
		Grad[1] += sqrt3 * Gradp[0] / 2 + Gradp[1] / 2;
	}

	xp = P[0] / 2 - sqrt3 * P[1] / 2, yp = sqrt3 * P[0] / 2 + P[1] / 2;
	x1 = xp / a1, y1 = yp / b1;
	double F31 = sqrt(x1 * x1 + y1 * y1);
	F[3] = F31 - 1;
	if (F[3] > 0) {
		double Gradp[2];
		Gradp[0] = x1 * F[3] * 2 / (a1*F31);
		Gradp[1] = y1 * F[3] * 2 / (b1*F31);
		Grad[0] += Gradp[0] / 2 + sqrt3 * Gradp[1] / 2;
		Grad[1] += -sqrt3 * Gradp[0] / 2 + Gradp[1] / 2;
	}
	
	F[4] = -P[2] - h2;
	if (F[4] > 0) {
		Grad[2] -= 2 * F[4];
	}

}

void DDFuncSqHess_Obj2_PC(const double* P, double* Hess) {
	int k = 5;
	double F[5];

	double h = 2, r0 = 0.6;
	double l = sqrt(h * h + r0 * r0);
	double t0 = r0 / h;
	double c0 = h / l;
	double s0 = r0 / l;
	double r_sq = P[0] * P[0] + P[1] * P[1];
	double r = sqrt(r_sq);
	double r_cu = r * r_sq;
	double sqrt3 = sqrt(3);
	double a1 = 0.2, b1 = 0.4;
	double h2 = 1;

	if (P[2] - t0 * r > h) {
		Hess[0] = 2;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 2;
		Hess[4] = 0;
		Hess[5] = 2;
	}
	else if ((P[2] - t0 * r <= h) && (P[2] + r / t0 > h)) {
		double b_c2 = c0 * r + s0 * P[2] - c0 * r0;
		Hess[0] = 2 * c0*c0*P[0] * P[0] / r_sq - 2 * c0*P[0] * P[0] * b_c2 / r_cu + 2 * c0*b_c2 / r;
		Hess[1] = 2 * c0*c0*P[0] * P[1] / r_sq - 2 * c0*P[0] * P[1] * b_c2 / r_cu;
		Hess[2] = 2 * c0*s0*P[0] / r;
		Hess[3] = 2 * c0*c0*P[1] * P[1] / r_sq - 2 * c0*P[1] * P[1] * b_c2 / r_cu + 2 * c0*b_c2 / r;
		Hess[4] = 2 * c0*s0*P[1] / r;
		Hess[5] = 2 * s0*s0;
	}
	else {
		Hess[0] = 0;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 0;
		Hess[4] = 0;
		Hess[5] = 0;
	}

	double x1, y1;
	double xp, yp;
	
	x1 = P[0] / a1, y1 = P[1] / b1;
	double F10 = x1 * x1 + y1 * y1;
	double F11 = sqrt(F10);
	double F12 = F10 * F11;
	F[1] = F11 - 1;
	double F13 = 1 / F10 - F[1] / F12;
	if (F[1] > 0) {
		Hess[0] += (x1*x1*F13 + F[1] / F11) * 2 / (a1*a1);
		Hess[1] += 2 * x1*y1*F13 / (a1*b1);
		Hess[3] += (y1*y1*F13 + F[1] / F11) * 2 / (b1*b1);
	}
	
	xp = P[0] / 2 + sqrt3 * P[1] / 2, yp = -sqrt3 * P[0] / 2 + P[1] / 2;
	x1 = xp / a1, y1 = yp / b1;
	double F20 = x1 * x1 + y1 * y1;
	double F21 = sqrt(F20);
	double F22 = F20 * F21;
	F[2] = F21 - 1;
	double F23 = 1 / F20 - F[2] / F22;
	if (F[2] > 0) {
		double Hessp[3];
		Hessp[0] = (x1*x1*F23 + F[2] / F21) * 2 / (a1*a1);
		Hessp[1] = 2 * x1*y1*F23 / (a1*b1);
		Hessp[2] = (y1*y1*F23 + F[2] / F21) * 2 / (b1*b1);

		Hess[0] += Hessp[0] / 4 - sqrt3 * Hessp[1] / 2 + 3 * Hessp[2] / 4;
		Hess[1] += sqrt3 * Hessp[0] / 4 - Hessp[1] / 2 - sqrt3 * Hessp[2] / 4;
		Hess[3] += 3 * Hessp[0] / 4 + sqrt3 * Hessp[1] / 2 + Hessp[2] / 4;
	}

	xp = P[0] / 2 - sqrt3 * P[1] / 2, yp = sqrt3 * P[0] / 2 + P[1] / 2;
	x1 = xp / a1, y1 = yp / b1;
	double F30 = x1 * x1 + y1 * y1;
	double F31 = sqrt(F30);
	double F32 = F30 * F31;
	F[3] = F31 - 1;
	double F33 = 1 / F30 - F[3] / F32;
	if (F[3] > 0) {
		double Hessp[3];
		Hessp[0] = (x1*x1*F33 + F[3] / F31) * 2 / (a1*a1);
		Hessp[1] = 2 * x1*y1*F33 / (a1*b1);
		Hessp[2] = (y1*y1*F33 + F[3] / F31) * 2 / (b1*b1);

		Hess[0] += Hessp[0] / 4 + sqrt3 * Hessp[1] / 2 + 3 * Hessp[2] / 4;
		Hess[1] += -sqrt3 * Hessp[0] / 4 - Hessp[1] / 2 + sqrt3 * Hessp[2] / 4;
		Hess[3] += 3 * Hessp[0] / 4 - sqrt3 * Hessp[1] / 2 + Hessp[2] / 4;
	}
	
	F[4] = -P[2] - h2;
	if (F[4] > 0) {
		Hess[5] += 2;
	}
}
/*-------------------------------  Object 2 Pencil --------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 3 Diamond  ------------------------------------------*/
void DDFuncSqValu_Obj3_DM(const double* P, double *val) {
	double a = 0.6, b = 1.3, c = 1.8;
	double n1 = 1.4, n2 = 1.6, n3 = 1.5;
	
	double x1 = fabs(P[0] / a), y1 = fabs(P[1] / b), z1 = fabs(P[2] / c);
	double F = pow(x1, n1) + pow(y1, n2) + pow(z1, n3) - 1;
	if (F > 0) {
		*val = F * F;
	}
	else {
		*val = 0;
	}
}

void DDFuncSqGrad_Obj3_DM(const double* P, double* Grad) {
	double a = 0.6, b = 1.3, c = 1.8;
	double n1 = 1.4, n2 = 1.6, n3 = 1.5;
	double x1 = fabs(P[0] / a), y1 = fabs(P[1] / b), z1 = fabs(P[2] / c);
	double sx = Sign(P[0]), sy = Sign(P[1]), sz = Sign(P[2]);
	double x11 = pow(x1, n1 - 1), y11 = pow(y1, n2 - 1), z11 = pow(z1, n3 - 1);
	double F = x11 * x1 + y11 * y1 + z11 * z1 - 1;
	if (F > 0) {
		Grad[0] = 2 * sx * n1*F*x11 / a;
		Grad[1] = 2 * sy * n2*F*y11 / b;
		Grad[2] = 2 * sz * n3*F*z11 / c;
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}
}

void DDFuncSqHess_Obj3_DM(const double* P, double* Hess) {
	double a = 0.6, b = 1.3, c = 1.8;
	double n1 = 1.4, n2 = 1.6, n3 = 1.5;
	double x1 = fabs(P[0] / a), y1 = fabs(P[1] / b), z1 = fabs(P[2] / c);
	double sx = Sign(P[0]), sy = Sign(P[1]), sz = Sign(P[2]);
	double x13 = pow(x1, n1 - 2), y13 = pow(y1, n2 - 2), z13 = pow(z1, n3 - 2);
	double x11 = x13 * x1, y11 = y13 * y1, z11 = z13 * z1;
	double F = x11 * x1 + y11 * y1 + z11 * z1 - 1;
	if (F > 0) {
		Hess[0] = 2 * n1*(n1*x11*x11 + (n1 - 1)*F*x13) / (a*a);
		Hess[1] = 2 * sx*sy * n1*n2*x11*y11 / (a*b);
		Hess[2] = 2 * sx*sz * n1*n3*x11*z11 / (a*c);
		Hess[3] = 2 * n2*(n2*y11*y11 + (n2 - 1)*F*y13) / (b*b);
		Hess[4] = 2 * sy*sz * n2*n3*y11*z11 / (b*c);
		Hess[5] = 2 * n3*(n3*z11*z11 + (n3 - 1)*F*z13) / (c*c);
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
/*-------------------------------   Object 3 Diamond  ------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 4 Piece of cake  ------------------------------------------*/
void DDFuncSqValu_Obj4_CK(const double* P, double *val) {
	double a = 2.1, b = 2.1, c = 0.5;
	double s = 4;
	double theta = PI / 3;
	double s0 = sin(theta), c0 = cos(theta), t0 = tan(theta);

	double x1 = P[0] / a, y1 = P[1] / b, z1 = P[2] / c;
	double Fxy0 = x1 * x1 + y1 * y1;

	double Fse = pow(pow(Fxy0, s) + pow(z1*z1, s), 0.5 / s) - 1;
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
	F[2] = c0 * P[0] - s0 * P[1];

	for (int i = 0;i < k;i++) {
		if (F[i] > 0) {
			*val += F[i] * F[i];
		}
	}

}

void DDFuncSqGrad_Obj4_CK(const double* P, double* Grad) {
	double a = 2.1, b = 2.1, c = 0.5;
	double s = 4;
	double theta = PI / 3;
	double s0 = sin(theta), c0 = cos(theta), t0 = tan(theta);

	double x1 = P[0] / a, y1 = P[1] / b, z1 = P[2] / c;
	double Fxy0 = x1 * x1 + y1 * y1;
	double Fxy1 = pow(Fxy0, s - 1);
	double z13 = pow(z1*z1, s - 1);
	double z11 = z13 * z1;
	double Ftr0 = Fxy1 * Fxy0 + z11 * z1;
	double Ftr1 = pow(Ftr0, 0.5 / s - 1);
	double Ftr = Ftr0 * Ftr1 - 1;

	if (Ftr > 0) {
		Grad[0] = 2 * x1 * Fxy1*Ftr1*Ftr / a;
		Grad[1] = 2 * y1 * Fxy1*Ftr1*Ftr / b;
		Grad[2] = 2 * z11 * Ftr1*Ftr / c;
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}

	double F[3];
	int k = 3;
	F[0] = -P[2];
	F[1] = -P[0];
	F[2] = c0 * P[0] - s0 * P[1];
	if (F[0] > 0) {
		Grad[2] += 2 * P[2];
	}
	if (F[1] > 0) {
		Grad[0] += 2 * P[0];
	}
	if (F[2] > 0) {
		Grad[0] += 2 * F[2] * c0;
		Grad[1] -= 2 * F[2] * s0;
	}

}

void DDFuncSqHess_Obj4_CK(const double* P, double* Hess) {
	double a = 2.1, b = 2.1, c = 0.5;
	double s = 4;
	double theta = PI / 3;
	double s0 = sin(theta), c0 = cos(theta), t0 = tan(theta);

	double x1 = P[0] / a, y1 = P[1] / b, z1 = P[2] / c;
	double Fxy0 = x1 * x1 + y1 * y1;
	double Fxy3 = pow(Fxy0, s - 2);
	double Fxy1 = Fxy3 * Fxy0;
	double Fxy2 = Fxy1 * Fxy1;
	double z13 = pow(z1*z1, s - 1);
	double z11 = z13 * z1;
	double z12 = z11 * z11;
	double Ftr0 = Fxy1 * Fxy0 + z11 * z1;
	double Ftr3 = pow(Ftr0, 0.5 / s - 2);
	double Ftr1 = Ftr3 * Ftr0;
	double Ftr2 = Ftr1 * Ftr1;
	double Ftr = Ftr0 * Ftr1 - 1;

	double s1 = -2 * s + 1;
	double s2 = 2 * s - 2;

	double T1 = Ftr2 + s1 * Ftr3*Ftr;
	double T2 = s2 * Fxy3*Ftr1*Ftr;

	if (Ftr > 0) {
		Hess[0] = (x1*x1*(Fxy2*T1 + T2) + Fxy1 * Ftr1*Ftr) * 2 / (a*a);
		Hess[1] = x1 * y1*(Fxy2*T1 + T2) * 2 / (a*b);
		Hess[2] = x1 * Fxy1*z11*T1 * 2 / (a*c);
		Hess[3] = (y1*y1*(Fxy2*T1 + T2) + Fxy1 * Ftr1*Ftr) * 2 / (b*b);
		Hess[4] = y1 * Fxy1*z11*T1 * 2 / (b*c);
		Hess[5] = (z12*T1 - s1 * z13*Ftr1*Ftr) * 2 / (c*c);
	}
	else {
		Hess[0] = 0;
		Hess[1] = 0;
		Hess[2] = 0;
		Hess[3] = 0;
		Hess[4] = 0;
		Hess[5] = 0;
	}

	double F[3];
	int k = 3;
	F[0] = -P[2];
	F[1] = -P[0];
	F[2] = c0 * P[0] - s0 * P[1];

	if (F[0] > 0) {
		Hess[5] += 2;
	}
	if (F[1] > 0) {
		Hess[0] += 2;
	}
	if (F[2] > 0) {
		Hess[0] += 2 * c0 * c0;
		Hess[1] -= 2 * c0 * s0;
		Hess[3] += 2 * s0 * s0;
	}
}
/*-------------------------------   Object 4 Piece of cake  ------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 5 Pillow  ------------------------------------------*/
void DDFuncSqValu_Obj5_PL(const double* P, double *val) {
	double a = 0.7, b = 1.2, c = 0.4;
	double s = 6;

	double x1 = P[0] / a, y1 = P[1] / b, z1 = P[2] / c;
	double Fxy0 = pow(x1*x1, s) + pow(y1*y1, s);
	double Fpl = sqrt(pow(Fxy0, 1 / s) + z1 * z1) - 1;

	if (Fpl > 0) {
		*val = Fpl * Fpl;
	}
	else {
		*val = 0;
	}
}

void DDFuncSqGrad_Obj5_PL(const double* P, double* Grad) {
	double a = 0.7, b = 1.2, c = 0.4;
	double s = 6;

	double x1 = P[0] / a, y1 = P[1] / b, z1 = P[2] / c;
	double x11 = pow(x1*x1, s - 1), y11 = pow(y1*y1, s - 1);
	double Fxy0 = pow(x1*x1, s) + pow(y1*y1, s);
	double Fxy1 = pow(Fxy0, 1/s - 1);

	double Fpl0 = Fxy0 * Fxy1 + z1 * z1;
	double Fpl1 = sqrt(Fpl0);
	double Fpl = Fpl1 - 1;

	if (Fpl > 0) {
		Grad[0] = 2 * x1*x11*Fxy1*Fpl / (a*Fpl1);
		Grad[1] = 2 * y1*y11*Fxy1*Fpl / (b*Fpl1);
		Grad[2] = 2 * z1*Fpl / (c*Fpl1);
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}
}

void DDFuncSqHess_Obj5_PL(const double* P, double* Hess) {
	double a = 0.7, b = 1.2, c = 0.4;
	double s = 6;

	double x1 = P[0] / a, y1 = P[1] / b, z1 = P[2] / c;
	double x13 = pow(x1*x1, s - 2), y13 = pow(y1*y1, s - 2);
	double x11 = x13 * x1*x1, y11 = y13 * y1*y1;
	double x12 = x11 * x11, y12 = y11 * y11;
	double Fxy0 = pow(x1*x1, s) + pow(y1*y1, s);
	double Fxy3 = pow(Fxy0, 1 / s - 2);
	double Fxy1 = Fxy3 * Fxy0;
	double Fxy2 = Fxy1 * Fxy1;

	double Fpl0 = Fxy0 * Fxy1 + z1 * z1;
	double Fpl1 = sqrt(Fpl0);
	double Fpl2 = Fpl0 * Fpl1;
	double Fpl = Fpl1 - 1;

	double s1 = -2 * s + 2;

	double T0 = 1 / Fpl0 - Fpl / Fpl2;
	double T1 = Fxy2 * T0 + s1 * Fxy3 * Fpl / Fpl1;

	if (Fpl > 0) {
		Hess[0] = (x1*x1*(x12*T1 - s1*x13 * Fxy1*Fpl/Fpl1) + x11 * Fxy1*Fpl / Fpl1) * 2 / (a*a);
		Hess[1] = x1 * y1*x11*y11*T1 * 2 / (a*b);
		Hess[2] = x1 * z1*x11*Fxy1*T0 * 2 / (a*c);
		Hess[3] = (y1*y1*(y12*T1 - s1*y13 * Fxy1*Fpl/Fpl1) + y11 * Fxy1*Fpl / Fpl1) * 2 / (b*b);
		Hess[4] = y1 * z1*y11*Fxy1*T0 * 2 / (b*c);
		Hess[5] = (z1*z1*T0 + Fpl / Fpl1) * 2 / (c*c);
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
/*-------------------------------   Object 5 Pillow  ------------------------------------------*/

/**************************************************************************************************/

/*-------------------------------   Object 6 Hexagonal nut  ------------------------------------------*/
void DDFuncSqValu_Obj6_HN(const double* P, double *val) {
	double c = 0.4;
	double s = 16;

	double p1 = fabs(P[0] - 0.5*P[1]), q1 = fabs(P[0] + 0.5*P[1]), y1 = fabs(P[1]), z1 = fabs(P[2] / c);

	double Fhn0 = pow(p1, s) + pow(q1, s) + pow(y1, s) + pow(z1, s);
	double Fhn = pow(Fhn0, 1 / s) - 1;

	if (Fhn > 0) {
		*val = Fhn * Fhn;
	}
	else {
		*val = 0;
	}
}

void DDFuncSqGrad_Obj6_HN(const double* P, double* Grad) {
	double c = 0.4;
	double s = 16;

	double p1 = fabs(P[0] - 0.5*P[1]), q1 = fabs(P[0] + 0.5*P[1]), y1 = fabs(P[1]), z1 = fabs(P[2] / c);
	double p11 = pow(p1, s - 1), q11 = pow(q1, s - 1), y11 = pow(y1, s - 1), z11 = pow(z1, s - 1);
	double sgp = Sign(P[0] - 0.5*P[1]), sgq = Sign(P[0] + 0.5*P[1]), sgy = Sign(P[1]), sgz = Sign(P[2]);
	
	double Fhn0 = p11 * p1 + q11 * q1 + y11 * y1 + z11 * z1;
	double Fhn1 = pow(Fhn0, -1 + 1 / s);
	double Fhn = Fhn0 * Fhn1 - 1;

	double w1 = p11 * sgp + q11 * sgq;
	double w2 = -0.5*p11*sgp + 0.5*q11*sgq + y11 * sgy;
	double w3 = z11 * sgz / c;

	if (Fhn > 0) {
		Grad[0] = 2 * Fhn1*Fhn*w1;
		Grad[1] = 2 * Fhn1*Fhn*w2;
		Grad[2] = 2 * Fhn1*Fhn*w3;
	}
	else {
		Grad[0] = 0;
		Grad[1] = 0;
		Grad[2] = 0;
	}
}

void DDFuncSqHess_Obj6_HN(const double* P, double* Hess) {
	double c = 0.4;
	double s = 16;

	double p1 = fabs(P[0] - 0.5*P[1]), q1 = fabs(P[0] + 0.5*P[1]), y1 = fabs(P[1]), z1 = fabs(P[2] / c);
	double p13 = pow(p1, s - 2), q13 = pow(q1, s - 2), y13 = pow(y1, s - 2), z13 = pow(z1, s - 2);
	double p11 = p13 * p1, q11 = q13 * q1, y11 = y13 * y1, z11 = z13 * z1;
	double sgp = Sign(P[0] - 0.5*P[1]), sgq = Sign(P[0] + 0.5*P[1]), sgy = Sign(P[1]), sgz = Sign(P[2]);

	double Fhn0 = p11 * p1 + q11 * q1 + y11 * y1 + z11 * z1;
	double Fhn3 = pow(Fhn0, -2 + 1 / s);
	double Fhn1 = Fhn3 * Fhn0;
	double Fhn2 = Fhn1 * Fhn1;
	double Fhn = Fhn0 * Fhn1 - 1;

	double w1 = p11 * sgp + q11 * sgq;
	double w2 = -0.5*p11*sgp + 0.5*q11*sgq + y11 * sgy;
	double w3 = z11 * sgz / c;

	double s1 = -s + 1;

	double T0 = Fhn2 + s1 * Fhn3*Fhn;

	if (Fhn > 0) {
		Hess[0] = (w1*w1*T0 - s1 * Fhn1*Fhn*(p13 + q13)) * 2;
		Hess[1] = (w1*w2*T0 - s1 * Fhn1*Fhn*0.5*(-p13 + q13)) * 2;
		Hess[2] = w1 * w3*T0 * 2;
		Hess[3] = (w2*w2*T0 - s1 * Fhn1*Fhn*(0.25*(p13 + q13) + y13)) * 2;
		Hess[4] = w2 * w3*T0 * 2;;
		Hess[5] = (w3*w3*T0 - s1 * z13*Fhn1*Fhn / (c*c)) * 2;
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