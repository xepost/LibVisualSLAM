/*
 * SL_Distortion.cpp
 *
 *  Created on: 2010-11-6
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */

#include "SL_Distortion.h"
#include "math/SL_LinAlg.h"
#include <cmath>

double maxRangeNormPlane(int w , int h , const double* iK) {

	double r_maxw1 = fabs(iK[2]);
	double r_maxw2 = fabs(w * iK[0] + h * iK[1] + iK[2]);
	double r_maxw = r_maxw1 > r_maxw2 ? r_maxw1 : r_maxw2;

	double r_maxh1 = fabs(iK[5]);
	double r_maxh2 = fabs(w * iK[3] + h * iK[4] + iK[5]);
	double r_maxh = r_maxh1 > r_maxh2 ? r_maxh1 : r_maxh2;
	double r_max = r_maxw > r_maxh ? r_maxw : r_maxh;
	return r_max;
}
void invDistorParam(double r , const double *kc , double *k_ud) {

	const int NUM_SAMPLES = 20;

	int num_eqns = NUM_SAMPLES;
	int num_vars = 7;

	double *A = new double[num_eqns * num_vars];
	double *b = new double[num_eqns];

	for (int i = 0; i < NUM_SAMPLES; i++) {

		double t = i * r / (NUM_SAMPLES - 1);

		double tt = t * t;
		double ttt = tt * t;
		double ttttt = ttt * tt;
		double ttttttt = ttttt * tt;
		double t1 = t + kc[0] * ttt + kc[1] * ttttt + kc[4] * ttttttt;

		double tmp = t1;
		for (int j = 0; j < num_vars; j++) {
			A[i * num_vars + j] = tmp;
			tmp = tmp * t1;
		}
		b[i] = t;
	}
	dgelsyFor(num_eqns, num_vars, 1, A, b, k_ud);

	delete[] A;
	delete[] b;
}
void invDistorParam(int w , int h , const double* iK , const double* kc , double* k_ud) {
	double r_max = maxRangeNormPlane(w, h, iK);
	invDistorParam(r_max, kc, k_ud);
}
void undistorNormPoints(const double* iK , const double* k_ud , int n , const double* a , double* an) {

	int i;
	const double* pa = a;
	double* pan = an;
	double xn, yn;

	for (i = 0; i < n; ++i) {
		double x = *pa++;
		double y = *pa++;

		//get the normalized points
		double zn0 = iK[8];
		double xn0 = (x * iK[0] + y * iK[1] + iK[2]) / zn0;
		double yn0 = (y * iK[4] + iK[5]) / zn0;

		//remove the distortion
		double r = sqrt(xn0 * xn0 + yn0 * yn0);
		if (r == 0.0) {
			xn = xn0;
			yn = yn0;
		} else {
			double t = r;
			double a = 0;
			a += t * k_ud[0];
			t = t * r;
			a += t * k_ud[1];
			t = t * r;
			a += t * k_ud[2];
			t = t * r;
			a += t * k_ud[3];
			t = t * r;
			a += t * k_ud[4];
			t = t * r;
			a += t * k_ud[5];
			t = t * r;
			a += t * k_ud[6];
			double factor = a / r;
			xn = factor * xn0;
			yn = factor * yn0;
		}
		*pan = xn;
		pan++;
		*pan = yn;
		pan++;
	}
}

void undistorPoint(const double* K , const double* k_ud , const double* pt , double * undisPt) {
	double x = pt[0];
	double y = pt[1];
	double ptNorm[2];

	//get the normalized points
	double iK[9];
	getInvK(K, iK);
	double xn0, yn0;
	normPoint(iK, x, y, xn0, yn0);
	//remove the distortion
	double r = sqrt(xn0 * xn0 + yn0 * yn0);
	if (r == 0.0) {
		ptNorm[0] = xn0;
		ptNorm[1] = yn0;
	} else {
		double t = r;
		double a = 0;
		a += t * k_ud[0];
		t = t * r;
		a += t * k_ud[1];
		t = t * r;
		a += t * k_ud[2];
		t = t * r;
		a += t * k_ud[3];
		t = t * r;
		a += t * k_ud[4];
		t = t * r;
		a += t * k_ud[5];
		t = t * r;
		a += t * k_ud[6];
		double factor = a / r;
		ptNorm[0] = factor * xn0;
		ptNorm[1] = factor * yn0;
	}
	imagePoint(K, ptNorm, undisPt);
}

void undistorPointNew(double* K, double* kc, double* x, double* x_undist){
	double invK[9];
	double xn[2];
	getInvK(K, invK);
	normPoint(invK, x[0], x[1], xn[0], xn[1]);
    /*
	k1 = k(1);
    k2 = k(2);
    k3 = k(5);
    p1 = k(3);
    p2 = k(4);

    x = xd; 				% initial guess

    for kk=1:20,

        r_2 = sum(x.^2);
        k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
        delta_x = [2*p1*x(1,:).*x(2,:) + p2*(r_2 + 2*x(1,:).^2);
        p1 * (r_2 + 2*x(2,:).^2)+2*p2*x(1,:).*x(2,:)];
        x = (xd - delta_x)./(ones(2,1)*k_radial);

    end;
    */
	double k1 = kc[0];
	double k2 = kc[1];
	double k3 = kc[4];
	double p1 = kc[2];
	double p2 = kc[3];
	double x_temp[2];
	x_temp[0] = xn[0];
	x_temp[1] = xn[1];
	for (int i = 0; i < 20; i++){
		double r_2 = pow(x_temp[0], 2) + pow(x_temp[1], 2);
		double k_radial =  1 + k1 * r_2 + k2 * r_2 * r_2 + k3 * pow(r_2, 3);
		double delta_x[2];
		delta_x[0] = 2*p1*x_temp[0]*x_temp[1] + p2*(r_2 + 2*x_temp[0] * x_temp[0]);
		delta_x[1] = p1* (r_2 + 2*x_temp[1] * x_temp[1])+2*p2*x_temp[0]*x_temp[1];
		x_temp[0] = (xn[0] - delta_x[0]) / k_radial;
		x_temp[1] = (xn[1] - delta_x[1]) / k_radial;
	}
	imagePoint(K, x_temp, x_undist);
}

void normPoints(const double* iK , int n , const double* a , double *an) {
	int i;
	const double* pa = a;
	double* pan = an;

	for (i = 0; i < n; ++i) {
		double x = *pa++;
		double y = *pa++;

		double zn = iK[8];
		double xn = (x * iK[0] + y * iK[1] + iK[2]) / zn;
		double yn = (y * iK[4] + iK[5]) / zn;

		*pan = xn;
		pan++;
		*pan = yn;
		pan++;
	}
}
void normPoint(const double* iK , const double* a , double *an) {
	double zn = iK[8];
	an[0] = (a[0] * iK[0] + a[1] * iK[1] + iK[2]) / zn;
	an[1] = (a[1] * iK[4] + iK[5]) / zn;
}

void normPoint(const double* iK , const double x , const double y , double& xNorm , double& yNorm) {
	double zn = iK[8];
	xNorm = (x * iK[0] + y * iK[1] + iK[2]) / zn;
	yNorm = (y * iK[4] + iK[5]) / zn;
}

void getInvK(const double* K , double* invK) {
	invK[0] = 1.0 / K[0];
	invK[1] = -K[1] / (K[0] * K[4]);
	invK[2] = -K[2] / K[0] + K[1] * K[5] / (K[0] * K[4]);
	invK[3] = 0;
	invK[4] = 1.0 / K[4];
	invK[5] = -K[5] / K[4];
	invK[6] = 0;
	invK[7] = 0;
	invK[8] = 1;
}
void imagePoints(const double* K , int n , const double* an , double* a) {
	int i;
	const double* pan = an;
	double* pa = a;

	for (i = 0; i < n; ++i) {
		double xn = *pan++;
		double yn = *pan++;

		double z = K[8];
		double x = (xn * K[0] + yn * K[1] + K[2]) / z;
		double y = (yn * K[4] + K[5]) / z;

		*pa = x;
		pa++;
		*pa = y;
		pa++;
	}
}

void imagePoint(const double* K , const double* an , double* a) {
	double xn = an[0];
	double yn = an[1];

	double z = K[8];
	double x = (xn * K[0] + yn * K[1] + K[2]) / z;
	double y = (yn * K[4] + K[5]) / z;

	a[0] = x;
	a[1] = y;
}


void imagePoint(const double* K , const double* kc , const double* an , double* a) {
	double xn = an[0];
	double yn = an[1];

	double r = sqrt(xn * xn + yn * yn);
	double rr = r * r;
	double rrrr = rr * rr;
	double rrrrrr = rrrr * rr;

	double factor_a = (1 + kc[0] * rr + kc[1] * rrrr + kc[4] * rrrrrr);
	double factor_bx = 2 * kc[2] * xn * yn + kc[3] * (rr + 2 * xn * xn);
	double factor_by = kc[2] * (rr + 2 * yn * yn) + 2 * kc[3] * xn * yn;

	double xnd = xn * factor_a + factor_bx;
	double ynd = yn * factor_a + factor_by;

	double z = K[8];
	double x = (K[0] * xnd + K[1] * ynd + K[2]) / z;
	double y = (K[4] * ynd + K[5]) / z;

	a[0] = x;
	a[1] = y;
}

void imagePoints(const double* K , const double* kc , int n , const double* an , double* a) {
	int i;
	const double* pan = an;
	double* pa = a;

	for (i = 0; i < n; ++i) {
		double xn = *pan++;
		double yn = *pan++;

		double r = sqrt(xn * xn + yn * yn);
		double rr = r * r;
		double rrrr = rr * rr;
		double rrrrrr = rrrr * rr;

		double factor_a = (1 + kc[0] * rr + kc[1] * rrrr + kc[4] * rrrrrr);
		double factor_bx = 2 * kc[2] * xn * yn + kc[3] * (rr + 2 * xn * xn);
		double factor_by = kc[2] * (rr + 2 * yn * yn) + 2 * kc[3] * xn * yn;

		double xnd = xn * factor_a + factor_bx;
		double ynd = yn * factor_a + factor_by;

		double z = K[8];
		double x = (K[0] * xnd + K[1] * ynd + K[2]) / z;
		double y = (K[4] * ynd + K[5]) / z;

		*pa = x;
		pa++;
		*pa = y;
		pa++;
	}
}
