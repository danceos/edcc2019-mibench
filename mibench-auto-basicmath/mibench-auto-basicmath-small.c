#include "snipmath.h"
#include <math.h>
#include "result.h"
#include "config.h"
//#include <stdio.h>


int main(void) {

	double  a1 = 1.0, b1 = -10.5, c1 = 32.0, d1 = -30.0;
	double  a2 = 1.0, b2 = -4.5, c2 = 17.0, d2 = -30.0;
	double  a3 = 1.0, b3 = -3.5, c3 = 22.0, d3 = -31.0;
	double  a4 = 1.0, b4 = -13.7, c4 = 1.0, d4 = -35.0;

	double  x[3];
	double X;
	int solutions;
	int i;
	unsigned long l = 0x3fed0169L;
	struct int_sqrt q;
	double y;
//	long n = 0;

	*FI_START = 1;
/*
	int f1 = 0xdead;
	int f2 = 0xf00;
	hashIncrResultData(&f1, 4);
	hashIncrResultData(&f2, 4);
	hashIncrResultData(0, 1016);
	f1 = 0xb00b0000;
	f2 = 0x311;
	hashIncrResultData(&f1, 4);
	hashIncrResultData(&f2, 4);


	return 0;
*/
	/* solve soem cubic functions */

	/* should get 3 solutions: 2, 6 & 2.5   */
	SolveCubic(a1, b1, c1, d1, &solutions, x);
	hashIncrResultData(x, sizeof(x));

	/* should get 1 solution: 2.5           */
	SolveCubic(a2, b2, c2, d2, &solutions, x);
	hashIncrResultData(x, sizeof(x));

	SolveCubic(a3, b3, c3, d3, &solutions, x);
	hashIncrResultData(x, sizeof(x));

	SolveCubic(a4, b4, c4, d4, &solutions, x);
	hashIncrResultData(x, sizeof(x));

	/* Now solve some random equations */

//	for(a1=1;a1<10;a1++) {
//		for(b1=10;b1>0;b1--) {
//			for(c1=5;c1<15;c1+=0.5) {
//				for(d1=-1;d1>-11;d1-=1) {

//		much less calculations!
	for(a1=1;a1<10;a1++) {
		for(b1=10;b1>0;b1-=2) {
			for(c1=5;c1<15;c1+=1.5) {
				for(d1=-1;d1>-11;d1-=3) {
					SolveCubic(a1, b1, c1, d1, &solutions, x);
					hashIncrResultData(x, sizeof(x));
				}
			}
		}
	}

	/* perform some integer square roots */
	for (i = 0; i < 1001; i++)
	{
		usqrt(i, &q);
		hashIncrResultData(&q, sizeof(q));
	}

	for (l = 0x3fed0169L; l < 0x3fed4169L; l++)
	{
		usqrt(l, &q);
		hashIncrResultData(&q, sizeof(q));
	}


	// ??? its in the original code...
	usqrt(l, &q);
	hashIncrResultData(&q, sizeof(q));

	/* convert some rads to degrees */
	for (X = 0.0; X <= 360.0; X += 1.0) {
		y = deg2rad(X);
		hashIncrResultData(&y, sizeof(y));
	}

	for (X = 0.0; X <= (2 * PI + 1e-6); X += (PI / 180)) {
		y = rad2deg(X);
		hashIncrResultData(&y, sizeof(y));
	}

	*FI_STOP = 1;

	*CPU_DONE = 1;

	return 0;
}
