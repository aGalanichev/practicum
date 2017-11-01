#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double alpha;
double PI = 3.14159265358979;
double E = 2.718281828;
double t0, t1;   //отрезок интегрирования
double A = 0, B = 0, C = 0;
double H[3];

double a[6] = { 0, 1./4, 3./8, 12./13, 1, 1./2 };

double b[6][5] = { { 0, 0, 0, 0, 0 },
{ 1. / 4, 0, 0, 0, 0 },
{ 3./32, 9./32, 0, 0, 0 },
{ 1932./2197, 7296./2197, 0, 0, 0 },
{ 439./216, -8, 3680./513, -845./4104, 0 },
{ -8./27, 2, -3544./2565, 1859./4104, -7./20 } };

double p[6] = { 25./216, 0, 1408./2565, 2197./4104, -1./5, 0 };

double f_x(double t, double x, double y, double z, double Px, double L1, double L2);
double f_y(double t, double x, double y, double z, double Px, double L1, double L2);
double f_z(double t, double x, double y, double z, double Px, double L1, double L2);
double f_Px(double t, double x, double y, double z, double Px, double L1, double L2);

double f_B0(double t, double x, double y, double z, double Px, double L1, double L2);


double f_Lmax(double t, double x, double y, double z, double Px, double L1, double L2);


void next(double *x, double *y, double *z, double *Px, double L1, double L2, double *B0, double t, double h);
void Runge(double *X, double *Y, double *Z, double *PX, double L1, double L2, double h, double tol, int print, double t1);
void RungeBack(double *X, double *Y, double *Z, double *PX, double L1, double L2, double h, double tol);
void Neuton(double a, double b, double c, double delta);

double F_a(double a, double b, double c);
double F_b(double a, double b, double c);
double F_c(double a, double b, double c);


void Matrix(double J[3][3], double F[3]);

void function(double x, double y, double z, double Px, double L1, double L2);

int main()
{
	
	double x = 0, y = 0, z = 0, Px = 2, L1 = 2, L2 = 2;
	double ALPHA = 10;	t0 = 0; t1 = 1;

	for(alpha = 0; alpha <= ALPHA; alpha += .1)
	{
		printf("alpha=%.1f   Ньютон:\n", alpha);
		Neuton(Px, L1, L2, .001);
		Px = A, L1 = B, L2 = C;

		function(x, y, z, Px, L1, L2);
	}

	

	printf("Нажмите клавишу\n");
	getchar();
	return 0;
}


void function(double X, double Y, double Z, double PX, double L1, double L2)
{
	double x = X, y = Y, z = Z, Px = PX;
	double x8  = X, y8  = Y, z8  = Z, Px8  = PX;
	double x10 = X, y10 = Y, z10 = Z, Px10 = PX;
	double x12 = X, y12 = Y, z12 = Z, Px12 = PX;
	double tol8 = pow(10, -8), tol10 = pow(10, -10), tol12 = pow(10, -12);
	double R_x, R_Px;	

	printf("alpha=%.1f\n", alpha);
	printf("\nНайденные значения:\n");
	printf("x(0)=%f    y(0)=%f    z(0)=%f    Px(0)=%f   L1=%f   L2=%f\n\n", x, y, z, Px, L1, L2);
	
	printf("\nЗначения, вычисленные с помощью Рунге-Кутта:\n");
	Runge( &x, &y, &z, &Px, L1, L2, .00000001, tol12, 1, t1);
	printf("x(PI)=%f   y(PI)=%f   z(PI)=%f   Px(PI)=%f\n\n", x, y, z, Px);

	RungeBack( &x, &y, &z, &Px, L1, L2, -.00000001, tol12);
	printf("Проход назад:\n");
	printf("x(0)=%f    y(0)=%f    z(0)=%f    Px(0)=%f\n\n", x, y, z, Px);

	
	Runge( &x8 , &y8 , &z8 , &Px8 , L1, L2, .00000001, tol8 , 2, t1);
	Runge( &x10, &y10, &z10, &Px10, L1, L2, .00000001, tol10, 2, t1);
	Runge( &x12, &y12, &z12, &Px12, L1, L2, .00000001, tol12, 2, t1);
	R_x = fabs((x8 - x10)/(x10 - x12));	R_Px = fabs((Px8 - Px10)/(Px10 - Px12));
	printf("Числа Рунге:\n");
	printf("R_x = %f   R_Px=%f\n\n", R_x, R_Px);
}



double f_x(double t, double x, double y, double z, double Px, double L1, double L2)  { return Px*exp(alpha*t); }
double f_y(double t, double x, double y, double z, double Px, double L1, double L2)  { return x; }
double f_z(double t, double x, double y, double z, double Px, double L1, double L2)  { return (x*x*x*t)/(1 + alpha*t*t); }
double f_Px(double t, double x, double y, double z, double Px, double L1, double L2) { return L1 + L2*(x*x*t)/(1 + alpha*t*t); }

double f_B0(double t, double x, double y, double z, double Px, double L1, double L2) { return Px*Px*exp(alpha*t); }

// считает |макс. соб. зн.|
double f_Lmax(double t, double x, double y, double z, double Px, double L1, double L2) 
{ return fabs( exp(alpha*t)/2 + 3*L2*x*t/(1 + alpha*t*t) ); }

double F_a(double a, double b, double c)
{
	double x = 0, y = 0, z = 0, Px = a, L1 = b, L2 = c;
	Runge( &x, &y, &z, &Px, L1, L2, .0000001, .00000001, 0, t1);
	return x - 1;
}

double F_b(double a, double b, double c)
{
	double x = 0, y = 0, z = 0, Px = a, L1 = b, L2 = c;
	Runge( &x, &y, &z, &Px, L1, L2, .0000001, .00000001, 0, t1);
	return y;
}

double F_c(double a, double b, double c)
{
	double x = 0, y = 0, z = 0, Px = a, L1 = b, L2 = c;
	Runge( &x, &y, &z, &Px, L1, L2, .0000001, .00000001, 0, t1);
	return z;
}

void Neuton(double a, double b, double c, double delta)
{
	double eps = 0.00001;
	double F[3];
	double J[3][3];
	int k;
 
	F[0] = F_a(a, b, c); F[1] = F_b(a, b, c); F[2] = F_c(a, b, c);
	do
	{
		J[0][0] = (F_a(a + eps, b, c) - F_a(a - eps, b, c)) / (2 * eps);
		J[0][1] = (F_a(a, b + eps, c) - F_a(a, b - eps, c)) / (2 * eps);
		J[0][2] = (F_a(a, b, c + eps) - F_a(a, b, c - eps)) / (2 * eps);
		J[1][0] = (F_b(a + eps, b, c) - F_b(a - eps, b, c)) / (2 * eps);
		J[1][1] = (F_b(a, b + eps, c) - F_b(a, b - eps, c)) / (2 * eps);
		J[1][2] = (F_b(a, b, c + eps) - F_b(a, b, c - eps)) / (2 * eps);
		J[2][0] = (F_c(a + eps, b, c) - F_c(a - eps, b, c)) / (2 * eps);
		J[2][1] = (F_c(a, b + eps, c) - F_c(a, b - eps, c)) / (2 * eps);
		J[2][2] = (F_c(a, b, c + eps) - F_c(a, b, c - eps)) / (2 * eps);
		
		Matrix( J, F); k = 1;
		//Модификация Исаева-Сонина
		/*do
		{
			G[0] = F_a(a + H[0]/k, b + H[1]/k, c + H[2]/k); 
			G[1] = F_b(a + H[0]/k, b + H[1]/k, c + H[2]/k); 
			G[2] = F_c(a + H[0]/k, b + H[1]/k, c + H[2]/k);
			k++;
			printf("	k=%d	a=%f  b=%f  c=%f      F_a=%f  F_b=%f  F_c=%f\n", k, a, b, c, G[0], G[1], G[2]);
		}
		while(fabs(F[0]) + fabs(F[1]) + fabs(F[2]) < fabs(G[0]) + fabs(G[1]) + fabs(G[2]));*/
		a += H[0]/k; b += H[1]/k; c += H[2]/k;
		F[0] = F_a(a, b, c); F[1] = F_b(a, b, c); F[2] = F_c(a, b, c);
		printf("a=%f  b=%f  c=%f      F_a=%f  F_b=%f  F_c=%f\n", a, b, c, F[0], F[1], F[2]); 

	} while( sqrt(pow(F[0], 2) + pow(F[1], 2) + pow(F[2], 2) ) > delta);
	A = a, B = b, C = c; 
}


void next(double *x, double *y, double *z, double *Px, double L1, double L2, double *B0, double t, double h)
{
	double K_x[6] = { 0, 0, 0, 0, 0, 0 };
	double K_y[6] = { 0, 0, 0, 0, 0, 0 };
	double K_z[6] = { 0, 0, 0, 0, 0, 0 };
	double K_Px[6] = { 0, 0, 0, 0, 0, 0 };
	double K_B0[6] = { 0, 0, 0, 0, 0, 0 };
	double SUM_X = 0, SUM_Y = 0, SUM_Z = 0, SUM_Px = 0, SUM_B0 = 0;
	double sum_x = 0, sum_y = 0, sum_z = 0, sum_Px = 0, sum_B0 = 0;
	int i, j;

	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < i; j++)
		{
			sum_x  += b[i][j] * K_x[j];
			sum_y  += b[i][j] * K_y[j];
			sum_z  += b[i][j] * K_z[j];
			sum_Px += b[i][j] * K_Px[j];
			sum_B0 += b[i][j] * K_B0[j];
		}
		K_x[i]  = h * f_x( t + a[i] * h, *x + sum_x, *y + sum_y, *z + sum_z, *Px + sum_Px, L1, L2);
		SUM_X  += p[i] * K_x[i];
		K_y[i]  = h * f_y( t + a[i] * h, *x + sum_x, *y + sum_y, *z + sum_z, *Px + sum_Px, L1, L2);
		SUM_Y  += p[i] * K_y[i];
		K_z[i]  = h * f_z( t + a[i] * h, *x + sum_x, *y + sum_y, *z + sum_z, *Px + sum_Px, L1, L2);
		SUM_Z  += p[i] * K_z[i];
		K_Px[i] = h * f_Px(t + a[i] * h, *x + sum_x, *y + sum_y, *z + sum_z, *Px + sum_Px, L1, L2);
		SUM_Px += p[i] * K_Px[i]; 
		K_B0[i] = h * f_B0(t + a[i] * h, *x + sum_x, *y + sum_y, *z + sum_z, *Px + sum_Px, L1, L2); 
		SUM_B0 += p[i] * K_B0[i];
	}

	*x  += SUM_X;
	*y  += SUM_Y;
	*z  += SUM_Z;
	*Px += SUM_Px;
	*B0 += SUM_B0;
}


void Runge(double *X, double *Y, double *Z, double *PX, double L1, double L2, double h, double tol, int print, double t1)
{
	double x = *X, y = *Y, z = *Z, Px = *PX;
	double X2, Y2, Z2, PX2;
	double Wx, Wy, Wz, WPx;
	double B0 = 0, B02, WB0;
	double t;
	double err, h0;
	double fac = 0.8, facmin = 0.95, facmax;
	double L, global = 0;
	char filename[100];

	if (print == 1) sprintf(filename, "track %.1f", alpha);
	FILE *file; if(print == 1) file = fopen( filename, "w");
	
	for (t = t0; t < t1; t += 2 * h)
	{
		facmax = 1.2;
		X2 = x, Y2 = y, Z2 = z, PX2 = Px, B02 = B0;
		Wx = x, Wy = y, Wz = z, WPx = Px, WB0 = B0;
		next(&X2, &Y2, &Z2, &PX2, L1, L2, &B02, t, h);
		next(&X2, &Y2, &Z2, &PX2, L1, L2, &B02, t + h, h);
		next(&Wx, &Wy, &Wz, &WPx, L1, L2, &WB0, t, 2 * h);
		err = sqrt(pow(X2 - Wx, 2) + pow(Y2 - Wy, 2) + pow(Z2 - Wz, 2) + pow(PX2 - WPx, 2)) / (pow(2, 5) - 1);
		do
		{
			h0 = fac * pow((tol / err), 1. / (5 + 1));
			h = h * fmin(facmax, fmax(facmin, h0));
			X2 = x, Y2 = y, Z2 = z, PX2 = Px, B02 = B0;
			Wx = x, Wy = y, Wz = z, WPx = Px, WB0 = B0;
			next(&X2, &Y2, &Z2, &PX2, L1, L2, &B02, t, h);
			next(&X2, &Y2, &Z2, &PX2, L1, L2, &B02, t + h, h);
			next(&Wx, &Wy, &Wz, &WPx, L1, L2, &WB0, t, 2 * h);
			err = sqrt(pow(X2 - Wx, 2) + pow(Y2 - Wy, 2) + pow(Z2 - Wz, 2) + pow(PX2 - WPx, 2)) / (pow(2, 5) - 1);
			facmax = 1;
			
		} while (err > tol);

		if(print)
		{ 
			L = f_Lmax(t, x, y, z, Px, L1, L2);
			global = err + global * pow(E, 2*h*L);
			if(print == 1) fprintf(file,"%f\t%f\t%f\t%f\n", t, x, Px, B0);
		}

		x = X2; y = Y2, z = Z2, Px = PX2, B0 = B02;
		if(isnan(x) || isnan(y) || isnan(z) || isnan(Px) || isnan(B0) ) { printf("Рунге сломался\n"); return; }
	}

	*X = x; *Y = y, *Z = z, *PX = Px;
	if(print) 
	{ 
		printf("B0=%f   global error= %.10f\n", B0, global); 
		if(print == 1) fclose(file); 
	}
}

void RungeBack(double *X, double *Y, double *Z, double *PX, double L1, double L2, double h, double tol)
{
	double x = *X, y = *Y, z = *Z, Px = *PX;
	double X2, Y2, Z2, PX2;
	double Wx, Wy, Wz, WPx;
	double B0 = 0, B02, WB0;
	double t;
	double err, h0;
	double fac = 0.8, facmin = 0.95, facmax;
	double L, global = 0;
	
	for (t = t1; t > t0; t += 2 * h)
	{
		facmax = 1.2;
		X2 = x, Y2 = y, Z2 = z, PX2 = Px, B02 = B0;
		Wx = x, Wy = y, Wz = z, WPx = Px, WB0 = B0;
		next(&X2, &Y2, &Z2, &PX2, L1, L2, &B02, t, h);
		next(&X2, &Y2, &Z2, &PX2, L1, L2, &B02, t + h, h);
		next(&Wx, &Wy, &Wz, &WPx, L1, L2, &WB0, t, 2 * h);
		err = sqrt(pow(X2 - Wx, 2) + pow(Y2 - Wy, 2) + pow(Z2 - Wz, 2) + pow(PX2 - WPx, 2)) / (pow(2, 5) - 1);
		do
		{
			h0 = fac * pow((tol / err), 1. / (5 + 1));
			h = h * fmin(facmax, fmax(facmin, h0));
			X2 = x, Y2 = y, Z2 = z, PX2 = Px, B02 = B0;
			Wx = x, Wy = y, Wz = z, WPx = Px, WB0 = B0;
			next(&X2, &Y2, &Z2, &PX2, L1, L2, &B02, t, h);
			next(&X2, &Y2, &Z2, &PX2, L1, L2, &B02, t + h, h);
			next(&Wx, &Wy, &Wz, &WPx, L1, L2, &WB0, t, 2 * h);
			err = sqrt(pow(X2 - Wx, 2) + pow(Y2 - Wy, 2) + pow(Z2 - Wz, 2) + pow(PX2 - WPx, 2)) / (pow(2, 5) - 1);
			facmax = 1;
			
		} while (err > tol);

		L = f_Lmax(t, x, y, z, Px, L1, L2);
		global = err + global * pow(E, 2*h*L);

		x = X2; y = Y2, z = Z2, Px = PX2, B0 = B02;
		if(isnan(x) || isnan(y) || isnan(z) || isnan(Px) || isnan(B0) ) { printf("Рунге сломался!!!\n"); return; }
	}

	*X = x; *Y = y, *Z = z, *PX = Px;
	printf("global error Back= %.10f\n", global);
}

// считает H= -J^(-1)*F
void Matrix(double J[3][3], double F[3])
{
	int i, j, n, m, q, p;
	double DET;
	double M[2][2], J_o[3][3];

	DET = J[0][0]*J[1][1]*J[2][2] + J[0][1]*J[1][2]*J[2][1] + J[1][0]*J[2][1]*J[0][2] - J[2][0]*J[1][1]*J[0][2] - J[0][1]*J[1][0]*J[2][2] - J[0][0]*J[1][2]*J[2][1];
		
		// Считаем обратную матрицу: J^(-1)	
		for( i = 0; i < 3; i++)
		{
			for( j = 0; j < 3; j++)
				{
					
					for( n = 0; n < 3; n++)
					{
						
						for( m = 0; m < 3; m++)
						{
							if(n == i || m == j) continue;
							q = n; p = m;
							if(n > i) q--;	if(m > j) p--;
							M[q][p] = J[n][m];
						}
						
					} 
					J_o[j][i] = pow(-1, i + j)*(M[0][0]*M[1][1] - M[1][0]*M[0][1]) / DET;
				}
		}
		
		// считаем H= -J^(-1)*F
		H[0] = 0; H[1] = 0; H[2] = 0;
		for( i = 0; i < 3; i++)
			for( j = 0; j < 3; j++)
			{
				H[i] -= J_o[i][j] * F[j]; 
			}		
}

