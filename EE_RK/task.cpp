#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <limits.h>
#include "linalg.h"
#include "linalg.cpp"

#define _USE_MATH_DEFINES

#define FINPNAME "vals.txt"
#define FOUTNAME "res.txt"
#define Q_PARAM 7
#define DEFQSTEPS 10


/************************************************
 * l*u''(t)+2*lambda*u'(t)+g*sin(u(t)) = b*sin(omega*t)
 *
 * equals (substitution v=u'(t)):
 * 
 * // v' = -2*lambda/l-g/l*sin(u(t))+b/l*sin(omega*t)
 * \\ u' = v
************************************************/



vector<double> f(double t, vector<double> y, double l, double lambda, double g, double b, 
									double omega);


int main(int argc, char ** argv)
	{
	try
		{
		const char* finpname = FINPNAME;
		if(argc > 2) finpname = argv[2];
		FILE* finp = fopen(finpname, "r");
		if(!finp)
			{
			printf("Can't open %s\n", finpname);
			throw;
			}
		double l = 0, lambda = 0, g = 0, b = 0, omega = 0, u_o = 0, v_o;
		int retval = fscanf(finp, "%lf %lf %lf %lf %lf %lf %lf", &l, &lambda, &g, &b, &omega, 
					&u_o, &v_o);
		if(retval != Q_PARAM)
			{
			printf("Wrong input file format\n");
			throw;
			}
		
		printf("========\n Parameters of the task:\n l = %lf\n lambda = %lf\n g = %lf\n \
b = %lf\n omega = %lf\n u_0 = %lf\n v_0 = %lf\n========\n",  l, lambda, g, b, omega, u_o, v_o);
		
		long steps = DEFQSTEPS;
		
		if (argc > 1)
			{
			char* endptr;
			steps = strtol(argv[1], &endptr, 10);
			if ((errno == ERANGE) && (steps == LONG_MAX || steps == LONG_MIN))
				{
				printf("Bad argument: out of range\n");
				throw;
				}
			if(endptr == argv[1])
				{
				printf("Bad argument: First argument must be number of steps. No digits found\n");
				throw;
				}
			if(steps < 1)
				{
				printf("Bad argument: number of steps can't be less than 1\n");
				throw;
				}
			}
		

		double T = 2 * M_PI * sqrt(l / g);

		double tau = 3 * T / steps;

		printf("Period T = %0.10f\n", T);
		printf("All time: 3 * T = %0.10f\n", 3 * T);
		printf("Numsteps = %ld\n", steps);
		printf("tau = %0.10f\n", tau);
		
		matrix<double> result (steps, 2);	

/************************************************
 * Euler method
************************************************/
		vector<double> y(2);
		y[0] = v_o;
		y[1] = u_o;

		double t = 0;

		printf("\nGauss:\n");
		
		int percent = 0;
		
		printf("00%%");	

		for(long i = 1; i < steps; i++)
			{
			percent = int(i * 1.0 / steps * 100);
			printf("\b\b\b\b%02d%%", percent);
			fflush(stdout);
			y = y + f(t, y, l, lambda, g, b, omega) * tau;
			t += tau;
			result[i][0] = y[1];
			}
		percent = int(steps * 1.0 / steps * 100);
		printf("\b\b\b\b%02d%%", percent);
/************************************************
 * Runge–Kutta method
************************************************/
		printf("\nRunge-Kutta:\n");

		y[0] = v_o;
		y[1] = u_o;

		t = 0;

		vector<double> k1 (2), k2 (2), k3 (2), k4(2);
		
		for(long i = 1; i < steps; i++)
			{
			percent = int(i * 1.0 / steps * 100);

			k1 = f(t, y, l, lambda, g, b, omega) * tau;
			k2 = f(t + tau / 2, y + k1 * 0.5, l, lambda, g, b, omega) * tau;
			k3 = f(t + tau / 2, y + k2 * 0.5, l, lambda, g, b, omega) * tau;
			k4 = f(t + tau, y + k3, l, lambda, g, b, omega) * tau;
			y = y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (1.0 / 6.0);
			t += tau;
			result[i][1] = y[1];
			printf("\b\b\b\b%02d%%", percent);
			fflush(stdout);
			}
		percent = int(steps * 1.0 / steps * 100);
		printf("\b\b\b\b%02d%%\n\n", percent);
		


		const char* foutname = FOUTNAME;
		if(argc > 3)foutname = argv[3];

		FILE* fout = fopen(foutname, "w");

		t = 0;
		for(long i = 1; i < steps; i++)
			{
			fprintf(fout, "t = %.10f EE: u = %.10f   RK: u = %.10f\n", t, result[i][0], 
									result[i][1]);
			t += tau;
			}
		
		printf("Result is written in %s\n", foutname);

		fclose(finp);
		fclose(fout);
		return 0;
		}
	catch(...)
		{
		std::cerr << "Some errors happened. Aborting\n";
		}
	}



		
//y[0] = v
//y[1] = u
vector<double> f(double t, vector<double> y, double l, double lambda, double g, double b, 
									double omega)
	{
	vector<double> res(2);
	res[0] = -2 * lambda / l * y[0] - g / l * sin(y[1])+ b / l * sin(omega * t);
	res[1] = y[0];
	return res;
	}

