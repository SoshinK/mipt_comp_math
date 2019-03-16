#include <iostream>
#include <cstdio>
#include <math.h>
#include "linalg.h"
#include "linalg.cpp"

/************************************************
 * -(g(x)y'(x))' + p(x)y(x) = f(x)
 *  y'(0) = a
 *  y'(1) = b
************************************************/


enum CONSTANTS
	{
	n = 20,
	POISON = -666,
	};

/*
double y(double x){return sin(x);} 
double dy(double x){return cos(x);}
double p(double x){return exp(x);}
double g(double x){return x * x + 1;}
double f(double x){return (x * x + 1) * sin(x) + exp(x) * sin(x) - 2 * x * cos(x);}
*/




double y(double x){return exp(x);} 
double dy(double x){return exp(x);}
double p(double x){return 1;}
double g(double x){return -1;}
double f(double x){return 0;}


int main()
	{
	double A = dy(0.0);
	double B = dy(1.0);
	printf("y'(0) = A = %f y'(1) = B = %f\n", A, B);	
	vector<double> a(n + 1);
	vector<double> b(n + 1);
	vector<double> c(n + 1);
	vector<double> d(n + 1);
	vector<double> t(n + 1);

	double tau = 1.0 / n;
	
	b[0] = 1.0 / tau;
	c[0] = 1.0 / tau;
	d[0] = A;

	for(int i = 0; i < n + 1; i++) t[i] = tau * double(i);
	
	for(int i = 1; i < n; i++)
		{
		a[i] = g(t[i] - tau / 2) / tau / tau;
		c[i] = g(t[i] + tau / 2) / tau / tau;
		b[i] = a[i] + c[i] - p(t[i]);
		d[i] = f(t[i]);
		}
	a[n] = -1.0 / tau;
	b[n] = -1.0 / tau;
	d[n] = B;

	std::cout << "a: " << a << '\n';
	std::cout << "b: " << b << '\n';
	std::cout << "c: " << c << '\n';
	std::cout << "d: " << d << '\n';
	
	vector<double> p(n + 1);
	vector<double> q(n + 1);
	
	p[1] = c[0] / b[0];
	q[1] = -d[0] / b[0];

	for(int i = 1; i < n; i++)
		{
		p[i + 1] = c[i] / (b[i] - a[i] * p[i]);
		q[i + 1] = (a[i] * q[i] - d[i]) / (b[i] - a[i] * p[i]);
		}

	std::cout << "!!!\n" << p << '\n' << q << '\n';

	vector<double> answer(n + 1);

	answer[n] = (d[n] - a[n] * q[n]) / (a[n] * p[n] - b[n]);

	for(int i = n; i > 0; i--) answer[i - 1] = p[i] * answer[i] + q[i];	

	for(int i = 0; i < n + 1; i++) printf("alg:%f real:%f\n", answer[i], y(t[i]));

	vector<double> diff(n + 1);
	for(int i = 0; i < n + 1; i++)diff[i] = fabs(answer[i] - y(t[i]));

	printf("max difference: %f\n", vmax(diff));

	return 0;
	}



