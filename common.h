#pragma once
#include <ostream>
#include <chrono>
#include <vector>
#include <queue>
#include <functional>
#include <numeric>
#include <unordered_set>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <array>
#include <iomanip>
#include <matplot/matplot.h>
using arr_t = std::vector<double>;
using namespace std::string_literals;
namespace plt = matplot;
struct data
{
  // max time
  double T;
  // max space
  double X;
  // double C;     // if p (rho) = C*rho
  double gamma; // if p (rho) = rho^gamma
  double mu;
  // number of steps in time
  int N;
  // number of steps in space
  int M;
  // size of time step
  double tau;
  // size of space step
  double h;
	// shift of space net
	double squeeze_h = 1;
	// shift of time net
	double squeeze_t = 1;
	// K
	double K = 0;
	// v on 0 space step
	double v_0 = 1;
	// rho on 0 space step
	double rho_0 = 1;
	// T_0
	double T_0 = 0;
};
inline std::ostream& operator<<(std::ostream& stream, const data& d) {
	return stream << "T: "<< d.T << '\n'
		<<"X: "<<d.X<<'\n'
		<<"gamma: " << d.gamma << '\n'
		<< "mu: " << d.mu << '\n'
		<< "N: " << d.N << '\n'
		<< "M: " << d.M << '\n'
		<< "tau: " << d.tau << '\n'
		<< "h: " << d.h << '\n'
		<< "squeeze_h: " << d.squeeze_h << '\n'
		<< "squeeze_t: " << d.squeeze_t << '\n'
		<< "v_0: " << d.v_0 << '\n'
		<< "rho_0: " << d.rho_0 << '\n'
		<< "T_0: " << d.T_0 << '\n'
		<< "K: " << d.K << '\n';
}
#define PARAM_T 1
#define PARAM_X 10
inline bool cmp(double a, double b)
{
	return (a - 1e-12) < b && (a+1e-12 > b);
}
void init_H (arr_t &H_curr, const data & data, double K, double (*h_init)(double t, double x, double K));
void init_V (arr_t &V_curr, const data & data, double K, double (*v_init)(double t, double x, double K));
void fill_system_for_V(arr_t &a, arr_t &b, arr_t &c, arr_t &V_next, arr_t &V_curr, arr_t &H_next, arr_t &H_curr,
		int n, const  data &datax, double (*f) (double, double, const data&));
void fill_system_for_V_4 (arr_t &a, arr_t &b, arr_t &c, arr_t &V_next, arr_t &V_curr, arr_t &H_next, arr_t &H_curr,
		int n, const data &d, double (*f) (double, double, const data&));
void fill_system_for_V_orig(arr_t &a, arr_t &b, arr_t &c, arr_t &V_next, arr_t &V_curr, arr_t &H_next, arr_t &H_curr,
		int n, const data &datax, double (*f) (double, double, const data&));
void tridiagonal_solve(const arr_t &a, const arr_t &b, arr_t &c, arr_t &d, int n);
