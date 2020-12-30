#pragma once
#include "common.h"
struct result_t1 {
	double v_resid;
	double h_resid;
	double elapsed;
};
double rho(double t, double x, double K = 0);
double rho_t(double t, double x);
double rho_x(double t, double x);
double u(double t, double x, double K = 0);
double u_t(double t, double x);
double u_x(double t, double x);
double u_xx(double t, double x);
double f0(double t, double x, const data &d);
double f(double t, double x, const data &d);
double f_orig(double t, double x, const data &d);

void scheme_1(const data &data, arr_t &V_curr, arr_t &G_curr, arr_t &V_next, arr_t &G_next, bool use_original_system = false);
std::array<double, 3> task_1_route (data d, bool use_original_system);
void task_1_gen_tables (const data& d, bool use_original_system);
void dump_tables_task_1 (const arr_t& mu_vec, const arr_t& m_vec, const arr_t & n_vec, std::vector<result_t1>& results);
