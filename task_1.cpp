#include "common.h"
#include "task_1.h"
#include "ctpl/ctpl.h"
#include <fstream>
#include <iomanip>
#include <ios>
void scheme_1(const data &data, arr_t &V_curr, arr_t &H_curr, arr_t &V_next, arr_t &H_next, bool use_original_system)
{
  int N = data.N;
  int M = data.M;
  double tau = data.tau;
  double h = data.h;

  init_H (H_curr, data, 0, rho);
  init_V (V_curr, data, 0, u);
  arr_t a(M + 1);
  arr_t b(M + 1);
  arr_t c(M + 1);

  for (int n = 0; n < N; n++)
  {
    double t = n * tau;
    for (int m = 0; m <= M; m++)
    {
      double x = m * h + h / 2.;
      a[m] = -(V_curr[m] + fabs(V_curr[m])) / (2. * h);
      if (V_curr[m] >= 0)
      {
        b[m] = V_curr[m + 1] / h + 1. / tau;
        c[m] = 0;
      }
      else
      {
        b[m] = -V_curr[m] / h + 1. / tau;
        c[m] = V_curr[m + 1] / h;
      }

      H_next[m] = H_curr[m] / tau + f0(t, x, data);
    }

    tridiagonal_solve(a, b, c, H_next, M + 1);

    //for (int m = 0; m < M; m++)
    //{
    //  double t = (n + 1) * tau;
    //  double x = m * h + h / 2.;
    //  H_next[m] = rho (t, x);
    //}

		if (use_original_system) {
			fill_system_for_V_orig(a, b, c, V_next, V_curr, H_next, H_curr, n, data, f_orig);
		}
		else {
			fill_system_for_V(a, b, c, V_next, V_curr, H_next, H_curr, n, data, f);

		} 

    tridiagonal_solve(a, b, c, V_next, M + 1);

    //for (int m = 0; m < M; m++)
    //{
    //  double t = (n + 1) * tau;
    //  double x = m * h;
    //  V_next[m] = u(t, x);
    //}

    std::swap(V_curr, V_next);
    std::swap(H_curr, H_next);
  }
}
// returns residual of v and h
std::array<double, 3> task_1_route (data d, bool use_original_system) {
	arr_t V_curr(d.M + 1);
	arr_t H_curr(d.M + 1);
	arr_t V_next(d.M + 1);
	arr_t H_next(d.M + 1);
	auto time_1 = std::chrono::high_resolution_clock::now ();
	scheme_1(d, V_curr, H_curr, V_next, H_next, use_original_system);
	auto time_2 = std::chrono::high_resolution_clock::now ();
	double t = d.N * d.tau;
	double residual_v = 0.;
	double residual_h = 0.;

	for (int m = 0; m < d.M; m++)
	{
		double x = m * d.h;
		double app_g = H_curr[m];
		double ori_g = rho(t, x + d.h / 2.);
		double val_g = fabs(app_g - ori_g);

		double app_v = V_curr[m];
		double ori_v = u(t, x);
		double val_v = fabs(app_v - ori_v);

		residual_v = std::max(residual_v, val_v);
		residual_h = std::max(residual_h, val_g);
	}

	printf("residual_V = %e ; residual_H = %e\n", residual_v, residual_h);
	std::cerr << "Elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>( time_2 - time_1 ).count()
		<< " ms" << '\n';
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( time_2 - time_1 ).count();
	return {residual_v, residual_h, static_cast<double>(elapsed)};
}
void task_1_gen_tables (const data& d, bool use_original_system)
{
	ctpl::thread_pool p(8);
	std::vector<data> var_params;
	std::vector<double> mu_vec{0.1, 0.01, 0.001};
	std::vector<double> m_vec{10, 100, 1000, 10000};
	std::vector<double> n_vec{100, 1000, 10000, 100000};
	for (auto mu: mu_vec)
	{
		for (auto m: m_vec) 
		{
			for (auto n: n_vec) 
			{
				data task_data = d;
				task_data.mu = mu;
				task_data.M = m;
				task_data.N = n;
				task_data.tau = task_data.T / task_data.N;
				task_data.h = task_data.X  / task_data.M;
				var_params.push_back (task_data);
			}
		}
	}
	std::vector<std::future<result_t1>> results_f (var_params.size());
	std::vector<result_t1> results (var_params.size());
	for (int i = 0; i < results_f.size(); i++) 
	{
		results_f[i] = p.push ([&var_params, i, use_original_system](int) -> result_t1
			{
				auto res = task_1_route(var_params[i], use_original_system);
				return {res[0], res[1], res[2]};
			});
	}
	for (int i = 0; i < results_f.size(); i++)
	{
		auto param = results_f[i].get ();
		results[i] = param;
		std::cerr << var_params[i] << param.v_resid << ' ' << param.h_resid << ' ' << param.elapsed << '\n';
		std::cerr << "--------------------------------------" << '\n';
	}
	dump_tables_task_1 (mu_vec, m_vec, n_vec, results);
}

 void dump_tables_task_1 (const arr_t& mu_vec, const arr_t& m_vec, const arr_t & n_vec, std::vector<result_t1>& results)
{
	std::ofstream ostream;
	ostream.open("tables/task_1.tex");
	ostream << std::scientific;
	int i = 0;
  for (auto mu: mu_vec) {
		ostream << "$$V,\\ \\mu = " << std::setprecision(1) << mu << ":$$\n";
    ostream << "\\begin{NiceTabular}{ccccc}[hvlines]\n";
		ostream << "\\diagbox{$\\tau$}{$h$}";
		for (auto n: n_vec) {
			ostream << std::setprecision(1)<< "& " << PARAM_X / n;
		}
		ostream << "\\\\ \n";
  	for (auto m: m_vec) {
			ostream << std::setprecision(1)<< PARAM_T / m;
  		for (auto n: n_vec) {
				ostream << std::setprecision(6)<< "& " << results[i].v_resid;
				i++;
  		}
			ostream << "\\\\ \n";
  	}
		ostream << "\\end{NiceTabular}\n";
  }
	 i = 0;
  for (auto mu: mu_vec) {
		ostream << "$$H,\\ \\mu = " << std::setprecision(1) << mu << ":$$\n";
    ostream << "\\begin{NiceTabular}{ccccc}[hvlines]\n";
		ostream << "\\diagbox{$\\tau$}{$h$}";
		for (auto n: n_vec) {
			ostream << std::setprecision(1)<< "& " << PARAM_X / n;
		}
		ostream << "\\\\ \n";
  	for (auto m: m_vec) {
			ostream << std::setprecision(1)<< PARAM_T / m;
  		for (auto n: n_vec) {
				ostream << std::setprecision(6)<< "& " << results[i].h_resid;
				i++;
  		}
			ostream << "\\\\ \n";
  	}
		ostream << "\\end{NiceTabular}\n";
  }
	ostream.close();
}
double
u(double t, double x, double)
{
  return cos(2. * M_PI * t) * sin(M_PI * x * x / 100.);
}


double
u_t(double t, double x)
{
  return -2. * M_PI * sin(2. * M_PI * t) * sin(M_PI * x * x / 100.);
}

double
u_x(double t, double x)
{
  return M_PI * x * cos(2. * M_PI * t) * cos(M_PI * x * x / 100.) / 50.;
}

double
u_xx(double t, double x)
{
  double ang = M_PI * x * x / 100.;
  double tmp = M_PI * x * x * sin(ang) - 50. * cos(ang);
  return -M_PI * cos(2. * M_PI * t) * tmp / 2500.;
}

double
rho(double t, double x, double)
{
  return exp(t) * (cos(M_PI * x / 10.) + 1.5);
}

double
rho_t(double t, double x)
{
  return rho(t, x);
}

double
rho_x(double t, double x)
{
  return -M_PI * exp(t) * sin(M_PI * x / 10.) / 10.;
}
double
f0(double t, double x, const data & /*data*/)
{
  return rho_t(t, x) + u(t, x) * rho_x(t, x) + u_x(t, x) * rho(t, x);
}
double
f_orig(double t, double x, const data &data)
{
  return u_t(t, x) + u(t, x) * rho_t(t,x) / rho(t,x)
     + rho_x(t,x) * u(t,x) * u(t,x) / rho (t,x)
     + 2 * u_x(t,x) * u(t,x)
     + data.gamma * pow(rho(t, x), data.gamma - 2) * rho_x(t, x)
     - data.mu * u_xx(t, x)/rho(t,x);
}
double
f(double t, double x, const data &data)
{
  return u_t(t, x) + u(t, x) * u_x(t, x) + data.gamma * pow(rho(t, x), data.gamma - 2) * rho_x(t, x) - data.mu * u_xx(t, x) / rho(t, x);
}
