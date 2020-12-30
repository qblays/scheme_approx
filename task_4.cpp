#include "common.h"
#include "ctpl/ctpl.h"
#include "task_2.h"
#include "task_4.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <math.h>
#include <matplot/core/line_spec.h>
#include <matplot/freestanding/axes_functions.h>
#include <matplot/freestanding/plot.h>
#include <matplot/util/common.h>
#include <matplot/util/keywords.h>
#include <numeric>
#include <optional>
#include <queue>
result_t2 scheme_4(const data_t2 &d, arr_t &V_curr, arr_t &H_curr, arr_t &V_next,
		arr_t &H_next, bool use_original_system, int task_number, std::vector<double> timesteps_to_save, bool osci, double K,
		bool to_save_ts_info)
{
	result_t2 ret;
  int M = d.M;
  double tau = d.tau;
	ret.full_answer.timestep_size = tau;
  double h = d.h;
	init_H (H_curr, d, K, rho_4);
	init_V (V_curr, d, K, v_4);
  arr_t a(M + 1);
  arr_t b(M + 1);
  arr_t c(M + 1);
	double initial_mass = std::accumulate (H_curr.cbegin(), H_curr.cend(), 0);
	int k = 100;
	// queue of last k V vectors
	std::deque<arr_t> queue;
	auto enqueue = [k, &queue] (const arr_t& V)
	{
		queue.push_back(V);
		if (queue.size () > k)
		{
			queue.pop_front();
		}
	};
  auto to_save_timestep = [ tau, to_save_ts_info] (int ts)
	{
		return ts * tau < 1000 && to_save_ts_info;
	};
	int ts;
  for (ts = 0;;ts++)
  {
		// gas mass
		double mass = 0;
    /* === INIT SYSTEM FOR H === */
    // first equation
		{
			int m = 0;
			a[m] = 0;
			b[m] = 1;
			c[m] = 0;
			H_next[m] = d.rho_0;
			mass += H_curr[m];
		}
    for (int m = 1; m < M; m++)
    {
			mass += H_curr[m];
			a[m] = -(V_curr[m] + fabs(V_curr[m])) / (2. * h);
      //if (V_curr[m] >= 0)
      //{
        //b[m] = V_curr[m + 1] / h + 1. / tau;
        //c[m] = 0;
      //}
      //else
      //{
        //b[m] = -V_curr[m] / h + 1. / tau;
        //c[m] = V_curr[m + 1] / h;
      //}
			b[m] = 1. / tau + (V_curr[m+1] * (V_curr[m] >= 0) -V_curr[m] * (V_curr[m] < 0))/h;
			c[m] = (V_curr[m]< 0) * V_curr[m+1]/h;

      //double aa = -(V_curr[m] + fabs(V_curr[m])) / (2. * h);
			//double bb = 1. / tau + (V_curr[m+1] * (V_curr[m] >= 0) -V_curr[m] * (V_curr[m] < 0))/h;
			H_next[m] = H_curr[m] / tau;
			//if (m == 0) {
				//c[0] /= bb;
				//H_next[0] /= bb;
			//}
			//else if (m < M)
			//{
				//c[m] /= bb - aa * c[m - 1];
				//H_next[m] = (H_next[m] - aa * H_next[m - 1]) / (bb - aa * c[m - 1]);
			//}
			//else {
				//H_next[M] = (H_next[M] - aa * H_next[M - 1]) / (bb - aa * c[M - 1]);
			//}
    }
		{
			int m = M;
			a[m] = a[m-1];
			b[m] = b[m-1];
			c[m] = c[m-1];
			H_next[m] = H_next[m-1];
		}
		//for (int m = M-1; m-- > 0;)
		//{
			//H_next[m] -= c[m] * H_next[m + 1];
		//}

    /* === SOLVE SYSTEM FOR H === */
		tridiagonal_solve(a, b, c, H_next, M + 1);


    /* === INIT SYSTEM FOR V === */
		if (use_original_system)
			fill_system_for_V_orig(a, b, c, V_next, V_curr, H_next, H_curr, ts, d, [](double, double, const data&) {
					return 0.;
					});
		else
			fill_system_for_V_4(a, b, c, V_next, V_curr, H_next, H_curr, ts, d, [](double, double, const data&) {
					return 0.;
					});

    /* === SOLVE SYSTEM FOR V === */
    tridiagonal_solve(a, b, c, V_next, M + 1);
		bool to_break = 0;
		double R = 0;
		int failed_on = 0;
		if (queue.size() == k)
		{
			//for (const auto& V : queue)
			//{
				//R = find_diff_norm(V, V_next);
        //if (R > d.eps)
				//{
					//to_break = 0;
					//break;
				//}
				//to_break = 1;
				//failed_on++;
			//}
			auto &first = queue.front();
			R = find_diff_norm(first, V_next);
			if (R > d.eps)
			{
				to_break = 0;
			}
			else
			{
				to_break = 1;
				failed_on++;
			}
		}

		ret.R_and_m_vec.push_back ({R, (mass - initial_mass)/initial_mass});
		if (to_save_timestep(ts))
		{
      ret.full_answer.add(V_curr, H_curr, h, ts);
		}
		if (ts % 1000 == 0)
		{
			fprintf (stderr, "t = %f, R = %.3e\n", ts * tau, R);
		}
		if (to_break && ts > 1000) {
			ts++;
			break;
		}
		if (ts > 100000000) {
			ts++;
			break;
		}
		enqueue(V_next);
    std::swap(V_curr, V_next);
    std::swap(H_curr, H_next);
  }
	ret.full_answer.add(V_curr, H_curr, h, ts);
	return ret;
}
double find_diff_norm (const arr_t& a, const arr_t& b)
{
	double ret = 0;
	for (int i = 0; i < a.size(); i++)
	{
    ret = std::max (ret, fabs (a[i] - b[i]));
	}
	return ret;
}
result_t2 task_4_route (data_t2 d, bool use_original_system)
{
	arr_t V_curr(d.M + 1);
	arr_t H_curr(d.M + 1);
	arr_t V_next(d.M + 1);
	arr_t H_next(d.M + 1);
	//if (d.T_0 / d.tau < 1)
	//{
		//fprintf (stderr, "k = %.3e\n", d.T_0 / d.tau);
		//return result_t2{};
	//}
	fprintf (stderr, "k = %.3e\n", d.T_0 / d.tau);
	std::cerr << d;
	auto time_1 = std::chrono::high_resolution_clock::now ();
	auto result = scheme_4(d, V_curr, H_curr, V_next, H_next, use_original_system, 0, {}, false, 0, 0);
	auto time_2 = std::chrono::high_resolution_clock::now ();

	printf("N = %lu, N tau = %.4lf, R = %.8e, mass_err = %.8e\n,", result.R_and_m_vec.size()
			, result.R_and_m_vec.size() * d.tau, result.R_and_m_vec.back()[0]
			, result.R_and_m_vec.back()[1]
			);
	std::cerr << "Elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>( time_2 - time_1 ).count()
		<< " ms" << '\n';
	result.elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( time_2 - time_1 ).count();
	return result;
}
void task_4_gen_tables (const data_t2& d, bool use_original_system)
{
	fprintf (stderr, "generating tables for task 4\n");
	ctpl::thread_pool p(std::thread::hardware_concurrency ());
	fprintf (stderr, "using %i hardware threads\n", std::thread::hardware_concurrency());
	std::vector<data_t2> var_params;
	std::vector<double> mu_vec{0.1, 0.01, 0.001};
	//std::vector<double> mu_vec{0.1};
	std::vector<double> rho_vec{1, 2, 3, 4};
	std::vector<double> v_vec{1, 2, 3, 4};
	for (auto mu: mu_vec)
	{
		for (auto rho: rho_vec) 
		{
			for (auto v: v_vec)
			{
				data_t2 task_data = d;
				task_data.mu = mu;
				task_data.rho_0 = rho;
				task_data.v_0 = v;
				task_data.M = task_data.X / task_data.h;
				var_params.push_back (task_data);
			}
		}
	}
	std::vector<std::future<result_t2>> results_f (var_params.size());
	std::vector<result_t2> results (var_params.size());
	std::unordered_set<int> visited;
	//for (int i = 0; i < results_f.size(); i++)
	//{
    //if (var_params[i].mu <= mu_vec.back()) {
			//results_f[i] = p.push ([&var_params, i, use_original_system](int) -> result_t2
					//{
					//return task_4_route(var_params[i], use_original_system);
					//});
			//visited.insert (i);
    //}
	//}
	for (int i = 0; i < results_f.size() && !visited.count (i); i++) 
	{
		results_f[i] = p.push ([&var_params, i, use_original_system](int) -> result_t2
				{
				return task_4_route(var_params[i], use_original_system);
				});
	}
	for (int i = 0; i < results_f.size(); i++)
	{
		auto res = results_f[i].get ();
		results[i] = res;
		std::cerr << var_params[i] << "K: " << var_params[i].K << " R: " << res.R_and_m_vec.back()[0]
			<< " tau N: " << res.R_and_m_vec.size() * var_params[i].tau
			<< " N: " << res.R_and_m_vec.size()
			<< " mat error: " << res.R_and_m_vec.back()[1] << " elapsed: " << res.elapsed << '\n';
		std::cerr << "-------------------------------------------------------" << '\n';
	}
	dump_tables_task_4 (mu_vec, rho_vec, v_vec, results, var_params);
}
void dump_tables_task_4 (const arr_t& mu_vec, const arr_t& rho_vec,const arr_t& v_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec) 
{
  dump_summary_table_task_4 (mu_vec, rho_vec, v_vec, res_vec, var_params_vec);
	dump_graph_slice_V_task_4 (mu_vec, rho_vec, v_vec, res_vec, var_params_vec);
}
void dump_graph_slice_V_task_4 (const arr_t& mu_vec, const arr_t& rho_vec,const arr_t& v_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec)
{
	fprintf (stderr, "dump_graph_slice_V_task_4\n");
	std::vector <double> v1 {1, 2, 3};
	std::vector <double> v2 {2,4, 5};
	std::vector <double> xa {0.1, 0.2, 0.3};

	std::vector<double> V;
	std::vector<double> H;
	std::vector<double> x;
	int i = 0;
	x.reserve(res_vec[i].full_answer[0].size ());
  V.reserve(res_vec[i].full_answer[0].size ());
  H.reserve(res_vec[i].full_answer[0].size ());
	for (auto mu: mu_vec)
	{
		for (auto rho: rho_vec) 
		{
			for (auto v: v_vec)
			{
				double tau = var_params_vec[i].tau;
				int ts = res_vec[i].full_answer.timesteps.front().first;
				fprintf (stderr, "graphing time %.4e, full_answer size = %i\n", ts * tau, res_vec[i].full_answer.size());
				fprintf (stderr, "random v, h = %lf, %le\n", 
						res_vec[i].full_answer[0].get(30).second.v, res_vec[i].full_answer[0].get(30).second.v);
				auto f = plt::figure(true);
				auto ax = f->add_axes(true);
				char buf[1000];
				std::vector<std::string> legend;
				const auto& ts_info = res_vec[i].full_answer[0];
				fprintf (stderr, "ts_info.size = %i\n", ts_info.size());
				for (int ss = 0; ss < ts_info.size(); ss++)
				{
					x.push_back (ts_info.get(ss).first);
					V.push_back (ts_info.get(ss).second.v);
					H.push_back (ts_info.get(ss).second.h);
					fprintf (stderr, "xvh = {%lf, %lf, %lf}\n", ts_info.get(ss).first, ts_info.get(ss).second.v,
							ts_info.get(ss).second.h);
				}
				ax->plot (x, V)->line_width(1.5);
				ax->hold(true);
				ax->plot (x, H)->line_width(1.5);
				ax->xlabel ("x");

				double ts_time = tau * ts;
				sprintf (buf, "(V(t=%.2lf,x)", ts_time);
				legend.push_back (buf);
				sprintf (buf, "(H(t=%.2lf,x)", ts_time);
				legend.push_back (buf);
				V.clear ();
				x.clear ();
				H.clear ();
				ax->grid(matplot::on);
				sprintf (buf, "Срез для t=%.2lf, mu = %.1e, rho = %.1lf, v = %.1lf", ts_time, mu, rho, v);
				ax->title(buf);
				sprintf (buf, "teh_inc/graph_t4_HV_slice_%.3lf_rho_%.3lf_v_%.3lf.png", mu, rho, v);
				fprintf (stderr, "%s\n", buf);
				auto lgd = plt::legend(legend);
				//lgd->location(plt::legend::general_alignment::bottomleft);
				ax->legend(lgd);
				f->backend()->render_data();
				f->position(0, 0, 1000, 1000);
				f->save(buf);
				i++;
			}
		}
	}
}
void dump_summary_table_task_4 (const arr_t& mu_vec, const arr_t& rho_vec, const arr_t& v_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec)
{
  char buf[1000];
	int i = 0;
	for (auto mu: mu_vec)
	{
		sprintf (buf, "tables/task_4_%.3e_summary.tex", mu);
		std::string filename = buf;
		std::ofstream ostream;
		ostream.open(filename);
		std::string output;
		auto minibuf = std::string (1 + v_vec.size (), 'c') ;
		printf ("printf strange minibuf \"%s\"\n", minibuf.c_str ());

		// make header
		sprintf(buf, "\\begin{NiceTabular}{%s}[hvlines]\n", minibuf.c_str());
		minibuf.clear();
		output+=buf;
		sprintf(buf, "\\diagbox{$\\widetilde\\rho$}{$\\widetilde v$}");
		output += buf;
		for (auto v : v_vec)
		{
			sprintf(buf, "& %i ", (int)v);
			output += buf;
		}
		output += "\\\\ \n";
		for (auto rho: rho_vec) 
		{
			sprintf (buf, "\\hspace*{3mm}%i ", (int) rho);
			output += buf;
			for (auto v: v_vec)
			{
				sprintf(buf, "& %.4lf ", res_vec[i].R_and_m_vec.size() * var_params_vec[i].tau);
				output += buf;
				i++;
			}
			output += "\\\\ \n";
		}
		sprintf(buf, "\\end{NiceTabular}\n");
		output+=buf;
		ostream << output;
		ostream.close();
	}
}
double rho_4 (double t, double x, double)
{
	return 1;
}
double v_4 (double t, double x, double)
{
	return 0;
}
