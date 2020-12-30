#include "common.h"
#include "ctpl/ctpl.h"
#include "task_2.h"
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
result_t2 scheme_2(const data_t2 &d, arr_t &V_curr, arr_t &H_curr, arr_t &V_next,
		arr_t &H_next, bool use_original_system, int task_number, std::vector<double> timesteps_to_save, bool osci, double K,
		bool to_save_all_ts_info)
{
	result_t2 ret;
  int M = d.M;
  double tau = d.tau;
	ret.full_answer.timestep_size = tau;
  double h = d.h;
	if (!osci)
	{
		if (task_number == 1)
		{
			init_H (H_curr, d, K, rho_2_1);
			init_V (V_curr, d, K, v_2_1);
		}
		else
		{
			init_H (H_curr, d, K, rho_2_2);
			init_V (V_curr, d, K, v_2_2);
		}
	}
	else {
		if (task_number == 1)
		{
			init_H (H_curr, d, K, rho_3_1);
			init_V (V_curr, d, K, v_3_1);
		}
		else
		{
			init_H (H_curr, d, K, rho_3_2);
			init_V (V_curr, d, K, v_3_2);
		}
	}
  arr_t a(M + 1);
  arr_t b(M + 1);
  arr_t c(M + 1);
	double initial_mass = std::accumulate (H_curr.cbegin(), H_curr.cend(), 0);
	
  auto to_save_timestep = [ tau, to_save_all_ts_info, &timesteps_to_save] (int ts)
	{
		double time = ts* tau;
		if (time < 10000 && to_save_all_ts_info)
		{
			return true;
		}
		if (std::binary_search(timesteps_to_save.begin(), timesteps_to_save.end(), time))
		{
			return true;
		}
		//if (std::find_if (timesteps_to_save.begin(), timesteps_to_save.end(), [time](double a)
					//{
						//return cmp(a, time);
					//})
				//!= timesteps_to_save.end())
		//{
			//return true;
		//}
		return false;
	};
	int ts;
  for (ts = 0;;ts++)
  {
		// gas mass
		double mass = 0;
    /* === INIT SYSTEM FOR H === */
    // first equation
    for (int m = 0; m <= M; m++)
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
			fill_system_for_V(a, b, c, V_next, V_curr, H_next, H_curr, ts, d, [](double, double, const data&) {
					return 0.;
					});

    /* === SOLVE SYSTEM FOR V === */
    tridiagonal_solve(a, b, c, V_next, M + 1);

    double R = find_R (V_curr, V_next, H_curr, H_next, M + 1);
		ret.R_and_m_vec.push_back ({R, (mass - initial_mass)/initial_mass});
		if (to_save_timestep(ts))
		{
      ret.full_answer.add(V_curr, H_curr, h, ts);
			//fprintf (stderr, "saving time %.4e\n", ts * tau);
		}
    std::swap(V_curr, V_next);
    std::swap(H_curr, H_next);
		if (ts % 10000 == 0)
		{
			printf ("n = %i, t = %lf, R = %.6e\n", ts, ts * tau, R);
		}
		if (R <= d.eps && ts > 1000) {
			break;
		}
		if (ts > 100000000) {
			break;
		}
  }
	ret.full_answer.add(V_curr, H_curr, h, ts);
	return ret;
}
// return vector of R
result_t2 task_2_route (data_t2 d, bool use_original_system, int task_number,
		const std::vector<double> ts_to_save) {
	arr_t V_curr(d.M + 1);
	arr_t H_curr(d.M + 1);
	arr_t V_next(d.M + 1);
	arr_t H_next(d.M + 1);
	auto time_1 = std::chrono::high_resolution_clock::now ();
	auto result = scheme_2(d, V_curr, H_curr, V_next, H_next, use_original_system, task_number,
			ts_to_save, false, 0, false);
	auto time_2 = std::chrono::high_resolution_clock::now ();

	printf("N = %lu, N tau = %.lf, R = %.8e, mass_err = %.8e\n, task_number = %i,"
			,	result.R_and_m_vec.size()
			, result.R_and_m_vec.size() * d.tau, result.R_and_m_vec.back()[0]
			, result.R_and_m_vec.back()[1]
			, task_number);
	std::cerr << "Elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>( time_2 - time_1 ).count()
		<< " ms" << '\n';
	result.elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( time_2 - time_1 ).count();
	return result;
}
void task_2_gen_tables (const data_t2& d, bool use_original_system, int task_number)
{
	fprintf (stderr, "generating tables for task 2(%i)\n", task_number);
	ctpl::thread_pool p(std::thread::hardware_concurrency ());
	std::vector<data_t2> var_params;
	//std::vector<double> mu_vec{0.1, 0.01};
	//std::vector<double> sq_t {1, 2};
	//std::vector<double> sq_h {1, 2};
	//std::vector<double> mu_vec{0.1, 0.01, 0.001};
	//std::vector<double> sq_t {1, 2};
	//std::vector<double> sq_h {1, 2};
	std::vector<double> mu_vec{0.1};
	std::vector<double> sq_t {1};
	std::vector<double> sq_h {1};
	std::vector<double> ts_to_save;
	for (int i = 15; i > -1; i--)
	{
		ts_to_save.push_back(i);
	}
	std::vector<double> ts_to_save_part2 {20., 30., 50., 100};
	ts_to_save.insert(ts_to_save.begin(), ts_to_save_part2.begin(), ts_to_save_part2.end());
	std::sort (ts_to_save.begin(), ts_to_save.end());
	for (auto mu: mu_vec)
	{
		for (auto squeeze_t: sq_t) 
		{
			for (auto squeeze_h: sq_h) 
			{
				data_t2 task_data = d;
				task_data.mu = mu;
				task_data.squeeze_t = squeeze_t;
				task_data.squeeze_h = squeeze_h;
				task_data.tau /= squeeze_t;
				task_data.h /= squeeze_h;
				task_data.M = task_data.X / task_data.h;
				var_params.push_back (task_data);
			}
		}
	}
	std::vector<std::future<result_t2>> results_f (var_params.size());
	std::vector<result_t2> results (var_params.size());
	std::unordered_set<int> visited;
	// first push tasks with lowest mu
	for (int i = 0; i < results_f.size(); i++)
	{
    if (var_params[i].mu <= mu_vec.back()) {
			results_f[i] = p.push (
					[&var_params, i, use_original_system, task_number, &ts_to_save](int) -> result_t2
					{
					return task_2_route(var_params[i], use_original_system, task_number, ts_to_save);
					});
			visited.insert (i);
    }
	}
	for (int i = 0; i < results_f.size() && !visited.count (i); i++) 
	{
		results_f[i] = p.push (
				[&var_params, i, use_original_system, task_number, &ts_to_save](int) -> result_t2
				{
				return task_2_route(var_params[i], use_original_system, task_number, ts_to_save);
				});
	}
	for (int i = 0; i < results_f.size(); i++)
	{
		auto res = results_f[i].get ();
		results[i] = res;
		std::cerr << var_params[i] << "R: " << res.R_and_m_vec.back()[0]
			<< " tau N: " << res.R_and_m_vec.size() * var_params[i].tau
			<< " N: " << res.R_and_m_vec.size()
			<< " mat error: " << res.R_and_m_vec.back()[1] << " elapsed: " << res.elapsed << '\n';
		std::cerr << "-------------------------------------------------------" << '\n';
	}
	dump_tables_task_2 (mu_vec, sq_h, sq_t, results, var_params, task_number, ts_to_save);
}

void dump_tables_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number, const std::vector<double>& ts_to_save) 
{
	dump_summary_table_task_2 (mu_vec, sq_h_vec, sq_t_vec, res_vec, var_params_vec, task_number);
	dump_summary_table_task_2_mass_error (mu_vec, sq_h_vec, sq_t_vec, res_vec,  
			var_params_vec, task_number);
	dump_graph_R_task_2 (mu_vec, sq_h_vec, sq_t_vec, res_vec, var_params_vec, task_number);
	dump_graph_m_task_2 (mu_vec, sq_h_vec, sq_t_vec, res_vec, var_params_vec, task_number);
	dump_graph_slice_V_task_2 (mu_vec, sq_h_vec, sq_t_vec, res_vec, var_params_vec, 
			task_number, ts_to_save);
	//dump_graph_colormap_task_2 (mu_vec, sq_h_vec, sq_t_vec, res_vec, var_params_vec, task_number, std::nullopt, 'H');
	//dump_graph_colormap_task_2 (mu_vec, sq_h_vec, sq_t_vec, res_vec, var_params_vec, task_number, std::nullopt, 'V');
	//dump_graph_colormap_task_2 (mu_vec, sq_h_vec, sq_t_vec, res_vec, var_params_vec, task_number, std::nullopt, 'H');
	//dump_graph_colormap_task_2 (mu_vec, sq_h_vec, sq_t_vec, res_vec, var_params_vec, task_number, std::nullopt, 'V');
}

void dump_summary_table_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec, 
		int task_number)
{
  char buf[1000];
	sprintf (buf, "tables/task_2_%i_summary_R.tex", task_number);
	fprintf (stderr, "%s\n", buf);
	auto filename = buf;
	std::ofstream ostream;
	ostream.open(filename);
	std::string output;
	// make header
	sprintf(buf, "\\begin{NiceTabular}{ccccccc}[hvlines]\n");
  output+=buf;
	sprintf(buf, "$\\mu$ & $\\Omega$ & $n = N_0/4$ & $n = N_0/2$ & $N = 3 N_0/4$ & $n=N_0$ & $N_0\\tau$ \\\\ \n");
  output+=buf;
	int i = 0;
	for (double mu : mu_vec)
	{
		for (double sq_t : sq_t_vec)
		{
			for (double sq_h : sq_h_vec)
			{
				double N_0 = res_vec[i].R_and_m_vec.size();
				auto get_R = [&](int n)
				{
					return res_vec[i].R_and_m_vec[n - 1][0];
				};
				sprintf(buf, 
						"%.1e & $\\Omega_{\\tau/%.1lf,h/%.1lf}$ & %.5e & %.5e & %.5e & %.5e & %.6lf\\\\ \n"
						, mu, sq_t, sq_h, get_R(N_0/4.), get_R(N_0/ 2.), get_R(0.75 * N_0), get_R(N_0) 
						, N_0* var_params_vec[i].tau);
				output+=buf;
				i ++;
			}
		}
		sprintf(buf, "\\hline \n");
		output += buf;
	}
	sprintf(buf, "\\end{NiceTabular}\n");
	output+=buf;
	ostream << output;
	ostream.close();
}
void dump_summary_table_task_2_mass_error (const arr_t& mu_vec, const arr_t& sq_h_vec,
		const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec, int task_number)
{
  char buf[1000];
	sprintf (buf, "tables/task_2_%i_summary_mass.tex", task_number);
	fprintf (stderr, "%s\n", buf);
	auto filename = buf;
	std::ofstream ostream;
	ostream.open(filename);
	std::string output;
	// make header
	sprintf(buf, "\\begin{NiceTabular}{ccccccc}[hvlines]\n");
  output+=buf;
	sprintf(buf, 
			"$\\mu$ & $\\Omega$ & $n = N_0/5$ & $n = 2N_0/5$ & $n = 3N_0/5$& $N = 4 N_0/5$ & $N =  N_0$ \\\\ \n");
  output+=buf;
	int i = 0;
	for (double mu : mu_vec)
	{
		for (double sq_t : sq_t_vec)
		{
			for (double sq_h : sq_h_vec)
			{
				double N_0 = res_vec[i].R_and_m_vec.size();
				auto get_m_err = [&](int n)
				{
					return res_vec[i].R_and_m_vec[n - 1][1];
				};
				sprintf(buf, "%.1e & $\\Omega_{\\tau/%.1lf,h/%.1lf}$ & %.5e & %.5e & %.5e & %.5e & %.5e\\\\ \n"
						, mu, sq_t, sq_h, get_m_err(0.2 * N_0), get_m_err(0.4 * N_0), get_m_err(0.6 * N_0)
						, get_m_err(0.8 * N_0), get_m_err(N_0));
				output+=buf;
				i ++;
			}
		}
		sprintf(buf, "\\hline \n");
		output += buf;
	}
	sprintf(buf, "\\end{NiceTabular}\n");
	output+=buf;
	ostream << output;
	ostream.close();
}
void dump_graph_R_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number)
{
	fprintf (stderr, "dump_graph_R_task_2\n");
	std::vector <double> v1 {1, 2, 3};
	std::vector <double> v2 {2,4, 5};
	std::vector <double> xa {0.1, 0.2, 0.3};
	int i = 0;
	std::vector<double> R;
	std::vector<double> x;
	x.reserve(res_vec[0].R_and_m_vec.size());
  R.reserve(res_vec[0].R_and_m_vec.size());
	for (double mu : mu_vec)
	{
		auto f = plt::figure(true);
		int w = 800;
		int h = 800;
		f->size(w, h);
		f->x_position(w);
		f->y_position(h);
		f->color({1, 0.94, 0.94, 0.94});
		printf ("w,h = %i, %i\n", f->width(), f->height());
		auto ax = f->add_axes(true);
		char buf[1000];
		std::vector<std::string> legend;
		int k = i;
		double max_y = 0;

		for (double sq_t : sq_t_vec)
		{
			for (double sq_h : sq_h_vec)
			{
				res_vec[k].reduce_and_save_R_and_m_ves_size(var_params_vec[k].tau);
				double loc_y_max = std::min (0.1, 
						3 * res_vec[k].reduced_R_and_m_vec[res_vec[k].reduced_R_and_m_vec.size()/2.][0]);
				max_y = std::max (max_y, loc_y_max);
				k++;
			}
		}
		for (double sq_t : sq_t_vec)
		{
			for (double sq_h : sq_h_vec)
			{
				for (int f = 0; f < res_vec[i].reduced_R_and_m_vec.size(); f++)
				{
					x.push_back (f*res_vec[i].reduced_tau);
					double y = std::min (res_vec[i].reduced_R_and_m_vec[f][0], max_y);
					R.push_back (y);
				}
				ax->hold(true);
				ax->plot (x, R)->line_width(1.5);
				ax->xlabel ("n tau");
				ax->ylabel ("R");
				
				sprintf (buf, "(Q(tau/%.0lf,h/%.0lf)", sq_t, sq_h);
				legend.push_back (buf);
				R.clear ();
				x.clear ();
				i++;
			}
		}
		//ax->colororder(newcolors);
		ax->grid(matplot::on);
		auto lgd = plt::legend(legend);
		ax->legend(lgd);
		sprintf (buf,"График R, mu = %.1e", mu);
		//f->title(buf);
		ax->title(buf);
		ax->title_enhanced(false);
		sprintf (buf, "teh_inc/graph_t2_%i_R_mu%.3lf.png", task_number, mu);
		f->save(buf);
		fprintf (stderr, "%s\n", buf);
		//sprintf (buf, 
				//"inkscape -w %i -h %i teh_inc/graph_t2_%i_R_mu%.3lf.svg -y 1 --export-filename teh_inc/graph_t2_%i_R_mu%.3lf.png"
				//,	w, h, task_number, mu, task_number, mu);
		//sprintf (buf, 
				//"svgexport teh_inc/graph_t2_%i_R_mu%.3lf.svg teh_inc/graph_t2_%i_R_mu%.3lf.png"
				//" \"svg{background:white;}\""
				//, task_number, mu, task_number, mu);
		//int ret = std::system (buf);
		//if (ret)
			//fprintf(stderr, "svgexport failed: %i\n", ret);
	}
}
void dump_graph_m_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number)
{
	fprintf (stderr, "dump_graph_m_task_2\n");
	std::vector <double> v1 {1, 2, 3};
	std::vector <double> v2 {2,4, 5};
	std::vector <double> xa {0.1, 0.2, 0.3};
	const double max_y = var_params_vec[0].eps * 10;

	int i = 0;
	std::vector<double> R;
	std::vector<double> x;
	x.reserve(res_vec[0].R_and_m_vec.size());
  R.reserve(res_vec[0].R_and_m_vec.size());
	for (double mu : mu_vec)
	{
		auto f = plt::figure(true);
		auto ax = f->add_axes(true);
		char buf[1000];
		std::vector<std::string> legend;
		for (double sq_t : sq_t_vec)
		{
			for (double sq_h : sq_h_vec)
			{
				res_vec[i].reduce_and_save_R_and_m_ves_size(var_params_vec[i].tau);
				for (int f = 0; f < res_vec[i].reduced_R_and_m_vec.size(); f++)
				{
					x.push_back (f*res_vec[i].reduced_tau);
					double y = std::min (res_vec[i].reduced_R_and_m_vec[f][1], max_y);
					R.push_back (y);
				}
				ax->hold(true);
				ax->plot (x, R)->line_width(1.5);
				ax->xlabel ("n tau");
				ax->ylabel ("Delta");
				
				sprintf (buf, "(Q(tau/%.0lf,h/%.0lf)", sq_t, sq_h);
				legend.push_back (buf);
				R.clear ();
				x.clear ();
				i++;
			}
		}
		ax->grid(matplot::on);
		sprintf (buf, "teh_inc/graph_t2_%i_m_mu%.3lf.png", task_number, mu);
		fprintf (stderr, "%s\n", buf);
		auto lgd = plt::legend(legend);
		lgd->location(plt::legend::general_alignment::bottomleft);
		ax->legend(lgd);
		f->backend()->render_data();
		f->position(0, 0, 1000, 1000);
		f->save(buf);
	}
}

void dump_graph_slice_V_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number, const std::vector<double>& ts_to_graph)
{
	fprintf (stderr, "dump_graph_slice_V_task_2\n");
	for (auto x: ts_to_graph){
		fprintf (stderr, "%lf ", x);
	}
	fprintf (stderr , "\n");
	auto to_save_timestep = [&ts_to_graph] (int ts, double tau, bool is_last_step)
	{
		if (is_last_step)
			return true;
		return std::binary_search (
				ts_to_graph.begin(), ts_to_graph.end(), tau * ts, [] (double a, double b)
				{
				return a < (b + 1e-12);
				});
	};
	std::vector <double> v1 {1, 2, 3};
	std::vector <double> v2 {2,4, 5};
	std::vector <double> xa {0.1, 0.2, 0.3};

	int i = 0;
	std::vector<double> V;
	std::vector<double> H;
	std::vector<double> x;
	x.reserve(res_vec[i].full_answer[0].size ());
  V.reserve(res_vec[i].full_answer[0].size ());
  H.reserve(res_vec[i].full_answer[0].size ());
	double mu = mu_vec.front();
		double sq_t = sq_t_vec.front();
		double sq_h = sq_h_vec.front();
		double tau = var_params_vec[i].tau;
    fprintf (stderr, "max n = %i, max t = %lf\n", res_vec[i].full_answer.size(), res_vec[i].full_answer.size() * tau);
		for (int ts = 0; ts < res_vec[i].full_answer.size(); ts++)
		{
			fprintf (stderr, "graphing time %.4e\n", ts * tau);
			auto f = plt::figure(true);
			auto ax = f->add_axes(true);
			char buf[1000];
			std::vector<std::string> legend;
			const auto& ts_info = res_vec[i].full_answer[ts];
			for (int ss = 0; ss < ts_info.size(); ss++)
			{
				x.push_back (ts_info.get(ss).first);
				V.push_back (ts_info.get(ss).second.v);
				H.push_back (ts_info.get(ss).second.h);
			}
			ax->plot (x, V)->line_width(1.5);
			ax->hold(true);
			ax->plot (x, H)->line_width(1.5);
			ax->xlabel ("x");

			double ts_time = tau * res_vec[i].full_answer.timesteps[ts].first;
			sprintf (buf, "(V(t=%.2lf,x)", ts_time);
			legend.push_back (buf);
			sprintf (buf, "(H(t=%.2lf,x)", ts_time);
			legend.push_back (buf);
			V.clear ();
			x.clear ();
			H.clear ();
			ax->grid(matplot::on);
			sprintf (buf, "Срез для t=%.2lf, mu = %.1e", ts_time, mu);
			ax->title(buf);
			sprintf (buf, "teh_inc/graph_t2_%i_HV_slice_%.3lf.png", task_number, ts_time);
		  fprintf (stderr, "%s\n", buf);
			auto lgd = plt::legend(legend);
			//lgd->location(plt::legend::general_alignment::bottomleft);
			ax->legend(lgd);
			f->backend()->render_data();
			f->position(0, 0, 1000, 1000);
			f->save(buf);
		}
}
void dump_graph_colormap_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number, std::optional<int> step_limit, char V_or_H)
{
	char filename[1000];
	if (V_or_H != 'H' && V_or_H != 'V')
	{
		fprintf (stderr, "exit: it is not V or H\n");
		return;
	}
	fprintf (stderr, "dump_graph_colormap_%c_task_2\n", V_or_H);
	auto gen_axe = [] (double a, double step, int size)
	{
		std::vector<double> ret;
		for (int i = 0; i < size; i++)
		{
      ret.push_back(a + i * step);
		}
		return ret;
	};

	int i = 0;
	for (double mu : mu_vec)
	{
		char buf[1000];
    sprintf (filename, "teh_inc/graph_t2_%i_%c_colormap_%.3lf%s.png", task_number, V_or_H, mu, step_limit? "_limit" : "");
		fprintf(stderr, "creating file %s\n", filename);
		double sq_t = sq_t_vec.front();
		double sq_h = sq_h_vec.front();
	  result_on_some_timesteps r = res_vec[i].full_answer.erase_timesteps (step_limit.value_or(10000000));
		int ts_number = r.size ();
		result_on_some_timesteps rt = r.reduce(500);
		result_on_some_timesteps rv = rt.reduce_ss();

		double tau = rv.timestep_size;
    fprintf (stderr, "max n = %i, max t = %lf\n", rv.size(), rv.size() * tau);
		auto f = plt::figure(true);
		auto ax = f->add_axes(true);
		std::vector<std::string> legend;
		std::vector<double> xx;
		std::vector<double> yy;
		yy = gen_axe(0, rv[0].spacestep_size, rv[0].size());
		printf ("tau = %.6e, size = %i\n", tau, rv.size());
		xx = gen_axe(0, tau, rv.size());
		auto [XX, YY] = plt::meshgrid(xx, yy);
		std::vector<std::vector<double>> ZZ;
		for (int ss = 0; ss < rv[0].size(); ss++)
		{
			std::vector<double> z;
			for (int ts = 0; ts < rv.size(); ts++)
			{
				if (V_or_H == 'V')
					z.push_back(rv[ts].get(ss).second.v);
				else
					z.push_back(rv[ts].get(ss).second.h);
			}
			ZZ.push_back(z);
		}
		printf ("x,y,z = %zu, %zu, %zu, ss size = %e, n of ss = %d\n", XX.size(), YY.size(), ZZ.size(),
				rv[0].spacestep_size, rv[0].size());
		ax->color_box(true);
		plt::line_spec ls;
		ls.line_width(20);
		ax->mesh(XX, YY, ZZ)->surface_in_2d(true).line_width(30);
		ax->minor_grid(false);
		ax->grid(false);
		ax->xlabel ("t");
		ax->ylabel ("x");
		sprintf (buf, "Динамика первых %d шагов %c при mu = %.1e", ts_number, V_or_H, mu);
		ax->title(buf);

		f->save(filename);
		i += sq_h_vec.size() * sq_t_vec.size();
	}
}
double find_R (const arr_t &a, const arr_t &a_next, const arr_t &b, const arr_t &b_next, int n) {
   double max = 0;
   for (int i = 0; i < n; i++) {
		 max = std::max (fabs (a_next[i]), max);
   }
   return max;
 }
double
rho_2_1 (double t, double x, double) {
  if (x >= 4.5 && x <= 5.5) {
    return 2;
  }
  return 1;
}
double
rho_2_2 (double t, double x, double) {
  return 1;
}
double
rho_3_1 (double t, double x, double K) {
  return 2 + sin(K * M_PI * x);
}
double
rho_3_2 (double t, double x, double K) {
  return 1;
}
double
v_2_1 (double t, double x, double) {
  return 0;
}
double
v_2_2 (double t, double x, double) {
  if (x >= 4.5 && x <= 5.5) {
    return 1;
  }
  return 0;
}
double
v_3_1 (double t, double x, double) {
  return 0;
}
double
v_3_2 (double t, double x, double K) {
  return sin(K * M_PI * x);
}
