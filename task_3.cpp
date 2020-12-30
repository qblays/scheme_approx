#include "task_3.h"
#include "ctpl/ctpl.h"
#include <cstdio>
result_t2 task_3_route (data_t2 d, bool use_original_system, int task_number)
{
	arr_t V_curr(d.M + 1);
	arr_t H_curr(d.M + 1);
	arr_t V_next(d.M + 1);
	arr_t H_next(d.M + 1);
	std::cout << "starting task 3" << std::endl;
	auto time_1 = std::chrono::high_resolution_clock::now ();
	auto result = scheme_2(d, V_curr, H_curr, V_next, H_next, use_original_system, task_number, {}, true, d.K, false);
	auto time_2 = std::chrono::high_resolution_clock::now ();

	printf("K = %.2lf, N = %lu, N tau = %lf, R = %.8e, mass_err = %.8e\n, task_number = %i,"
			, d.K
			, result.R_and_m_vec.size()
			, result.R_and_m_vec.size() * d.tau, result.R_and_m_vec.back()[0]
			, result.R_and_m_vec.back()[1]
			, task_number);
	std::cerr << "Elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>( time_2 - time_1 ).count()
		<< " ms" << '\n';
	result.elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( time_2 - time_1 ).count();
	return result;
}
void task_3_gen_tables (const data_t2& d, bool use_original_system, int task_number)
{
	fprintf (stderr, "generating tables for task 3(%i)\n", task_number);
	ctpl::thread_pool p(std::thread::hardware_concurrency ());
	fprintf (stderr, "using %i hardware threads\n", std::thread::hardware_concurrency());
	std::vector<data_t2> var_params;
	std::vector<double> mu_vec{0.1, 0.01, 1e-3};
	std::vector<double> K_vec{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	for (auto mu: mu_vec)
	{
		for (auto K: K_vec) 
		{
			data_t2 task_data = d;
			task_data.mu = mu;
			task_data.K = K;
			task_data.M = task_data.X / task_data.h;
			var_params.push_back (task_data);
		}
	}
	std::vector<std::future<result_t2>> results_f (var_params.size());
	std::vector<result_t2> results (var_params.size());
	std::unordered_set<int> visited;
	for (int i = 0; i < results_f.size(); i++)
	{
    if (var_params[i].mu <= mu_vec.back()) {
			results_f[i] = p.push ([&var_params, i, use_original_system, task_number](int) -> result_t2
					{
					return task_3_route(var_params[i], use_original_system, task_number);
					});
			visited.insert (i);
    }
	}
	for (int i = 0; i < results_f.size() && !visited.count (i); i++) 
	{
		results_f[i] = p.push ([&var_params, i, use_original_system, task_number](int) -> result_t2
				{
				return task_3_route(var_params[i], use_original_system, task_number);
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
	dump_tables_task_3 (mu_vec, K_vec, results, var_params, task_number);
}
void dump_tables_task_3 (const arr_t& mu_vec, const arr_t& K_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number) 
{
  dump_summary_table_task_3 (mu_vec, K_vec, res_vec, var_params_vec, task_number);
}
void dump_summary_table_task_3 (const arr_t& mu_vec, const arr_t& K_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec, int task_number)
{
  char buf[1000];
	sprintf (buf, "tables/task_3_%i_summary.tex", task_number);
	fprintf (stderr, "%s\n", buf);
	auto filename = buf;
	std::ofstream ostream;
	ostream.open(filename);
	std::string output;
  auto minibuf = std::string (1 + mu_vec.size (), 'c') ;
	printf ("printf strange minibuf \"%s\"\n", minibuf.c_str ());
	
	// make header
	sprintf(buf, "\\begin{NiceTabular}{%s}[hvlines]\n", minibuf.c_str());
	minibuf.clear();
  output+=buf;
	sprintf(buf, "\\diagbox{K}{$\\mu$}");
	output += buf;
	for (auto mu : mu_vec)
	{
		sprintf(buf, "& %.e ", mu);
		output += buf;
	}
	output += "\\\\ \n";

	int j = 0;
	int n = K_vec.size();
	for (double K : K_vec)
	{
		sprintf (buf, "\\hspace*{3mm}%i ", (int) K);
		output += buf;
		int i = 0;
    for (double mu : mu_vec)
		{
			sprintf(buf, "& %.4lf ", res_vec[i * n + j].R_and_m_vec.size() * var_params_vec[i*n+j].tau);
			output += buf;
			i++;
		}
	  output += "\\\\ \n";
		j++;
	}
	sprintf(buf, "\\end{NiceTabular}\n");
	output+=buf;
	ostream << output;
	ostream.close();
}
