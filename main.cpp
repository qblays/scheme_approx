#include "common.h"
#include "task_1.h"
#include "task_2.h"
#include "task_3.h"
#include "task_4.h"
#include <boost/program_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <iostream>
#include <optional>
namespace po = boost::program_options;
int main(int argc, char *argv[])
{

  //if (argc != 5)
  //{
    //fprintf(stderr, "Use: %s mu gamma M N\n", argv[0]);
    //return -1;
  //}
  po::options_description general("Allowed options");
  general.add_options()
	("help", "produce help message")
	("task_1_help", "help message for task 1")
	("task_2_help", "help message for task 2")
	("task_3_help", "help message for task 3")
	("task_4_help", "help message for task 4")
	("gen_tables", "generate latex tables")

	("task_1", "do task 1")
	("task_2", "do task 2")
	("task_3", "do task 3")
	("task_4", "do task 4")
	("compression", po::value<int>(), "set compression level")
	;
	po::options_description opts_t1 ("task 1 options");
	opts_t1.add_options()
		("M", po::value<int>(), "set number of steps in space")
		("N", po::value<int>(), "set number of steps in time")
		("mu", po::value<double>(), "set viscosity")
		("gamma", po::value<double>(), "set gamma")
		("orig", "use original system filling")
		
	;
	po::options_description opts_t2 ("task 2 options");
	opts_t2.add_options()
		("tau", po::value<double>(), "step in time")
		("h", po::value<double>(), "step in space")
		("eps", po::value<double>(), "R < eps")
		("squeeze_t", "shift time net")
		("squeeze_h", "shift space net")
		("task_number", po::value<int>(), "task 1 or 2, both if not specified")
	;
	po::options_description opts_t3 ("task 3 options");
	opts_t3.add_options()
		("K", po::value<double>(), "K")
	;
	po::options_description opts_t4 ("task 4 options");
	opts_t4.add_options()
		("v_0", po::value<double>(), "v on first space step")
		("rho_0", po::value<double>(), "rho on first space step")
		("T_0", po::value<double>(), "T_0")
	;
  
  po::variables_map vm;
  po::options_description all("Allowed option");
	all.add(general).add(opts_t1).add(opts_t3).add(opts_t2).add(opts_t4);
  po::store(po::parse_command_line(argc, argv, all), vm);
  po::notify(vm);    

  if (vm.count("help")) {
		std::cout << general << '\n';
		return 0;
  }
	if (vm.count ("task_1_help")) {
		std::cout << opts_t1 << '\n';
		return 0;
	}
	if (vm.count ("task_2_help")) {
		std::cout << opts_t2 << '\n';
		return 0;
	}
	if (vm.count ("task_3_help")) {
		std::cout << opts_t2 << '\n';
		std::cout << opts_t3 << '\n';
		return 0;
	}
	if (vm.count ("task_4_help")) {
		std::cout << opts_t2 << '\n';
		std::cout << opts_t4 << '\n';
		return 0;
	}
	auto presence_check = [&] (const std::vector<std::string>& variables) -> bool
	{
		bool ret = true;
		for (auto &var: variables) {
			if (!vm.count (var))
			{
				ret = false;
				fprintf (stderr, "%s is not defined\n", var.c_str());
			}
		}
		return ret;
	};
  bool use_original_system = vm.count("orig");
	if (vm.count ("task_1")) {
		data data;
		data.T = PARAM_T;
		data.X = PARAM_X;
		if (!presence_check ({ "M", "N", "mu", "gamma" }))
		{
			return 1;
		}
		data.M = vm["M"].as<int> ();
		data.N = vm["N"].as<int> ();
		data.gamma = vm["gamma"].as<double> ();
		data.mu = vm["mu"].as<double> ();
		data.tau = data.T / data.N;
		data.h = data.X / data.M;
		if (vm.count ("gen_tables"))
		{
			task_1_gen_tables (data, use_original_system);
		}
		else {
			std::cerr << data << '\n';
			task_1_route (data, use_original_system);
		}
		return 0;
	}
	else if (vm.count ("task_2")) {
		data_t2 d;
		d.X = PARAM_X;
		if (!presence_check ({ "mu", "gamma", "tau", "eps", "h"}))
		{
			return 1;
		}
		d.gamma = vm["gamma"].as<double> ();
		d.mu = vm["mu"].as<double> ();
		d.tau = vm["tau"].as<double> ();
		d.h = vm["h"].as<double> ();
		d.eps = vm["eps"].as<double> ();
		d.M = d.X / d.h;
		if (vm.count("N") || vm.count("M")) {
			std::cerr << "specify h tau directly instead of setting M, N\n";
			return 1;
		}
		if (vm.count ("squeeze_t")) {
			d.squeeze_t = vm["squeeze_t"].as<double> ();
		}
		if (vm.count ("squeeze_h")) {
			d.squeeze_h = vm["squeeze_h"].as<double> ();
		}
		std::cerr << d << '\n';
		std::optional<int> task_number;
		if (vm.count ("task_number"))
		{
			task_number = vm["task_number"].as<int> ();
			if (task_number != 1 && task_number != 2)
			{
				fprintf (stderr, "task_number can be 1 or 2\n");
				return 1;
			}
		}
		if (vm.count ("gen_tables"))
		{
			if (task_number)
			{
				task_2_gen_tables (d, use_original_system, task_number.value());
			}
			else {
				task_2_gen_tables (d, use_original_system, 1);
				task_2_gen_tables (d, use_original_system, 2);
			}
		}
		else {
			if (task_number)
				task_2_route (d, use_original_system, task_number.value(), {});
			else
			{
				task_2_route (d, use_original_system, 1, {});
				task_2_route (d, use_original_system, 2, {});
			};
		}
		return 0;
	}
	else if (vm.count ("task_3")) {
		data_t2 d;
		d.X = 1;
		if (!presence_check ({ "mu", "gamma", "tau", "eps", "h"}))
		{
			return 1;
		}
		d.gamma = vm["gamma"].as<double> ();
		d.mu = vm["mu"].as<double> ();
		d.tau = vm["tau"].as<double> ();
		d.h = vm["h"].as<double> ();
		d.eps = vm["eps"].as<double> ();
		d.M = d.X / d.h;
		if (vm.count("N") || vm.count("M")) {
			std::cerr << "specify h tau directly instead of setting M, N\n";
			return 1;
		}
		std::cerr << d << '\n';
		std::optional<int> task_number;
		if (vm.count ("task_number"))
		{
			task_number = vm["task_number"].as<int> ();
			if (task_number != 1 && task_number != 2)
			{
				fprintf (stderr, "task_number can be 1 or 2\n");
				return 1;
			}
		}
		if (vm.count ("gen_tables"))
		{
			if (task_number)
			{
				task_3_gen_tables (d, use_original_system, task_number.value());
			}
			else {
				task_3_gen_tables (d, use_original_system, 1);
				task_3_gen_tables (d, use_original_system, 2);
			}
		}
		else {
			if (!presence_check ({"K"}))
			{
				return 1;
			}
		  d.K = vm["K"].as<double> ();
			if (task_number)
				task_3_route (d, use_original_system, task_number.value());
			else
			{
				task_3_route (d, use_original_system, 1);
				task_3_route (d, use_original_system, 2);
			};
		}
		return 0;
		
	}
	else if (vm.count ("task_4")) {
		data_t2 d;
		d.X = 10;
		if (!presence_check ({ "mu", "gamma", "tau", "eps", "h"}))
		{
			return 1;
		}
		d.gamma = vm["gamma"].as<double> ();
		d.mu = vm["mu"].as<double> ();
		d.tau = vm["tau"].as<double> ();
		d.h = vm["h"].as<double> ();
		d.eps = vm["eps"].as<double> ();
		d.M = d.X / d.h;
		if (vm.count("N") || vm.count("M")) {
			std::cerr << "specify h tau directly instead of setting M, N\n";
			return 1;
		}
		if (vm.count("v_0"))
		{
			d.v_0 = vm["v_0"].as<double> ();
		}
		if (vm.count("rho_0"))
		{
			d.rho_0 = vm["rho_0"].as<double> ();
		}
		std::cerr << d << '\n';
		if (vm.count ("gen_tables"))
		{
			task_4_gen_tables(d, use_original_system);
		}
		else {
			task_4_route (d, use_original_system);
		}
		return 0;
		
	}
  double mu = atof(argv[1]);
  double gamma = atof(argv[2]);
  int M = atoi(argv[3]);
  double N = atof(argv[4]);

  if (mu <= 0. || gamma <= 0.)
  {
    fprintf(stderr, "Input param (mu, gamma) incorrect (mu anc C must be great 0)\n");
    return -1;
  }

  if (M <= 0)
  {
    fprintf(stderr, "Input param M incorrect (M must be great 0)\n");
    return -1;
  }

  if (N <= 0)
  {
    fprintf(stderr, "Input param tau incorrect (N must be great 0)\n");
    return -1;
  }

  data data;
  data.T = PARAM_T;
  data.X = PARAM_X;
  data.gamma = gamma;
  data.mu = mu;

  data.N = N;
  data.M = M;
  data.tau = data.T / N;
  data.h = data.X / M;
  printf("h=%.3e, tau=%.3e, gamma=%.3e, mu = %.3e\n", data.h, data.tau, data.gamma, data.mu);

  arr_t V_curr(M + 1);
  arr_t H_curr(M);
  arr_t V_next(M + 1);
  arr_t H_next(M);


  scheme_1(data, V_curr, H_curr, V_next, H_next);
  //scheme_2 (data, V_curr, H_curr, V_next, H_next);
  double t = N * data.tau;
  double residual_v = 0.;
  double residual_h = 0.;

  for (int m = 0; m < M; m++)
  {

    double x = m * data.h;

    double app_g = H_curr[m];
    double ori_g = rho(t, x + data.h / 2.);
    double val_g = fabs(app_g - ori_g);

    double app_v = V_curr[m];
    double ori_v = u(t, x);
    double val_v = fabs(app_v - ori_v);

    residual_v = std::max(residual_v, val_v);
    residual_h = std::max(residual_h, val_g);
    //if (m < 30)
  }

  printf("residual_V = %e ; residual_H = %e\n", residual_v, residual_h);
}


