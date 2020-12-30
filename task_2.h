#pragma once
#include "common.h"

#include "task_1.h"
#include <math.h>
struct data_t2: public data
{
	// some const
	double eps;
};
inline std::ostream& operator<<(std::ostream& stream, const data_t2& d) {
	return stream << (const data&)d
		<< "eps: " << d.eps << '\n';
}

// answer on space step
struct v_h
{
	// velocity
	double v = 0;
	// density
	double h = 0;

  v_h& operator+=(const v_h& rhs)
	{
		v+=rhs.v;
		h+=rhs.h;
		return *this;
	}
  v_h& operator/=(double c)
	{
		v/=c;
		h/=c;
		return *this;
	}
};

struct result_on_timestep
{
	std::vector<v_h> v_h_spacestep;
	// space coord and ans
	std::pair<double, v_h> get (int i) const
	{
		return {spacestep_size * i, v_h_spacestep[i]};
	}
	int size () const
	{
		return v_h_spacestep.size ();
	}
	double spacestep_size;
};

struct result_on_some_timesteps
{
	//map timestep number to results
	std::vector<std::pair<int, result_on_timestep>> timesteps;
	//std::vector<std::pair<int, result_on_timestep>> timesteps_reduced;
	double timestep_size;


	int size () const
	{
		return timesteps.size();
	}
	// time and ans
	std::pair<double, result_on_timestep&> get (int i)
	{
		return {timesteps[i].first * timestep_size, timesteps[i].second};
	}
	const result_on_timestep& operator [] (int i) const
	{
		return timesteps[i].second;
	}
	void add (const arr_t& V, const arr_t& H, double spacestep_size, int timestep)
	{
    result_on_timestep r;
		r.spacestep_size = spacestep_size;
		r.v_h_spacestep.reserve(V.size());
    for (int i = 0; i < V.size(); i++)
		{
			r.v_h_spacestep.push_back({V[i], H[i]});
		}
		timesteps.push_back ({timestep, r});
	}
	result_on_some_timesteps erase_timesteps (int limit) const
	{
		result_on_some_timesteps ret;
		ret = *this;
		if (timesteps.size() <= limit)
			return ret;
		int old_size = timesteps.size();
    ret.timesteps.erase (ret.timesteps.begin () + limit, ret.timesteps.end ());
		fprintf (stderr, "ts from %i erased, size changed %i -> %d\n", limit, old_size, ret.size ());
		return ret;
	}
	result_on_some_timesteps reduce (int max_size = 1000) const
	{
		result_on_some_timesteps ret;
		if (timesteps.size() <= max_size) {
			ret.timesteps = timesteps;
			return ret;
		}
		double c = timesteps.size() / (double)max_size;
		int cc = ceil(c);
		int new_size = timesteps.size() / cc + (timesteps.size() % cc > 0);
		fprintf (stderr, "timesteps size is greater then %i, reducing it's size to %d by averaging\n", max_size, new_size);
		ret.timesteps.resize(new_size);
		ret.timestep_size = timestep_size * cc;
		for (auto &ts_info : ret.timesteps)
		{
			ts_info.second.v_h_spacestep.resize (timesteps[0].second.v_h_spacestep.size());
			ts_info.first = timesteps[0].first;
			ts_info.second.spacestep_size = timesteps[0].second.spacestep_size;
		}

		for (int i = 0; i < timesteps.size(); i++)
		{
			int j = i / cc;
			for (int k = 0; k < timesteps[i].second.size(); k++)
			{
        ret.timesteps[j].second.v_h_spacestep[k] += timesteps[i].second.v_h_spacestep[k];
			}
	    if ((i + 1) % cc == 0)
			{
				for (int k = 0; k < timesteps[i].second.size(); k++)
				{
					ret.timesteps[j].second.v_h_spacestep[k] /= cc;
				}
			}
			else if (i == timesteps.size() - 1)
			{
				for (int k = 0; k < timesteps[i].second.size(); k++)
				{
					ret.timesteps[j].second.v_h_spacestep[k] /= timesteps.size() % cc;
				}
			}
		}
		fprintf (stderr, "size reduced from %d to %d\n", size(), ret.size());
		return ret;
	}
	result_on_some_timesteps reduce_ss (int max_size = 500) const
	{
		result_on_some_timesteps ret;
		ret.timestep_size = timestep_size;
		int old_size = timesteps[0].second.size();
		if (old_size <= max_size) {
			ret.timesteps = timesteps;
			return ret;
		}
		double c = old_size / (double)max_size;
		int cc = ceil(c);
		int new_size = old_size / cc + (old_size % cc > 0);
		fprintf (stderr, "spacesteps size is greater then %i, reducing it's size to %d by averaging,cc=%i,ss=%e\n", max_size, new_size,
				cc, timesteps[0].second.spacestep_size * cc);
		ret.timesteps.resize(timesteps.size());
		for (auto &ts_info : ret.timesteps)
		{
			ts_info.second.v_h_spacestep.resize (new_size);
			ts_info.first = timesteps[0].first;
			ts_info.second.spacestep_size = timesteps[0].second.spacestep_size * cc;
		}

		for (int i = 0; i < timesteps.size(); i++)
		{
			for (int k = 0; k < timesteps[i].second.size(); k++)
			{
				int j = k / cc;
				ret.timesteps[i].second.v_h_spacestep[j] += timesteps[i].second.v_h_spacestep[k];
				if ((k + 1) % cc == 0)
				{
					ret.timesteps[i].second.v_h_spacestep[j] /= cc;
				}
				else if (k == timesteps[i].second.size() - 1)
				{
					ret.timesteps[i].second.v_h_spacestep[j] /= old_size % cc;
				}
			}
		}
		fprintf (stderr, "ss size reduced from %d to %d\n", old_size, ret.timesteps[0].second.size());
		return ret;
	}
};
struct result_t2 {
	std::vector<std::array<double, 2>> R_and_m_vec;
	mutable std::vector<std::array<double, 2>> reduced_R_and_m_vec;
	mutable double reduced_tau = 0.;

	// for some steps sace full answer
	result_on_some_timesteps full_answer;
  double elapsed; // in ms
	std::vector<std::array<double, 2>> reduce_R_and_m_vec_size (int wanted_size,
			double initial_tau, double &new_tau) const
	{
    if (R_and_m_vec.size() <= wanted_size)
		{
			new_tau = initial_tau;
			return R_and_m_vec;
		}
		double c = R_and_m_vec.size () / double (wanted_size);
		int ic = floor(c);
		new_tau = initial_tau * c;
		printf ("defined c = %le\n", c);
		std::vector<std::array<double, 2>> ret;
		ret.reserve(wanted_size);
		int j = 0;
    for (int i = 0; j < wanted_size; i+=ic, j++)
		{
			ret.push_back(R_and_m_vec[i]);
		}
		return ret;
	}
	void reduce_and_save_R_and_m_ves_size (double initial_tau) const
	{
		reduced_R_and_m_vec = 
			reduce_R_and_m_vec_size(6000, initial_tau, reduced_tau);
	}
};
double
rho_2_1 (double t, double x, double K);
double
rho_2_2 (double t, double x, double K);
double
rho_3_1 (double t, double x, double K);
double
rho_3_2 (double t, double x, double K);
double
v_2_1 (double t, double x, double K);
double
v_2_2 (double t, double x, double K);
double
v_3_1 (double t, double x, double K);
double
v_3_2 (double t, double x, double K);
double find_R (const arr_t &a, const arr_t &a_next, const arr_t &b, const arr_t &b_next, int n); 
result_t2 scheme_2(const data_t2 &d, arr_t &V_curr, arr_t &H_curr, arr_t &V_next,
		arr_t &H_next, bool use_original_system, int task_number, std::vector<double> timesteps_to_save, bool osci, double K,
		bool to_save_ts_info = false);
result_t2 task_2_route (data_t2 d, bool use_original_system, int task_number,
		const std::vector<double> ts_to_save);
void task_2_gen_tables (const data_t2& d, bool use_original_system, int task_number);
void dump_tables_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number, const std::vector<double>& ts_to_save);
void dump_summary_table_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec, int task_number);
void dump_summary_table_task_2_mass_error (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec, int task_number);
void dump_graph_R_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number);
void dump_graph_m_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number);
void dump_graph_slice_V_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number, const std::vector<double>& ts_to_graph);
void dump_graph_slice_H_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number, const std::vector<double>& ts_to_graph);
void dump_graph_colormap_task_2 (const arr_t& mu_vec, const arr_t& sq_h_vec, const arr_t& sq_t_vec, 
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number, std::optional<int> step_limit, char V_or_H);
