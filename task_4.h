#include "task_2.h"
result_t2 scheme_4(const data_t2 &d, arr_t &V_curr, arr_t &H_curr, arr_t &V_next,
		arr_t &H_next, bool use_original_system, int task_number, std::vector<double> timesteps_to_save, bool osci, double K,
		bool to_save_ts_info);
result_t2 task_4_route (data_t2 d, bool use_original_system);
double rho_4 (double t, double x, double);
double v_4 (double t, double x, double);
double find_diff_norm (const arr_t& a, const arr_t& b);
void dump_tables_task_4 (const arr_t& mu_vec, const arr_t& rho_vec,const arr_t& v_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec); 


void task_4_gen_tables (const data_t2& d, bool use_original_system);
void dump_summary_table_task_4 (const arr_t& mu_vec, const arr_t& rho_vec, const arr_t& v_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec);
void dump_graph_slice_V_task_4 (const arr_t& mu_vec, const arr_t& rho_vec,const arr_t& v_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec);
