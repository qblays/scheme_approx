#pragma once
#include "common.h"
#include "task_2.h"
result_t2 task_3_route (data_t2 d, bool use_original_system, int task_number);
void dump_summary_table_task_3 (const arr_t& mu_vec, const arr_t& K_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec, int task_number);
void task_3_gen_tables (const data_t2& d, bool use_original_system, int task_number);
void dump_tables_task_3 (const arr_t& mu_vec, const arr_t& K_vec,
		const std::vector<result_t2>& res_vec, const std::vector<data_t2>& var_params_vec,
		int task_number);
