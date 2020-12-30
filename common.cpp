#include "common.h"
void fill_system_for_V(arr_t &a, arr_t &b, arr_t &c, arr_t &V_next, arr_t &V_curr, arr_t &H_next, arr_t &H_curr, int n, const data &d, double (*f) (double, double, const data&))
{
  int M = d.M;
  double tau = d.tau;
  double h = d.h;
  double gamma = d.gamma;
  double mu = d.mu;

  a[0] = 0;
  b[0] = 1.;
  c[0] = 0.;
  V_next[0] = 0;
  for (int m = 1; m < M; m++)
  {
    double Hs = H_next[m] + H_next[m - 1];
    if (fabs (Hs) < 1e-12) {
      a[m] = 0;
      b[m] = 1;
      c[m]=0;
      V_next[m] = 0;
      continue;
    }
    V_next[m] = V_curr[m] / tau - gamma / (gamma - 1) * (pow(H_next[m], gamma - 1) - pow(H_next[m - 1], gamma - 1)) / h + f(n * tau + tau, h * m, d);
    if (V_curr[m] >= 0.)
    {
      a[m] = -V_curr[m] / h - 2 * mu / (h * h * Hs);
      b[m] = 1. / tau + V_curr[m] / h + 4 * mu / (h * h * Hs);
      c[m] = -2 * mu / (h * h * Hs);
    }
    else
    {
      a[m] = -2 * mu / (h * h * Hs);
      b[m] = 1. / tau - V_curr[m] / h + 4 * mu / (h * h * Hs);
      c[m] = V_curr[m] / h - 2 * mu / (h * h * Hs);
    }
  }
  a[M] = 0;
  b[M] = 1.;
  c[M] = 0.;
  V_next[M] = 0;
}
void fill_system_for_V_4 (arr_t &a, arr_t &b, arr_t &c, arr_t &V_next, arr_t &V_curr, arr_t &H_next, arr_t &H_curr,
		int n, const data &d, double (*f) (double, double, const data&))
{
  int M = d.M;
  double tau = d.tau;
  double h = d.h;
  double gamma = d.gamma;
  double mu = d.mu;
	double v_0 = d.v_0;

  a[0] = 0;
  b[0] = 1.;
  c[0] = 0.;
  V_next[0] = v_0;
  for (int m = 1; m < M; m++)
  {
    double Hs = H_next[m] + H_next[m - 1];
    if (fabs (Hs) < 1e-12) {
      a[m] = 0;
      b[m] = 1;
      c[m]=0;
      V_next[m] = 0;
      continue;
    }
    V_next[m] = V_curr[m] / tau - gamma / (gamma - 1) * (pow(H_next[m], gamma - 1) - pow(H_next[m - 1], gamma - 1)) / h + f(n * tau + tau, h * m, d);
    if (V_curr[m] >= 0.)
    {
      a[m] = -V_curr[m] / h - 2 * mu / (h * h * Hs);
      b[m] = 1. / tau + V_curr[m] / h + 4 * mu / (h * h * Hs);
      c[m] = -2 * mu / (h * h * Hs);
    }
    else
    {
      a[m] = -2 * mu / (h * h * Hs);
      b[m] = 1. / tau - V_curr[m] / h + 4 * mu / (h * h * Hs);
      c[m] = V_curr[m] / h - 2 * mu / (h * h * Hs);
    }
  }
  a[M] = -1;
  b[M] = 1.;
  c[M] = 0.;
  V_next[M] = 0;
}

void fill_system_for_V_orig(arr_t &a, arr_t &b, arr_t &c, arr_t &V_next, arr_t &V_curr, arr_t &H_next, arr_t &H_curr, int n, const data &datax, double (*f) (double, double, const data&))
{
  int N = datax.N;
  int M = datax.M;
  double tau = datax.tau;
  double h = datax.h;
  double gamma = datax.gamma;
  double mu = datax.mu;

  a[0] = 0;
  b[0] = 1.;
  c[0] = 0.;
  V_next[0] = 0;
  for (int m = 1; m < M; m++)
  {
    double Hbb = m == 1? H_next[0]: H_next[m-2];
    double Hf = m== M-1? H_next[M-1] : H_next[m+1];
		Hf = H_next[m+1];
    double Hs = H_next[m] + H_next[m - 1];
    if (fabs (Hs) < 1e-12) {
      a[m] = 0;
      b[m] = 1;
      c[m]=0;
      V_next[m] = 0;
      continue;
    }
    V_next[m] = Hs * V_curr[m] / (2. *tau) - gamma * Hs / (gamma - 1) * (pow(H_next[m], gamma - 1) - pow(H_next[m - 1], gamma - 1)) / (2. *h) + f(n * tau  + tau, h * m, datax) * Hs / 2.;
    if (V_curr[m] >= 0.)
    {
      a[m] = -H_next[m-1] * V_curr[m] / (2.*h) - Hbb * V_curr[m-1] / (2. * h) - mu/ (h*h);
      b[m] = Hs / (2. * tau) + H_next[m] * V_curr[m+1]/(2. * h) + H_next[m-1]*V_curr[m] / (2. * h) + 2*mu/(h*h);
      c[m] = -mu / (h*h);
    }
    else
    {
      a[m] = -mu / (h*h);
      b[m] = Hs / (2. * tau) - H_next[m] * V_curr[m]/(2. * h) - H_next[m-1]*V_curr[m-1] / (2. * h) + 2*mu/(h*h);
      c[m] = Hf * V_curr[m+1] / (2 *h) + H_next[m] * V_curr[m] / (2.*h) - mu/(h*h);
    }
  }
  a[M] = 0;
  b[M] = 1.;
  c[M] = 0.;
  V_next[M] = 0;
}
void init_H (arr_t &H_curr, const data & data, double K, double (*h_init)(double t, double x, double K)) {
  auto M = H_curr.size();
  for (size_t i = 0; i< M; i++) {
    double x = i * data.h + data.h / 2.; 
    H_curr[i] = h_init (0, x, K);
  }
}
void init_V (arr_t &V_curr, const data & data, double K, double (*v_init)(double t, double x, double K)) {
  auto M = V_curr.size();
  for (size_t i = 0; i< M; i++) {
    double x = i * data.h; 
    V_curr[i] = v_init (0, x, K);
  }
}
void tridiagonal_solve(const arr_t &a, const arr_t &b, arr_t &c, arr_t &d, int n)
{
  n--; // since we start from x0 (not x1)
  c[0] /= b[0];
  d[0] /= b[0];
  for (int i = 1; i < n; i++)
  {
    c[i] /= b[i] - a[i] * c[i - 1];
    d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
  }
  d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);
  for (int i = n; i-- > 0;)
  {
    d[i] -= c[i] * d[i + 1];
  }
}
