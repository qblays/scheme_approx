#include <stdio.h>
#include <math.h>

#include <vector>

typedef struct
{
  // max time
  double T;
  // max space
  double X;
  double C;     // if p (rho) = C*rho
  double gamma; // if p (rho) = rho^gamma
  double mu;
  // number of steps in time
  int N;
  // number of steps in space
  int M;
  // size of time step
  double tau;
  // size of space step
  double h;
} data;

using arr_t = std::vector<double>;

double rho(double t, double x);
double rho_t(double t, double x);
double rho_x(double t, double x);

double g(double t, double x);
double g_t(double t, double x);
double g_x(double t, double x);

double u(double t, double x);
double u_t(double t, double x);
double u_x(double t, double x);
double u_xx(double t, double x);

double f0(double t, double x, const data &d);
double f(double t, double x, const data &d);

#define PARAM_T 1
#define PARAM_X 10
void fill_system_for_V(arr_t &a, arr_t &b, arr_t &c, arr_t &V_next, arr_t &V_curr, arr_t &H_next, arr_t &H_curr, int n, const  data &data);
double fabs(double d)
{
  return (d < 0) * (-d) + (d >= 0) * (d);
}

void tridiagonal_solve(arr_t &a, arr_t &b, arr_t &c, arr_t &d, int n);

void scheme(const data &data, arr_t &V_curr, arr_t &G_curr, arr_t &V_next, arr_t &G_next);
void scheme_1(const data &data, arr_t &V_curr, arr_t &G_curr, arr_t &V_next, arr_t &G_next);

int main(int argc, char *argv[])
{

  if (argc != 5)
  {
    fprintf(stderr, "Use: %s mu gamma M N\n", argv[0]);
    return -1;
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

  for (int m = 0; m <= M; m++)
  {
    double t = 0.;
    double x = m * data.h;
    V_curr[m] = u(t, x);
  }
  for (int m = 0; m < M; m++)
  {
    double t = 0.;
    double x = m * data.h + data.h / 2.;
    H_curr[m] = rho(t, x);
  }

  scheme_1(data, V_curr, H_curr, V_next, H_next);

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
    printf("%i: (%.3e, %.3e) ", m, val_v, val_g);
  }
  printf("\n");

  printf("residual_V = %e ; residual_H = %e\n", residual_v, residual_h);
}

/* ================================================================================================================== */
/* ================================================================================================================== */
/* ================================================================================================================== */

double
rho(double t, double x)
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

/* ================================================================================================================== */

double
g(double t, double x)
{
  return log(rho(t, x));
}

double
g_t(double /*t*/, double /*x*/)
{
  return 1.;
}

double
g_x(double t, double x)
{
  return rho_x(t, x) / rho(t, x);
}

/* ================================================================================================================== */

double
u(double t, double x)
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

/* ================================================================================================================== */

double
f0(double t, double x, const data & /*data*/)
{
  return rho_t(t, x) + u(t, x) * rho_x(t, x) + u_x(t, x) * rho(t, x);
  return g_t(t, x) + u(t, x) * g_x(t, x) + u_x(t, x);
}

double
f(double t, double x, const data &data)
{
  // return u_t(t, x) + u(t, x) * rho_t(t,x) / rho(t,x)
  //    + rho_x(t,x) * u(t,x) * u(t,x) / rho (t,x)
  //    + 2 * u_x(t,x)
  //    + data.gamma * pow(rho(t, x), data.gamma - 2) * rho_x(t, x)
  //    - data.mu * u_xx(t, x)/rho(t,x);
  return u_t(t, x) + u(t, x) * u_x(t, x) + data.gamma * pow(rho(t, x), data.gamma - 2) * rho_x(t, x) - data.mu * u_xx(t, x) / rho(t, x);
}

double
sigma(double F_next, double F_curr, double V)
{
  return (V < 0) * (F_next) + (V >= 0) * (F_curr);
}

double
back_diff_d(double f_curr, double f_prev, double h)
{
  return (f_curr - f_prev) / h;
}

/* ================================================================================================================== */
/* ================================================================================================================== */
/* ================================================================================================================== */

void tridiagonal_solve(arr_t &a, arr_t &b, arr_t &c, arr_t &d, int n)
{
  /*
  source: https://gist.github.com/nlw0/8edc1241bd05d5a9e5483bee763696a8
  // n is the number of unknowns
  |b0 c0 0 ||x0| |d0|
  |a1 b1 c1||x1|=|d1|
  |0  a2 b2||x2| |d2|
  1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->
      x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0
  2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
      from 1st it.: -| a1x0 + a1g0x1        = a1r0
                  -----------------------------
                        (b1 - a1g0)x1 + c1x2 = d1 - a1r0
      x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)
  3rd iteration:      | a2x1 + b2x2   = d2
      from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                     -----------------------
                     (b2 - a2g1)x2 = d2 - a2r2
      x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
  Finally we have a triangular matrix:
  |1  g0 0 ||x0| |r0|
  |0  1  g1||x1|=|r1|
  |0  0  1 ||x2| |r2|
  Condition: ||bi|| > ||ai|| + ||ci||
  in this version the c matrix reused instead of g
  and             the d matrix reused instead of r and x matrices to report results
  Written by Keivan Moradi, 2014
  */
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

/* ================================================================================================================== */
/* ================================================================================================================== */
/* ================================================================================================================== */

void scheme(const data &data, arr_t &V_curr, arr_t &H_curr, arr_t &V_next, arr_t &H_next)
{
  int N = data.N;
  int M = data.M;
  double tau = data.tau;
  double h = data.h;
  double gamma = data.gamma;

  double mu = data.mu;
  printf("gamma= %.3e, mu = %.3e\n", data.gamma, data.mu);
  arr_t a(M + 1);
  arr_t b(M + 1);
  arr_t c(M + 1);

  for (int n = 0; n < N; n++)
  {
    double t = n * tau;
    /* === INIT SYSTEM FOR H === */

    // first equation
    for (int m = 0; m < M; m++)
    {
      double x = m * h + h / 2;
      a[m] = -(V_curr[m] + fabs(V_curr[m])) / (2. * h);
      b[m] = 1. / tau + (V_curr[m + 1] + fabs(V_curr[m + 1]) - V_curr[m] + fabs(V_curr[m])) / (2. * h);
      c[m] = (V_curr[m + 1] - fabs(V_curr[m + 1])) / (2. * h);
      H_next[m] = H_curr[m] / tau + f0(t, x, data);
    }

    /* === SOLVE SYSTEM FOR H === */
    tridiagonal_solve(a, b, c, H_next, M);

    // для отладки (точное значение)
    for (int m = 0; m < M; m++)
    {
      double t = (n + 1) * tau;
      double x = m * h + h / 2;
      H_next[m] = rho(t, x);
    }

    /* === INIT SYSTEM FOR V === */

    // first equation
    b[0] = 1.;
    c[0] = 0.;
    V_next[0] = 0;

    for (int m = 1; m < M; m++)
    {
      double x = m * h;
      if (fabs(H_next[m - 1] + H_next[m]) > 1e-12)
      {
        if (m > 1)
        {
          a[m] = -(
                     (fabs(V_curr[m - 1]) + V_curr[m - 1]) * H_next[m - 2] +
                     (fabs(V_curr[m]) + V_curr[m]) * H_next[m - 1]) /
                     (4. * h) -
                 (mu / (h * h));
        }
        else
        {
          //printf("and this happened ..\n");
          a[m] = -((fabs(V_curr[m]) + V_curr[m]) * H_next[m - 1]) / (2. * h) - (mu / (h * h));
          // a[m] = 0;
        }
        b[m] = ((fabs(V_curr[m - 1]) - V_curr[m - 1] + fabs(V_curr[m]) + V_curr[m]) * H_next[m - 1] +
                (fabs(V_curr[m + 1]) + V_curr[m + 1] + fabs(V_curr[m]) - V_curr[m]) * H_next[m]) /
                   (4 * h) +
               (H_next[m - 1] + H_next[m]) / (2. * tau) + 2 * (mu / (h * h));
        if (m < M - 1)
        {
          c[m] = -(
                     (fabs(V_curr[m]) - V_curr[m]) * H_next[m] +
                     (fabs(V_curr[m + 1]) - V_curr[m + 1]) * H_next[m + 1]) /
                     (4. * h) -
                 (mu / (h * h));
        }
        else
        {
          c[m] = -(
                     (fabs(V_curr[m]) - V_curr[m]) * H_next[m]) /
                     (2. * h) -
                 (mu / (h * h));
        }
        V_next[m] = (H_next[m - 1] + H_next[m]) / 2 * f(t + tau, x, data) - gamma / (gamma - 1) * (sigma(H_next[m], H_curr[m], V_curr[m])) * (pow(H_next[m], gamma - 1.) - pow(H_next[m - 1], gamma - 1)) / (h) + (H_curr[m - 1] + H_curr[m]) * V_curr[m] / (2 * tau);
      }
      else
      {
        static int count = 0;
        a[m] = 0;
        b[m] = 1.;
        c[m] = 0.;
        V_next[m] = 0;
        printf("it happened %i times\n", count++);
      }
    }

    // last equation
    a[M] = 0.;
    b[M] = 1.;
    V_next[M] = 0;
    auto print_arr = [=](arr_t &a) {
      for (int i = 0; i < M + 1; i++)
      {
        printf("%i: %.3e ", i, a[i]);
      }
      printf("\n");
    };
    if (n == 500)
    {
      printf("a: \n");
      print_arr(a);
      printf("b: \n");
      print_arr(b);
      printf("c: \n");
      print_arr(c);
      printf("rhs: \n");
      print_arr(V_next);
    }

    /* === SOLVE SYSTEM FOR V === */
    tridiagonal_solve(a, b, c, V_next, M + 1);

    //для отладки (точное значение)
    for (int m = 0; m <= M; m++)
    {
      double t = (n + 1) * tau;
      double x = m * h + h / 2;
      V_next[m] /= rho(t, x);
      // if (m%250 == 0)
      //   printf ("V[%i] = %.3e -> %.3e\n", m, V_next[m], u(t,x));
      //V_next[m] = u(t, x);
    }

    std::swap(V_curr, V_next);
    std::swap(H_curr, H_next);
  }
}

void scheme_1(const data &data, arr_t &V_curr, arr_t &H_curr, arr_t &V_next, arr_t &H_next)
{
  int N = data.N;
  int M = data.M;
  double tau = data.tau;
  double h = data.h;
  double gamma = data.gamma;

  double mu = data.mu;
  printf("gamma= %.3e, mu = %.3e\n", data.gamma, data.mu);
  arr_t a(M + 1);
  arr_t b(M + 1);
  arr_t c(M + 1);

  for (int n = 0; n < N; n++)
  {
    double t = n * tau;
    /* === INIT SYSTEM FOR H === */

    // first equation
    for (int m = 0; m < M; m++)
    {
      double x = m * h + h / 2.;
      //a[m] = -(V_curr[m] + fabs(V_curr[m])) / (2. * h);
      //b[m] = 1. / tau + (V_curr[m + 1] + fabs(V_curr[m + 1]) - V_curr[m] + fabs(V_curr[m])) / (2. * h);
      //c[m] = (V_curr[m + 1] - fabs(V_curr[m + 1])) / (2. * h);
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

    /* === SOLVE SYSTEM FOR H === */
    tridiagonal_solve(a, b, c, H_next, M);

    // для отладки (точное значение)
    for (int m = 0; m < M; m++)
    {
      double t = (n + 1) * tau;
      double x = m * h + h / 2.;
      // H_next[m] = rho (t, x);
    }

    /* === INIT SYSTEM FOR V === */
    fill_system_for_V(a, b, c, V_next, V_curr, H_next, H_curr, n, data);
    // first equation
    //    a[0] = 0;
    //    b[0] = 1.;
    //    c[0] = 0.;
    //    V_next[0] = 0;
    //    if (fabs (t - 0.75) < 1e-13) {
    //      printf ("we are here\n");
    //    }
    //    for (int m = 1; m < M; m++)
    //    {
    //      double x = m * h;
    //      if (fabs(H_next[m - 1] + H_next[m]) > 1e-14)
    //      {
    //        double v = V_curr[m];
    //        if (fabs (v) < 1e-13){
    //          printf("pisdetz v = %.5e, m = %i, t=%.3e\n", v, m, t);
    //          v = -1e-8;
    //        }
    //        double av = fabs (v);
    //        double vf = V_curr[m+1];
    //        if (fabs (vf) < 1e-16){
    //          vf = -1e-8;
    //        }
    //        double avf  = fabs (vf);
    //        double vb = V_curr[m-1];
    //        double H = H_next[m];
    //        double Hf = H_next[m+1 == M? m: m+1];
    //        double Hb = H_next[m-1];
    //        double pf = Hf * (avf - vf);
    //        double pb = Hb*(av + v);
    //        b[m] = (vf/(4*avf * h) * (pf + Hb * (avf + vf))
    //                -v/(4*av  * h) * (H * (av - v) + pb));
    //        a[m] = -v /(2*av  * h) * (H * (av - v) + pb);
    //        c[m] =  vf/(2*avf * h) * (pf + H * (avf + vf));
    //
    //        b[m] += 2 *mu/(h*h);
    //        a[m] +=-1 *mu/(h*h);
    //        c[m] +=-1 *mu/(h*h);
    //
    //        b[m] += (Hb + H)/ (2*tau);
    //
    //        V_next[m] = - gamma/(gamma - 1)* (H+Hb)*(pow(H, gamma-1) - pow(Hb, gamma-1))/(2.*h)
    //                    + (H + Hb) * f(t, x, data) / 2.
    //                    + (H_curr[m] + H_curr[m-1])* v/(2*tau);
    //      }
    //      else
    //      {
    //        static int count = 0;
    //        a[m] = 0;
    //        b[m] = 1.;
    //        c[m] = 0.;
    //        V_next[m] = 0;
    //        printf ("it happened %i times\n", count++);
    //      }
    //    }
    //
    //    // last equation
    //    a[M] = 0.;
    //    b[M] = 1.;
    //    c[M] = 0;
    //    V_next[M] = 0;
    //    auto print_arr = [=] (arr_t &a) {
    //      for (int i = 0; i< M+1; i++) {
    //        printf ("%i: %.3e ", i, a[i]);
    //      }
    //      printf ("\n");
    //    };
    //    if (n == 0) {
    //      printf ("a: \n");
    //       print_arr (a);
    //       printf ("b: \n");
    //       print_arr (b);
    //       printf ("c: \n");
    //       print_arr (c);
    //       printf("rhs: \n");
    //       print_arr(V_next);
    //    }

    /* === SOLVE SYSTEM FOR V === */
    tridiagonal_solve(a, b, c, V_next, M + 1);

    //для отладки (точное значение)
    for (int m = 0; m < M; m++)
    {
      double t = (n + 1) * tau;
      double x = m * h;
      // if (m%250 == 0)
      //   printf ("V[%i] = %.3e -> %.3e\n", m, V_next[m], u(t,x));
      //V_next[m] = u(t, x);
      // if (fabs (V_next[m]) < 1e-13){
      //     printf("pisdetz v = %.5e, t=%.3e, x=%.3e\n", V_next[m], t, x);
      //   }
    }

    std::swap(V_curr, V_next);
    std::swap(H_curr, H_next);
  }
}

void fill_system_for_V(arr_t &a, arr_t &b, arr_t &c, arr_t &V_next, arr_t &V_curr, arr_t &H_next, arr_t &H_curr, int n, const data &data)
{
  int N = data.N;
  int M = data.M;
  double tau = data.tau;
  double h = data.h;
  double gamma = data.gamma;
  double mu = data.mu;

  a[0] = 0;
  b[0] = 1.;
  c[0] = 0.;
  V_next[0] = 0;
  for (int m = 1; m < M; m++)
  {
    double Hs = H_next[m] + H_next[m - 1];
    V_next[m] = V_curr[m] / tau - gamma / (gamma - 1) * (pow(H_next[m], gamma - 1) - pow(H_next[m - 1], gamma - 1)) / h + f(n * tau, h * m, data);
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