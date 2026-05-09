/*!
 * \file
 * \brief Newton Search optimization algorithms - source file
 * \author Tony Ottosson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/optim/newton_search.h>
#include <itpp/base/specmat.h>
#include <itpp/stat/misc_stat.h>


namespace itpp
{


Newton_Search::Newton_Search()
{
  method = BFGS;

  initial_stepsize = 1.0;
  stop_epsilon_1 = 1e-4;
  stop_epsilon_2 = 1e-8;
  max_evaluations = 100;

  f = NULL;
  df_dx = NULL;

  no_feval = 0;
  init = false;
  finished = false;
  trace = false;
}

void Newton_Search::set_function(double(*function)(const vec&))
{
  // Add checks to see that function is OK???
  f = function;
}

void Newton_Search::set_gradient(vec(*gradient)(const vec&))
{
  // Add checks to see that function is OK???
  df_dx = gradient;
}

void Newton_Search::set_start_point(const vec &x, const mat &D)
{
  // check that parameters are valid???
  x_start = x;
  n = x.size();
  D_start = D;

  finished = false;
  init = true;
}

void Newton_Search::set_start_point(const vec &x)
{
  // check that parameters are valid???
  x_start = x;
  n = x.size();
  D_start = eye(n);

  finished = false;
  init = true;
}

bool Newton_Search::search()
{
  // Check parameters and function call ???
  // check that x_start is a valid point, not a NaN and that norm(x0) is not inf

  it_assert(f != NULL, "Newton_Search: Function pointer is not set");
  it_assert(df_dx != NULL, "Newton_Search: Gradient function pointer is not set");

  it_assert(init, "Newton_Search: Starting point is not set");


  F = f(x_start); // function initial value
  vec g = df_dx(x_start); // gradient initial value
  vec x = x_start;
  no_feval++;

  finished = false;

  // Initial inverse Hessian, D
  mat D = D_start;


  bool fst = true; // what is this???

  bool stop = false;

  // Finish initialization
  no_iter = 0;
  ng = max(abs(g)); // norm(g,inf)

  double Delta = initial_stepsize;
  nh = 0; // what is this???
  vec h;

  if (trace) { // prepare structures to store trace data
    x_values.set_size(max_evaluations);
    F_values.set_size(max_evaluations);
    ng_values.set_size(max_evaluations);
    Delta_values.set_size(max_evaluations);
  }

  Line_Search ls;
  ls.set_functions(f, df_dx);

  if (ng <= stop_epsilon_1)
    stop = true;
  else {
    h = zeros(n);
    nh = 0;
    ls.set_stop_values(0.05, 0.99);
    ls.set_max_iterations(5);
    ls.set_max_stepsize(2);
  }

  bool more = true; //???

  while (!stop && more) {
    vec h, w, y, v;
    double yh, yv, a;

    // Previous values
    vec xp = x, gp = g;
    // double Fp = F;           ### 2006-02-03 by ediap: Unused variable!
    double nx = norm(x);

    h = D * (-g);
    nh = norm(h);
    bool red = false;

    if (nh <= stop_epsilon_2*(stop_epsilon_2 + nx))  // stop criterion
      stop = true;
    else {
      if (fst || nh > Delta) {  // Scale to ||h|| = Delta
        h = (Delta / nh) * h;
        nh = Delta;
        fst = false;
        red = true;
      }
      //  Line search
      ls.set_start_point(x, F, g, h);
      more = ls.search(x, F, g);
      no_feval = no_feval + ls.get_no_function_evaluations();

      if (more == false) {  // something wrong in linesearch?
        x_end = x;
        return false;
      }
      else {
        if (ls.get_alpha() < 1)   // Reduce Delta
          Delta = .35 * Delta;
        else if (red && (ls.get_slope_ratio() > .7))   // Increase Delta
          Delta = 3 * Delta;

        //  Update ||g||
        ng = max(abs(g)); // norm(g,inf);

        if (trace) { // store trace
          x_values(no_iter) = x;
          F_values(no_iter) = F;
          ng_values(no_iter) = ng;
          Delta_values(no_iter) = Delta;
        }

        no_iter++;
        h = x - xp;
        nh = norm(h);

        //if  (nh == 0)
        //  found = 4;
        //else {
        y = g - gp;
        yh = dot(y, h);
        if (yh > std::sqrt(eps) * nh * norm(y)) {
          //  Update  D
          v = D * y;
          yv = dot(y, v);
          a = (1 + yv / yh) / yh;
          w = (a / 2) * h - v / yh;
          D += outer_product(w, h) + outer_product(h, w); //D = D + w*h' + h*w';
        }  // update D
        //  Check stopping criteria
        double thrx = stop_epsilon_2 * (stop_epsilon_2 + norm(x));
        if (ng <= stop_epsilon_1)
          stop = true; // stop = 1, stop by small gradient
        else if (nh <= thrx)
          stop = true; // stop = 2, stop by small x-step
        else if (no_feval >= max_evaluations)
          stop = true; // stop = 3, number of function evaluations exeeded
        else
          Delta = std::max(Delta, 2 * thrx);
        //} found =4
      }  // Nonzero h
    } // nofail
  }  // iteration

  //  Set return values
  x_end = x;
  finished = true;

  if (trace) { // trim size of trace output
    x_values.set_size(no_iter, true);
    F_values.set_size(no_iter, true);
    ng_values.set_size(no_iter, true);
    Delta_values.set_size(no_iter, true);
  }

  return true;
}

bool Newton_Search::search(vec &xn)
{
  bool state = search();
  xn = get_solution();
  return state;
}

bool Newton_Search::search(const vec &x0, vec &xn)
{
  set_start_point(x0);
  bool state = search();
  xn = get_solution();
  return state;
}

vec Newton_Search::get_solution()
{
  it_assert(finished, "Newton_Search: search is not run yet");
  return x_end;
}

double Newton_Search::get_function_value()
{
  if (finished)
    return F;
  else
    it_warning("Newton_Search::get_function_value, search has not been run");

  return 0.0;
}

double Newton_Search::get_stop_1()
{
  if (finished)
    return ng;
  else
    it_warning("Newton_Search::get_stop_1, search has not been run");

  return 0.0;
}

double Newton_Search::get_stop_2()
{
  if (finished)
    return nh;
  else
    it_warning("Newton_Search::get_stop_2, search has not been run");

  return 0.0;
}

int Newton_Search::get_no_iterations()
{
  if (finished)
    return no_iter;
  else
    it_warning("Newton_Search::get_no_iterations, search has not been run");

  return 0;
}

int Newton_Search::get_no_function_evaluations()
{
  if (finished)
    return no_feval;
  else
    it_warning("Newton_Search::get_no_function_evaluations, search has not been run");

  return 0;
}


void Newton_Search::get_trace(Array<vec> & xvalues, vec &Fvalues, vec &ngvalues, vec &dvalues)
{
  if (finished) {
    if (trace) { // trim size of trace output
      xvalues = x_values;
      Fvalues = F_values;
      ngvalues = ng_values;
      dvalues = Delta_values;
    }
    else
      it_warning("Newton_Search::get_trace, trace is not enabled");
  }
  else
    it_warning("Newton_Search::get_trace, search has not been run");
}

//================================== Line_Search =============================================

Line_Search::Line_Search()
{
  method = Soft;

  if (method == Soft) {
    stop_rho = 1e-3;
    stop_beta = 0.99;
  }

  max_iterations = 10;
  max_stepsize = 10;

  f = NULL;
  df_dx = NULL;
  no_feval = 0;
  init = false;
  finished = false;
  trace = false;
}

void Line_Search::set_function(double(*function)(const vec&))
{
  // Add checks to see that function is OK???
  f = function;
}

void Line_Search::set_gradient(vec(*gradient)(const vec&))
{
  // Add checks to see that function is OK???
  df_dx = gradient;
}


void Line_Search::set_stop_values(double rho, double beta)
{
  // test input values???
  stop_rho = rho;
  stop_beta = beta;
}


void Line_Search::set_start_point(const vec &x, double F, const vec &g, const vec &h)
{
  // check values ???
  x_start = x;
  F_start = F;
  g_start = g;
  h_start = h;
  n = x.size();

  finished = false;
  init = true;
}

void Line_Search::get_solution(vec &xn, double &Fn, vec &gn)
{
  it_assert(finished, "Line_Search: search is not run yet");

  xn = x_end;
  Fn = F_end;
  gn = g_end;
}

bool Line_Search::search()
{
  it_assert(f != NULL, "Line_Search: Function pointer is not set");
  it_assert(df_dx != NULL, "Line_Search: Gradient function pointer is not set");

  it_assert(init, "Line_search: Starting point is not set");

  // Default return values and simple checks
  x_end = x_start;
  F_end = F_start;
  g_end = g_start;

  // add some checks???
  finished = false;

  vec g;

  // return parameters
  no_feval = 0;
  slope_ratio = 1;



  // Check descent condition
  double dF0 = dot(h_start, g_end);

  if (trace) { // prepare structures to store trace data
    alpha_values.set_size(max_iterations);
    F_values.set_size(max_iterations);
    dF_values.set_size(max_iterations);
    alpha_values(0) = 0;
    F_values(0) = F_end;
    dF_values(0) = dF0;
  }


  if (dF0 >= -10*eps*norm(h_start)*norm(g_end)) {  // not significantly downhill
    if (trace) { // store trace
      alpha_values.set_size(1, true);
      F_values.set_size(1, true);
      dF_values.set_size(1, true);
    }
    return false;
  }

  // Finish initialization
  double F0 = F_start, slope0, slopethr;

  if (method == Soft) {
    slope0 = stop_rho * dF0;
    slopethr = stop_beta * dF0;
  }
  else { // exact line search
    slope0 = 0;
    slopethr = stop_rho * std::abs(dF0);
  }

  // Get an initial interval for am
  double a = 0, Fa = F_end, dFa = dF0;
  bool stop = false;
  double b = std::min(1.0, max_stepsize), Fb = 0, dFb = 0;


  while (!stop) {
    Fb = f(x_start + b * h_start);
    g = df_dx(x_start + b * h_start);
    // check if these values are OK if not return false???
    no_feval++;

    dFb = dot(g, h_start);
    if (trace) { // store trace
      alpha_values(no_feval) = b;
      F_values(no_feval) = Fb;
      dF_values(no_feval) = dFb;
    }

    if (Fb < F0 + slope0*b) {  // new lower bound
      alpha = b;
      slope_ratio = dFb / dF0; // info(2);

      if (method == Soft) {
        a = b;
        Fa = Fb;
        dFa = dFb;
      }

      x_end = x_start + b * h_start;
      F_end = Fb;
      g_end = g;

      if ((dFb < std::min(slopethr, 0.0)) && (no_feval < max_iterations) && (b < max_stepsize)) {
        // Augment right hand end
        if (method == Exact) {
          a = b;
          Fa = Fb;
          dFa = dFb;
        }
        if (2.5*b >= max_stepsize)
          b = max_stepsize;
        else
          b = 2 * b;
      }
      else
        stop = true;
    }
    else
      stop = true;
  } // phase 1: expand interval



  if (stop)  // OK so far.  Check stopping criteria
    stop = (no_feval >= max_iterations)
           || (b >= max_stepsize && dFb < slopethr)
           || (a > 0 && dFb >= slopethr);
  // Commented by ediap 2006-07-17: redundant check
  //  || ( (method == Soft) && (a > 0 & dFb >= slopethr) );  // OK


  if (stop && trace) {
    alpha_values.set_size(no_feval, true);
    F_values.set_size(no_feval, true);
    dF_values.set_size(no_feval, true);
  }

  // Refine interval
  while (!stop) {

    double c, Fc, dFc;

    //c = interpolate(xfd,n);
    double C = Fb - Fa - (b - a) * dFa;
    if (C >= 5*n*eps*b) {
      double A = a - 0.5 * dFa * (sqr(b - a) / C);
      c = std::min(std::max(a + 0.1 * (b - a), A), b - 0.1 * (b - a));  // % Ensure significant resuction
    }
    else
      c = (a + b) / 2;

    Fc = f(x_start + c * h_start);
    g = df_dx(x_start + c * h_start);
    dFc = dot(g, h_start);
    // check these values???
    no_feval++;

    if (trace) { // store trace
      alpha_values(no_feval) = c;
      F_values(no_feval) = Fc;
      dF_values(no_feval) = dFc;
    }

    if (method == Soft) {
      // soft line method
      if (Fc < F0 + slope0*c) {  // new lower bound
        alpha = c;
        slope_ratio = dFc / dF0;

        x_end = x_start + c * h_start;
        F_end = Fc;
        g_end = g;
        a = c;
        Fa = Fc;
        dFa = dFc; // xfd(:,1) = xfd(:,3);
        stop = (dFc > slopethr);
      }
      else { // new upper bound
        b = c;
        Fb = Fc;
        dFb = dFc; // xfd(:,2) = xfd(:,3);
      }

    }
    else { // Exact line search
      if (Fc < F_end) {  // better approximant
        alpha = c;
        slope_ratio = dFc / dF0;
        x_end = x_start + c * h_start;
        F_end = Fc;
        g_end = g;
      }
      if (dFc < 0) {  // new lower bound
        a = c;
        Fa = Fc;
        dFa = dFc; // xfd(:,1) = xfd(:,3);
      }
      else { //new upper bound
        b = c;
        Fb = Fc;
        dFb = dFc; // xfd(:,2) = xfd(:,3);
      }
      stop = (std::abs(dFc) <= slopethr) | ((b - a) < stop_beta * b);
    }

    stop = (stop | (no_feval >= max_iterations));
  } // refine

  finished = true;

  if (trace) { // store trace
    alpha_values.set_size(no_feval + 1, true);
    F_values.set_size(no_feval + 1, true);
    dF_values.set_size(no_feval + 1, true);
  }

  return true;
}

bool Line_Search::search(vec &xn, double &Fn, vec &gn)
{
  bool state = search();
  get_solution(xn, Fn, gn);
  return state;
}

bool Line_Search::search(const vec &x, double F, const vec &g, const vec &h,
                         vec &xn, double &Fn, vec &gn)
{
  set_start_point(x, F, g, h);
  bool state = search();
  get_solution(xn, Fn, gn);
  return state;
}


double Line_Search::get_alpha()
{
  if (finished)
    return alpha;
  else
    it_warning("Line_Search::get_alpha, search has not been run");

  return 0.0;
}

double Line_Search::get_slope_ratio()
{
  if (finished)
    return slope_ratio;
  else
    it_warning("Line_Search::get_slope_raio, search has not been run");

  return 0.0;
}

int Line_Search::get_no_function_evaluations()
{
  if (finished)
    return no_feval;
  else
    it_warning("Line_Search::get_no_function_evaluations, search has not been run");

  return 0;
}


void Line_Search::set_max_iterations(int value)
{
  it_assert(value > 0, "Line_Search, max iterations must be > 0");
  max_iterations = value;
}

void Line_Search::set_max_stepsize(double value)
{
  it_assert(value > 0, "Line_Search, max stepsize must be > 0");
  max_stepsize = value;
}

void Line_Search::set_method(const Line_Search_Method &search_method)
{
  method = search_method;

  if (method == Soft) {
    stop_rho = 1e-3;
    stop_beta = 0.99;
  }
  else { // exact line search
    method = Exact;
    stop_rho = 1e-3;
    stop_beta = 1e-3;
  }
}


void Line_Search::get_trace(vec &alphavalues, vec &Fvalues, vec &dFvalues)
{
  if (finished) {
    if (trace) { // trim size of trace output
      alphavalues = alpha_values;
      Fvalues = F_values;
      dFvalues = dF_values;
    }
    else
      it_warning("Line_Search::get_trace, trace is not enabled");
  }
  else
    it_warning("Line_Search::get_trace, search has not been run");
}

// =========================== functions ==============================================

vec fminunc(double(*function)(const vec&), vec(*gradient)(const vec&), const vec &x0)
{
  Newton_Search newton;
  newton.set_functions(function, gradient);

  vec xn;
  newton.search(x0, xn);

  return xn;
}



} // namespace itpp
