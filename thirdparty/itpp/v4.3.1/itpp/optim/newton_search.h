/*!
 * \file
 * \brief Newton Search optimization algorithms - header file
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

#ifndef NEWTON_SEARCH_H
#define NEWTON_SEARCH_H

#include <itpp/base/vec.h>
#include <itpp/base/array.h>
#include <limits>
#include <itpp/itexports.h>
#include <itpp/base/base_exports.h>

namespace itpp
{

/*!
  \brief Numerical optimization routines
  \addtogroup optimization
*/
//@{


//! Newton Search method
enum Newton_Search_Method {BFGS};

/*!
  \brief Newton Search

  Newton or Quasi-Newton optimization method that try to minimize the objective function \f$f(\mathbf{x})\f$
  given an initial guess \f$\mathbf{x}\f$.

  The search is stopped when either criterion 1:
  \f[
  \left\| \mathbf{f}'(\mathbf{x})\right\|_{\infty} \leq \varepsilon_1
  \f]
  or criterion 2:
  \f[
  \left\| d\mathbf{x}\right\|_{2} \leq \varepsilon_2 (\varepsilon_2 + \| \mathbf{x} \|_{2} )
  \f]
  is fulfilled. Another possibility is that the search is stopped when the number of function evaluations
  exceeds a threshold (100 per default).

  The default update rule for the inverse of the Hessian matrix is the BFGS algorithm with
  \f$\varepsilon_1 = 10^{-4}\f$ an \f$\varepsilon_2 = 10^{-8}\f$.

*/
class ITPP_EXPORT Newton_Search
{
public:
  //! Default constructor
  Newton_Search();
  //! Destructor
  ~Newton_Search() {};

  //! Set function pointer
  void set_function(double(*function)(const vec&));
  //! Set gradient function pointer
  void set_gradient(vec(*gradient)(const vec&));
  //! Set both function and gradient function pointers
  void set_functions(double(*function)(const vec&), vec(*gradient)(const vec&)) { set_function(function); set_gradient(gradient); }

  //! Set start point \c x for search and approx inverse Hessian at \c x
  void set_start_point(const vec &x, const mat &D);

  //! Set start point \c x for search
  void set_start_point(const vec &x);

  //! Get solution, function value and gradient at solution point
  vec get_solution();

  //! Do the line search
  bool search();
  //! Do the line search and return solution
  bool search(vec &xn);
  //! Set starting point, do the Newton search, and return the solution
  bool search(const vec &x0, vec &xn);

  //! Set stop criterion values
  void set_stop_values(double epsilon_1, double epsilon_2);
  //! Return stop value rho
  double get_epsilon_1() { return stop_epsilon_1; }
  //! Return stop value beta
  double get_epsilon_2() { return stop_epsilon_2; }

  //! Set max number of function evaluations
  void set_max_evaluations(int value);
  //! Return max number of function evaluations
  int get_max_evaluations() { return max_evaluations; }

  //! Set max stepsize
  void set_initial_stepsize(double value);
  //! Return max number of iterations
  double get_initial_stepsize() { return initial_stepsize; }

  //! Set Line search method
  void set_method(const Newton_Search_Method &method);

  //! get function value at solution point
  double get_function_value();
  //! get value of stop criterion 1 at solution point
  double get_stop_1();
  //! get value of stop criterion 2 at solution point
  double get_stop_2();
  //! get number of iterations used to reach solution
  int get_no_iterations();
  //! get number of function evaluations used to reach solution
  int get_no_function_evaluations();

  //! enable trace mode
  void enable_trace() { trace = true; }
  //! disable trace
  void disable_trace() { trace = false; }

  /*! get trace outputs
    \c xvalues are the solutions of every iteration
    \c Fvalues are the function values
    \c ngvalues are the norm(gradient,inf) values
    \c dvalues are the delta values
  */
  void get_trace(Array<vec> & xvalues, vec &Fvalues, vec &ngvalues, vec &dvalues);

private:
  int n; // dimension of problem, size(x)
  double(*f)(const vec&);  // function to minimize
  vec(*df_dx)(const vec&);  // df/dx, gradient of f

  // start variables
  vec x_start;
  mat D_start;

  // solution variables
  vec x_end;

  // trace variables
  Array<vec> x_values;
  vec F_values, ng_values, Delta_values;

  Newton_Search_Method method;

  // Parameters
  double initial_stepsize; // opts(1)
  double stop_epsilon_1; // opts(2)
  double stop_epsilon_2; // opt(3)
  int max_evaluations; // opts(4)

  // output parameters
  int no_feval; // number of function evaluations
  int no_iter; // number of iterations
  double F, ng, nh; // function value, stop_1, stop_2 values at solution point

  bool init, finished, trace;
};



//! Line Search method
enum Line_Search_Method {Soft, Exact};

/*!
  \brief Line Search

  The line search try to minimize the objective function \f$f(\mathbf{x})\f$
  along the direction \f$\mathbf{h}\f$ from the current position \f$\mathbf{x}\f$.

  Hence we look at
  \f[
  \varphi(\alpha) = f(\mathbf{x} + \alpha \mathbf{h})
  \f]
  and try to find an \f$\alpha_s\f$ that minimizes \f$f\f$.

  Two variants are used. Either the soft line search (default) or the exact line
  search.

  The soft line search stops when a point in the acceptable region is found, i.e.
  \f[
  \phi(\alpha_s) \leq \varphi(0) + \alpha_s \rho \varphi'(0)
  \f]
  and
  \f[
  \varphi'(\alpha_s) \geq \beta \varphi'(0),\: \rho < \beta
  \f]
  Default vales are \f$\rho = 10^{-3}\f$ and \f$\beta = 0.99\f$.

  The exact line search
  \f[
  \| \varphi(\alpha_s)\|  \leq \rho \| \varphi'(0) \|
  \f]
  and
  \f[
  b-a \leq \beta b,
  \f]
  where \f$\left[a,b\right]\f$ is the current interval for \f$\alpha_s\f$.
  Default vales are \f$\rho = 10^{-3}\f$ and \f$\beta = 10^{-3}\f$.

  The exact line search can at least in theory give the exact resutl, but it may require
  many extra function evaluations compared to soft line search.
*/
class ITPP_EXPORT Line_Search
{
public:
  //! Default constructor
  Line_Search();
  //! Destructor
  ~Line_Search() {};

  //! Set function pointer
  void set_function(double(*function)(const vec&));
  //! Set gradient function pointer
  void set_gradient(vec(*gradient)(const vec&));
  //! Set both function and gradient function pointers
  void set_functions(double(*function)(const vec&), vec(*gradient)(const vec&)) { set_function(function); set_gradient(gradient); }

  //! Set start point for search
  void set_start_point(const vec &x, double F, const vec &g, const vec &h);

  //! Get solution, function value and gradient at solution point
  void get_solution(vec &xn, double &Fn, vec &gn);

  //! Do the line search
  bool search();
  //! Do the line search and return solution
  bool search(vec &xn, double &Fn, vec &gn);
  //! Set starting point, do the line search, and return the solution
  bool search(const vec &x, double F, const vec &g, const vec &h, vec &xn,
              double &Fn, vec &gn);


  //! return alpha at solution point, xn = x + alpha h
  double get_alpha();
  //! return the slope ratio at solution poin, xn
  double get_slope_ratio();
  //! return number of function evaluations used in search
  int get_no_function_evaluations();


  //! Set stop criterion values
  void set_stop_values(double rho, double beta);
  //! Return stop value rho
  double get_rho() { return stop_rho; }
  //! Return stop value beta
  double get_beta() { return stop_beta; }

  //! Set max number of iterations
  void set_max_iterations(int value);
  //! Return max number of iterations
  int get_max_iterations() { return max_iterations; }

  //! Set max stepsize
  void set_max_stepsize(double value);
  //! Return max number of iterations
  double get_max_stepsize() { return max_stepsize; }

  //! Set Line search method
  void set_method(const Line_Search_Method &method);

  //! enable trace mode
  void enable_trace() { trace = true; }
  //! disable trace
  void disable_trace() { trace = false; }

  /*! get trace outputs
    \c alphavalues are the solutions of every iteration
    \c Fvalues are the function values
    \c dFvalues
  */
  void get_trace(vec &alphavalues, vec &Fvalues, vec &dFvalues);

private:
  int n; // dimension of problem, size(x)
  double(*f)(const vec&);  // function to minimize
  vec(*df_dx)(const vec&);  // df/dx, gradient of f

  // start variables
  vec x_start, g_start, h_start;
  double F_start;

  // solution variables
  vec x_end, g_end;
  double F_end;

  // trace variables
  vec alpha_values, F_values, dF_values;

  bool init; // true if functions and starting points are set
  bool finished; // true if functions and starting points are set
  bool trace; // true if trace is enabled

  // Parameters
  Line_Search_Method method;
  double stop_rho; // opts(2)
  double stop_beta; // opts(3)
  int max_iterations; // opts(4)
  double max_stepsize; // opts(5)

  // output parameters
  double alpha; // end value of alpha, info(1)
  double slope_ratio; // slope ratio at xn, info(2)
  int no_feval; // info(3)
};

/*!
  \brief Unconstrained minimization

  Unconstrained minimization using a Newton or Quasi-Newton optimization method
  that try to minimize the objective function \f$f(\mathbf{x})\f$ given an initial guess \f$\mathbf{x}\f$.

  The function and the gradient need to be known and supplied.

  The default algorithm is a Quasi-Newton search using BFGS updates of the inverse Hessian matrix.
*/
ITPP_EXPORT vec fminunc(double(*function)(const vec&), vec(*gradient)(const vec&), const vec &x0);

//@}

} // namespace itpp

#endif // #ifndef NEWTON_SEARCH_H
