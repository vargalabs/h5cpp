/*!
 * \file
 * \brief Definition of an argument parser class
 * \author Thomas Eriksson, Pal Frenger and Johan Bergman
 *
 * Thanks to Svante Signell for valuable input
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

#ifndef PARSER_H
#define PARSER_H

// #define MAX_STR_LEN 4096

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/array.h>
#include <iostream>
#include <itpp/itexports.h>

namespace itpp
{
//! \cond

#if (defined(_MSC_VER) && defined(ITPP_SHARED_LIB))
//MSVC explicitely instantiate required template while building the shared library
template  class ITPP_EXPORT Array<std::string>;
#endif

//! \endcond

/*!
  \addtogroup parser

  \brief Argument Parser
  \author Thomas Eriksson and Pal Frenger (Thanks to Svante Signell for
          valuable input), modifications by Johan Bergman

  This class parses strings to variables. The syntax is compatible with Matlab
  and Octave. It can be used in several different ways. The following test
  program and test data file gives several examples:

  <ul>
  <li> Use the Parser class on a parameter file </li>
  <li> Use the Parser class on the command line input </li>
  <li> Use the Parser class on a char pointer (usually the command line input) </li>
  <li> Use the Parser class on a parameter file and a char pointer </li>
  <li> Use the Parser class on an Array of strings </li>
  <li> Use the Parser::get() method on a parameter file </li>
  <li> Use the Parser::exist() method to check if variables exist </li>
  </ul>

  The test program looks as follows:
  \include parser_test.cpp

  The data file \c parser_test_data.txt looks as follows:
  \include parser_test_data.txt

  Beside the type-specific get methods, like get_int(), there is a templated
  get method called get(), that can handle all variable types that have an
  istream operator defined, e.g. Arrays. The get() method takes the variable
  value as an input/output parameter. If the variable name cannot be found,
  the old value is kept. Here is an example:

  \code
  // Declare and initialize a variable
  Array<ivec> var;
  set_array(var, "{[1] [2 3] [4 5 6]}");

  // Let the Parser p get the variable named my_var_name
  bool my_var_name_was_found = p.get(var, "my_var_name");
  \endcode

  In the above example, if \c my_var_name was not found, \c var keeps its old
  value {[1] [2 3] [4 5 6]}. In non-silent mode, the get() method echoes the
  variable values to the standard output in Matlab/Octave m-file style.
 */

/*!
  \ingroup parser
  \brief Argument Parser Class
  \author Thomas Eriksson and Pal Frenger (Thanks to Svante Signell for
  valuable input)

  This class parses strings to variables. The syntax is compatible with Matlab
  and Octave. It can be used in several different ways. See the Detailed Description
  in the \ref parser module.
*/

class ITPP_EXPORT Parser
{
public:

  //! Default Constructor
  Parser();

  //! Constructor. Sets input file name.
  Parser(const std::string &filename);

  //! Constructor. Uses argc and argv (command line arguments)
  Parser(int argc, char *argv[]);

  //! Constructor. Sets input file name and uses argc and argv (command line arguments)
  Parser(const std::string &filename, int argc, char *argv[]);

  //! Constructor. Sets and Array of strings
  Parser(const Array<std::string> &setup);

  //! Initialization function. Sets input file name.
  void init(const std::string &filename);

  //! Initialization function. Uses argc and argv (command line arguments)
  void init(int argc, char *argv[]);

  //! Initialization function. Sets input file name and uses argc and argv (command line arguments)
  void init(const std::string &filename, int argc, char *argv[]);

  //! Initialization function. Sets and Array of strings
  void init(const Array<std::string> &setup);

  //! Sets silent mode if true, or verbose mode if false
  void set_silentmode(bool v = true);

  //! Check is \a name exists in the file. Returns \c true if the \a name is found and \c false otherwise.
  bool exist(const std::string &name);

  //! Get variable value if \a name can be found (and return true), otherwise keep old value (and return false)
  template<class T>
  bool get(T &var, const std::string &name, int num = -1);

  //! Interpret variable \a name as a bool
  bool get_bool(const std::string &name,  int num = -1);

  //! Interpret variable \a name as an integer
  int get_int(const std::string &name,  int num = -1);

  //! Interpret variable \a name as a double
  double get_double(const std::string &name, int num = -1);

  //! Interpret variable \a name as a string
  std::string get_string(const std::string &name, int num = -1);

  //! Interpret variable \a name as a vec
  vec get_vec(const std::string &name, int num = -1);

  //! Interpret variable \a name as a ivec
  ivec get_ivec(const std::string &name, int num = -1);

  //! Interpret variable \a name as a svec
  svec get_svec(const std::string &name, int num = -1);

  //! Interpret variable \a name as a bvec
  bvec get_bvec(const std::string &name, int num = -1);

  //! Interpret variable \a name as a mat
  mat get_mat(const std::string &name, int num = -1);

  //! Interpret variable \a name as a imat
  imat get_imat(const std::string &name, int num = -1);

  //! Interpret variable \a name as a smat
  smat get_smat(const std::string &name, int num = -1);

  //! Interpret variable \a name as a bmat
  bmat get_bmat(const std::string &name, int num = -1);

protected:

private:

  //! Find the string \c name
  std::string findname(const std::string &name,
                       bool &error_flag,
                       bool &print_flag,
                       int num = 0,
                       bool keep_brackets = false);

  void pre_parsing(void);

  Array<std::string> SetupStrings;

  bool VERBOSE;
};

// ----------------------- Implementation starts here -----------------------

template<class T>
bool Parser::get(T &var, const std::string &name, int num)
{
  bool error_flag, print_flag;
  std::string str = findname(name, error_flag, print_flag, num, true);
  std::istringstream buffer(str);
  if (error_flag) {
    if (VERBOSE) {
      std::cout << name << " = " << var << ";" << std::endl;
    }
  }
  else {
    buffer >> var;
    if (print_flag) {
      std::cout << name << " = " << var << std::endl;
    }
    else if (VERBOSE) {
      std::cout << name << " = " << var << ";" << std::endl;
    }
  }
  return !error_flag;
}

//! Specialization or \c get() for std::string
template<>
ITPP_EXPORT bool Parser::get(std::string &var, const std::string &name, int num);
//! Specialization of \c get() for int
template<>
ITPP_EXPORT bool Parser::get(int &var, const std::string &name, int num);
//! Specialization of \c get() for bool
template<>
ITPP_EXPORT bool Parser::get(bool &var, const std::string &name, int num);

} // namespace itpp

#endif // #ifndef PARSER_H
