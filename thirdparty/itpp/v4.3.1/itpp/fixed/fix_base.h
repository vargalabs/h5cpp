/*!
 * \file
 * \brief Definitions of a base class for fixed-point data types
 * \author Johan Bergman
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

#ifndef FIX_BASE_H
#define FIX_BASE_H

#include <itpp/base/ittypes.h>
#include <itpp/stat/misc_stat.h>
#include <itpp/itexports.h>


namespace itpp
{

/*!
  \addtogroup fixed
  \brief Fixed-Point Data Types
  \author Johan Bergman

  \section fix_contents Contents

  <ul>
  <li> \ref fix_intro
  <li> \ref fix_base
  <ul>
  <li> \ref fix_base_shift
  <li> \ref fix_base_wordlen
  <li> \ref fix_base_emode
  <li> \ref fix_base_omode
  <li> \ref fix_base_qmode
  <li> \ref fix_base_stat
  <li> \ref fix_base_outputmode
  </ul>
  <li> \ref fix_real
  <li> \ref fix_complex
  <li> \ref fix_factory
  <li> \ref fix_ops
  <ul>
  <li> \ref fix_ops_asn
  <li> \ref fix_ops_add
  <li> \ref fix_ops_mult
  <li> \ref fix_ops_shift
  <li> \ref fix_ops_conv
  <li> \ref fix_ops_get
  <li> \ref fix_ops_io
  </ul>
  <li> \ref fix_fcn
  <ul>
  <li> \ref fix_fcn_isfix
  <li> \ref fix_fcn_setfix
  <li> \ref fix_fcn_shiftfix
  <li> \ref fix_fcn_assert
  <li> \ref fix_fcn_unfix
  <li> \ref fix_fcn_to
  <li> \ref fix_fcn_other
  </ul>
  </ul>

  \section fix_intro Introduction

  How to include the support for fixed-point data types in your program:
  \code
  #include <itpp/itbase.h>
  #include <itpp/itfixed.h>

  using namespace itpp;
  \endcode

  Fixed-point data types in IT++:
  <ul>
  <li> itpp::Fix (real-valued, restrictions specified as constructor arguments)
  <li> itpp::Fixed (real-valued, restrictions specified as template arguments)
  <li> itpp::CFix (complex-valued, restrictions specified as constructor arguments)
  <li> itpp::CFixed (complex-valued, restrictions specified as template arguments)
  </ul>

  These classes have a common base class called Fix_Base; see inheritance
  diagram in the itpp::Fix_Base documentation. The following data members are
  inherited from Fix_Base:
  <ul>
  <li> Shift factor
  <li> Word length
  <li> Sign encoding mode
  <li> Overflow mode
  <li> Quantization mode
  <li> Statistics object pointer
  <li> Output mode
  </ul>

  The term "fixed-point restrictions" refers to all these data members except
  for the shift factor, which is considered to be part of the "fixed-point
  number". The shift factor has some resemblance to a binary point. The value of
  the shift factor is set in initializations and assignments and it is modified
  by multiplications, divisions and bit-shifting operations. The shift factor is
  used for checking that both terms have been shifted the same amount in
  additions and subtractions. Also, it is used to "un-shift" the data when a
  fixed-point number is converted to floating point.

  Names of classes and enums have been aligned with the fixed-point data types
  in SystemC to some extent, but the fixed-point data types in IT++ and SystemC
  are quite different. In fact, the fixed-point data types in IT++ probably
  correspond better to the variable-precision integer types in SystemC (with
  one important difference: the fixed-point numbers in IT++ remember the amount
  of bit-shifting that has been applied to them, so that they can be converted
  back to "floating-point magnitude" easily if this is desired). The reason for
  this design choice in IT++ is to make the fixed-point simulations as fast as
  possible. If necessary, the core parts in itbase.h (e.g. Array, Vec and Mat)
  should be able to use some other data type than the ones presented here,
  assuming that a proper itpp::Factory is created for the data type, just like
  itpp::Fix_Factory has been created for these data types.

  Sometimes the documentation for the IT++ fixed-point data types states that
  a function is "useful in templated code". This means that the function
  facilitates writing templated code where the template argument is meant to be
  either a floating-point type (double or complex<double>) or a fixed-point type
  (Fix or CFix), i.e. code which is supposed to support both floating-point and
  fixed-point simulations. For example, the operator >>= is defined for Fix and
  CFix, but not for double and complex<double>, so it might be a better idea to
  use the function rshift_fix which is defined for Fix and CFix as well as
  double and complex<double>.

  For an example program, take a look at tests/fix_test.cpp.

  \section fix_base Fix_Base

  \subsection fix_base_shift Shift factor

  Supported shift factors: -64 ... +63 bit-shifts. 0 is \e default.

  An IT++ fixed-point number consists of a bit representation and a shift
  factor. The shift factor is a member of Fix_Base, while the bit representation
  is a member of the inherited class (Fix or CFix). The shift factor indicates
  the number of bit-shifts that have been performed on the data. A positive
  shift factor means that the data has been left-shifted while a negative shift
  factor means that the data has been right-shifted. For information about how
  the shift factor is affected by different operators, see section \ref fix_ops.

  \subsection fix_base_wordlen Word length

  Supported word lengths: 1 ... 64 bits. 64 is \e default.

  \warning Fix, Fixed, CFix and CFixed always use \e signed 64-bit integers to
  represent the fixed-point data. Therefore it is not recommended to declare
  variables with 64 bits and sign encoding mode US.

  \subsection fix_base_emode Sign encoding mode

  Supported sign encoding modes (itpp::e_mode):
  <ul>
  <li> TC (Two's complement)
  <li> US (Unsigned)
  </ul>

  TC is \e default.

  \warning Fix, Fixed, CFix and CFixed always use \e signed 64-bit integers to
  represent the fixed-point data. Therefore it is not recommended to declare
  variables with 64 bits and sign encoding mode US.

  \subsection fix_base_omode Overflow mode

  Supported overflow modes (itpp::o_mode):
  <ul>
  <li> SAT (Saturation)
  <li> WRAP (Wrap-around)
  </ul>

  WRAP is \e default.

  \note Fix, Fixed, CFix and CFixed apply this restriction during initialization
  and assignments only.

  \subsection fix_base_qmode Quantization mode

  Supported quantization modes (itpp::q_mode), with definitions borrowed from
  SystemC (see SystemC documentation for further details):
  <ul>
  <li> RND (Rounding to plus infinity): Add the most significant deleted bit to
  the remaining bits.
  <li> RND_ZERO (Rounding to zero): If the most significant deleted bit is 1,
  and either the sign bit or at least one other deleted bit is 1, add 1 to
  the remaining bits.
  <li> RND_MIN_INF (Rounding to minus infinity): If the most significant deleted
  bit is 1, and at least one other deleted bit is 1, add 1 to the remaining
  bits.
  <li> RND_INF (Rounding to infinity): If the most significant deleted bit is 1,
  and either the inverted value of the sign bit or at least one other
  deleted bit is 1, add 1 to the remaining bits.
  <li> RND_CONV (Convergent rounding with half-way value rounded to even value):
  If the most significant deleted bit is 1, and either the least
  significant of the remaining bits or at least one other deleted bit is 1,
  add 1 to the remaining bits.
  <li> RND_CONV_ODD (Convergent rounding with half-way value rounded to odd
  value): If the most significant deleted bit is 1, and either the least
  significant of the remaining bits is 0 or at least one other deleted bit
  is 1, add 1 to the remaining bits (not defined in SystemC).
  <li> TRN (Truncation): Just copy the remaining bits.
  <li> TRN_ZERO (Truncation to zero): If the sign bit is 1, and either the most
  significant deleted bit or at least one other deleted bit is 1, add 1 to
  the remaining bits.
  </ul>

  TRN is \e default. RND and TRN are usually the most implementation friendly.
  However, note that it is RND_INF that corresponds to "ordinary rounding" and
  TRN_ZERO that corresponds to "ordinary truncation".

  \note Fix, Fixed, CFix and CFixed apply this restriction during initialization
  and right-shift operations only.

  \subsection fix_base_stat Statistics object pointer

  Supported statistics object pointer values: either a pointer to an itpp::Stat
  object or 0. 0 is \e default.

  The sample method in the statistics object is called during initializations
  and assignments. A single statistics object can collect statistics from more
  than one fixed-point variable.

  \subsection fix_base_outputmode Output mode

  Supported output modes (itpp::output_mode), used by the output stream operator
  \<\<:
  <ul>
  <li> OUTPUT_FIX: Output fixed-point representation only
  <li> OUTPUT_FIX_SHIFT: Output fixed-point representation followed by \<shift\>
  <li> OUTPUT_FLOAT: Output floating-point value
  <li> OUTPUT_FLOAT_SHIFT: Output floating-point value followed by \<\<shift
  </ul>

  OUTPUT_FIX_SHIFT is \e default. Unlike the other modes, output_mode is a
  \e static data member of Fix_Base, i.e. the output_mode is common for all
  fixed-point variables. Use the following commands to change output_mode:
  \code
  Fix_Base::set_output_mode(OUTPUT_FIX);
  Fix_Base::set_output_mode(OUTPUT_FIX_SHIFT);
  Fix_Base::set_output_mode(OUTPUT_FLOAT);
  Fix_Base::set_output_mode(OUTPUT_FLOAT_SHIFT);

  // Alternative using a string parameter
  Fix_Base::set_output_mode("OUTPUT_FIX");
  Fix_Base::set_output_mode("OUTPUT_FIX_SHIFT");
  Fix_Base::set_output_mode("OUTPUT_FLOAT");
  Fix_Base::set_output_mode("OUTPUT_FLOAT_SHIFT");

  // Alternative using an ostream modifier
  cout << OUTPUT_FIX;
  cout << OUTPUT_FIX_SHIFT;
  cout << OUTPUT_FLOAT;
  cout << OUTPUT_FLOAT_SHIFT;
  \endcode

  \section fix_real Fix and Fixed

  Fix and Fixed are real-valued fixed-point data types primarily intended to
  replace \c double when a design is refined from floating- to fixed-point
  implementation. The data is stored in the least significant bits of a 64-bit
  integer variable.

  The following example shows how to declare a two's complement (i.e. a signed)
  20-bit variable with wrap-around as overflow handling with the initial value
  3.14 shifted up 10 bits:
  \code
  Fix a(3.14, 10, 20, TC, WRAP);
  Fixed<20, TC, WRAP> b(3.14, 10);
  \endcode

  Note that Fix takes the initial values as well as the fixed-point restrictions
  as constructor arguments. Fixed also takes the initial values as constructor
  arguments but it takes the fixed-point restrictions as template arguments.
  Choose Fix or Fixed depending on your needs. There are three main reasons why
  you would want to choose Fix instead of Fixed. First, if you want to change
  the fixed-point restrictions for a variable during run time, you have to use
  Fix, since the fixed-point restrictions for Fixed have been "fixed" at compile
  time. Second, if your code is using a lot of templating, you might end up with
  many more template arguments if you use Fixed than you would if you use Fix,
  since each set of fixed-point restrictions that you want to use will
  correspond to another type (based on the class template Fixed) instead of just
  different configurations of a single type (Fix). Third, the vector and matrix
  operations currently work better for Fix than for Fixed.

  \note Fixed is derived from Fix, which means that operators, methods and
  functions for Fix can be used for Fixed as well. However, the functions for
  Vec<Fix> (fixvec) and Mat<Fix> (fixmat) cannot be used for Vec<Fixed> and
  Mat<Fixed>.

  If you choose Fix, you should also read the section \ref fix_factory. If you
  choose Fixed, you may find it convenient to use the following typedefs:
  \code
  typedef Fixed<1, TC, WRAP> fixed1;  // for Fixed with 1 bit
  ...
  typedef Fixed<64, TC, WRAP> fixed64;  // for Fixed with 64 bits

  typedef Fixed<1, US, WRAP> ufixed1;  // for Unsigned Fixed with 1 bit
  ...
  typedef Fixed<64, US, WRAP> ufixed64;  // for Unsigned Fixed with 64 bits

  typedef Fixed<1, TC, SAT> sfixed1;  // for Saturated Fixed with 1 bit
  ...
  typedef Fixed<64, TC, SAT> sfixed64;  // for Saturated Fixed with 64 bits

  typedef Fixed<1, US, SAT> sufixed1;  // for Saturated Unsigned Fixed with 1 bit
  ...
  typedef Fixed<64, US, SAT> sufixed64;  // for Saturated Unsigned Fixed with 64 bits
  \endcode

  \note These typedefs use the default values for quantization mode (TRN) and
  statistics object pointer value (0). Also note that U stands for Unsigned but
  S stands for Saturated, NOT for Signed.

  Declaration corresponding to the above Fixed example but using one of the
  typedefs:
  \code
  fixed20 b(3.14, 10);
  \endcode

  \section fix_complex CFix and CFixed

  CFix and CFixed are complex-valued fixed-point data types primarily intended
  to replace <tt>complex<double></tt> when a design is refined from floating- to
  fixed-point implementation. The data is stored in the least significant bits
  of two 64-bit integer variables: one for the real part and one for the
  imaginary part. The two parts have a common shift factor (the one inherited
  from Fix_Base), so it is not possible to shift only one of them.

  The following example shows two ways to declare a two's complement (i.e. a
  signed) 20-bit variable with wrap-around as overflow handling with the initial
  value 1.11 + 2.22i shifted up 10 bits:
  \code
  CFix a(1.11, 2.22, 10, 20, TC, WRAP);
  CFixed<20, TC, WRAP> b(1.11, 2.22, 10);

  CFix c(complex<double>(1.11, 2.22), 0.0, 10, 20, TC, WRAP);
  CFixed<20, TC, WRAP> d(complex<double>(1.11, 2.22), 0.0, 10);
  \endcode

  \note The shift factor is passed as the third argument to the CFix/CFixed
  constructors. If the first argument is complex, the second argument is a
  dummy (that was set to 0.0 in the examples above).

  Choose CFix or CFixed depending on your needs; see section \ref fix_real. If
  you choose CFix, you should also read the section \ref fix_factory.
  \note CFixed is derived from CFix, which means that operators, methods and
  functions for CFix can be used for CFixed as well. However, the functions for
  Vec<CFix> (cfixvec) and Mat<CFix> (cfixmat) cannot be used for Vec<CFixed> and
  Mat<CFixed>.

  If you choose CFixed, you may find it convenient to use the following
  typedefs, which are predefined in IT++:
  \code
  typedef CFixed<1, TC, WRAP> cfixed1;  // for CFixed with 1 bit
  ...
  typedef CFixed<64, TC, WRAP> cfixed64;  // for CFixed with 64 bits

  typedef CFixed<1, TC, SAT> scfixed1;  // for Saturated CFixed with 1 bit
  ...
  typedef CFixed<64, TC, SAT> scfixed64;  // for Saturated CFixed with 64 bits
  \endcode

  \note These typedefs use the default values for sign encoding mode (TC),
  quantization mode (TRN) and statistics object pointer value (0). Also note
  that S stands for Saturated, NOT for Signed.

  Declarations corresponding to the above CFixed examples but using one of the
  typedefs:
  \code
  cfixed20 b(1.11, 2.22, 10);

  cfixed20 d(complex<double>(1.11, 2.22), 0.0, 10);
  \endcode

  \section fix_factory Fix_Factory

  IF you choose to use Fix/CFix (instead of Fixed/CFixed) AND you want to
  declare an Array/Vec/Mat with Fix/CFix elements AND you wish to specify
  some non-default fixed-point restrictions for these elements (i.e. something
  else than 64-bit word length, two's complement as sign encoding mode,
  wrap-around as overflow mode, truncation as quantization mode, and no
  statistics object), THEN you will need to use a Fix_Factory when declaring the
  Array/Vec/Mat of Fix/CFix. Here is how it works, somewhat simplified: you give
  the fixed-point restrictions as parameters to the Fix_Factory, then you give
  the Fix_Factory as a parameter to the Array/Vec/Mat, and finally the
  Array/Vec/Mat uses the Fix_Factory to create Fix/CFix elements with those
  fixed-point restrictions. All constructors for Array, Vec and Mat can take a
  Fix_Factory (or any other Factory for that matter) as their last, optional
  argument. It is assumed that all elements in the Array/Vec/Mat should have the
  same fixed-point restrictions (and use a common statistics object, if any).
  Note that a Fix_Factory can create both Fix and CFix objects. For information
  on Factory in general and Fix_Factory in particular, see the Detailed
  Descriptions for itpp::Factory and itpp::Fix_Factory, respectively.

  The following example shows how to declare a vector of length 7 with Fix
  elements that are two's complement 20-bit variables with wrap-around as
  overflow handling:
  \code
  Vec<Fix> a(7, FIX20);
  \endcode

  FIX20 is one of many predefined Fix_Factory; see the Detailed Description for
  itpp::Fix_Factory.

  \note One might wonder why the Array/Vec/Mat classes themselves cannot take
  the fixed-point restrictions as parameters directly and create the Fix/CFix
  elements without help from a Fix_Factory. The main benefit with the chosen
  solution is that the Array/Vec/Mat classes are not "contaminated" with
  knowledge (parameters, methods, etc) that is specific to the Fix/CFix types.
  If the user for some reason prefers to use some other type (i.e. a type not
  known by IT++) as the Array/Vec/Mat element type, this should work fine as
  long as he or she creates a corresponding Factory. And this is exactly the
  way that Fix/CFix and Fix_Factory work.

  Fix/CFix should not need to know about Fix_Factory, but for the sake of
  uniform syntax in declarations, an exception has been made:
  \code
  Fix a(FIX20);
  \endcode
  i.e. a Fix/CFix declaration can take a Fix_Factory as an argument, just like
  the Vec<Fix> declaration above did.

  \note All declarations with a Fix_Factory as a constructor argument also work
  in templated code; see the Detailed Description for itpp::Fix_Factory.

  \section fix_ops Operators and methods

  \subsection fix_ops_asn Initialization and assignment

  Fixed-point variables can be initialized with a fixed- or a floating-point
  value:
  \code
  // Initialize a with the floating-point value double(3.14*pow2(10))
  // and word length 20, two's complement, wrap-around and rounding
  Fix a(3.14, 10, 20, TC, WRAP, RND);

  // Initialize b with the fixed-point value a
  // and word length 7, two's complement, wrap-around and rounding
  Fix b(a, 7, TC, WRAP, RND);
  \endcode

  In this example, \c b was initialized with the same value as \c a but with
  smaller word length resulting in overflow, since round(3.14*pow2(10)) does not
  fit in the 7-bit variable \c b.

  \warning All fixed-point data types have default constructors. For Fix/CFix,
  the default constructor gives full word length (64 bits). If you call
  templated functions with Fix/CFix as the template argument, they might use
  these default constructors for declaration of temporary variables, which may
  then contain unrestricted (64-bit) temporary results. If you want these
  temporary results to be restricted, then you may have to modify the function
  to introduce fixed-point restrictions for the temporary results, and even
  introduce new temporary variables in order to avoid more than one operation
  per expression (demonstrated in an example below). For Fixed/CFixed, on the
  other hand, the default constructor gives certain fixed-point restrictions,
  but note that it may still be necessary to introduce other fixed-point
  restrictions as well as new temporary variables in the function.

  The assignment operators =, +=, -=, *=, /=, <<= and >>= are supported. For =,
  +=, -=, *= and /=, the right-hand operand can be another fixed-point variable
  or an integer value. If it is an integer value, it is interpreted as a
  fixed-point value with shift factor 0. The = operator simply copies the
  shift factor from the right side to the left side. For information on how the
  other operators treat the shift factor, see sections \ref fix_ops_add,
  \ref fix_ops_mult and \ref fix_ops_shift.

  If assignment to a scaled \c double is desired (when initialization has
  already taken place), the itpp::Fix::set method can be used.
  \code
  // Initialize c with the floating-point value double(3.14*pow2(10))
  // The shift factor is set to 10
  Fix c(3.14, 10);

  // Set c equal to 123. The shift factor is set to 0
  // Note that the old shift factor 10 is discarded
  c = 123;

  // Set c equal to the integer portion of double(3.14*pow2(10))
  // The shift factor is set to 10 (again)
  c.set(3.14, 10);

  // Achieve the same result using a temporary variable
  // Note that the assignment operator copies the shift factor
  c = Fix(3.14, 10);
  \endcode

  When the floating-point value is quantized, the quantization mode of the
  fixed-point variable (TRN in the example above, since \c c has this
  quantization mode) will be used, unless some other quantization mode (e.g.
  RND) is specified as a third argument to set:
  \code
  c.set(3.14, 10, RND);
  \endcode

  \note If you write templated code, you are better off if you use the set_fix
  function described in section \ref fix_fcn_setfix instead of the set method.

  There are also methods for setting data representation and shift directly:
  itpp::Fix::set_re, itpp::CFix::set_im and itpp::Fix_Base::set_shift. They are
  mainly intended for internal use.

  \subsection fix_ops_add Addition and subtraction

  Addition and subtraction between two fixed-point variables as well as between
  a fixed-point variable and an integer variable is supported. The unary minus
  operator is also defined. For Fix and CFix, several vector and matrix
  operations are also supported.
  \code
  // Declare a fixed-point vector with 7 elements
  // (using the predefined Fix_Factory FIX20)
  Vec<Fix> d(7, FIX20);

  // Set all 7 elements equal to 77 with shift factor 0
  d = Fix(77);

  // Declare an integer vector with 7 elements
  ivec e = "1 2 3 4 5 6 7";

  // Add fixed-point vector d and integer vector e. Both have shift factor 0
  Vec<Fix> f(d + e, FIX20);
  \endcode

  \note The addition and subtraction operators require that both operands have
  the same shift factor, unless at least one of the operands is zero. If \c d
  had been assigned with a different shift factor than 0 in the above example
  (and ASSERT_LEVEL > 0), the addition <tt>d + e</tt> would have failed,
  resulting in termination with the error message "assert_shifts: Different
  shifts not allowed!".

  As hinted earlier, the fixed-point restrictions are applied during
  initialization, assignment and bit-shifting operations only. This means that
  the result of an addition or subtraction is unrestricted (64 bits).
  \code
  Fix g(0, 0, 8, TC, SAT);
  Fix h(100, 0, 8, TC, SAT);
  Fix i(100, 0, 8, TC, SAT);
  Fix j(-100, 0, 8, TC, SAT);

  // The result of h + i is unrestricted (64 bits) but when it is assigned to g,
  // it is restricted according to the fixed-point restrictions of g (8 bits).
  // We get overflow, since 100+100=200 doesn't fit in an 8-bit signed variable.
  // The saturated result will be 127
  g = h + i;

  // But now we don't get overflow since 100+100-100=100 does fit!
  g = h + i + j;

  // If we do want the temporary result to be restricted, we have to make
  // an explicit temporary variable (with appropriate restrictions) for it
  Fix tmp(0, 0, 8, TC, SAT);
  // The first sum will be saturated to 127
  tmp = h + i;
  // The final sum will be 127-100=27, i.e. we got a different
  // result when we introduced a restricted temporary variable
  g = tmp + j;
  \endcode

  \subsection fix_ops_mult Multiplication and division

  Multiplication and division between two fixed-point variables as well as
  between a fixed-point variable and an integer variable is supported. For Fix
  and CFix, several vector and matrix operations are also supported.

  As stated earlier, the fixed-point restrictions are applied during
  initialization, assignment and bit-shifting operations only. This means that
  the result of a multiplication or division is unrestricted (64 bits) in the
  same way as for an addition or subtraction; see section \ref fix_ops_add.

  The resulting shift factor after a multiplication is the sum of the two shift
  factors, while the resulting shift factor after a division is the difference
  between the numerator shift factor and the denominator shift factor. The
  result of a division is always quantized using truncation, i.e. the
  quantization modes of the involved fixed-point variables do not matter. Note
  that sometimes divisions can be replaced with multiplications and/or
  bit-shifting operations; see section \ref fix_ops_shift.

  \warning Complex multiplications and divisions are supported, but they utilize
  temporary variables with full word length (64 bits).

  \subsection fix_ops_shift Bit-shifting

  The <<= and >>= operators are defined for the fixed-point data types. As an
  alternative, you can use the itpp::Fix::lshift and itpp::Fix::rshift methods.
  The appropriate fixed-point restrictions of the variable are applied, i.e.
  for left-shifting the overflow mode is applied and for right-shifting the
  quantization mode is applied. There is also a version of rshift that takes a
  quantization mode as the last argument (but there is no corresponding version
  of the >>= operator since it cannot take an extra argument).
  \code
  // Declare a fixed-point variable with the default quantization mode (TRN)
  Fix a(3.14, 10);

  // Right shift 5 bits using the quantization mode of a (i.e. TRN)
  a.rshift(5);

  // Right shift 5 bits using the specified quantization mode (i.e. RND)
  a.rshift(5, RND);
  \endcode

  \note If you write templated code, you are better off if you use the
  lshift_fix and rshift_fix functions described in section \ref fix_fcn_shiftfix
  instead of the <<= and >>= operators and the lshift and rshift methods.

  The << and >> operators are not defined for the fixed-point data types since
  it is not clear what quantization mode that should be applied for the >>
  operator.

  \subsection fix_ops_conv Conversion

  Fix and Fixed can be converted to double, while CFix and CFixed can be
  converted to complex<double>. The conversion operators "un-shift" the data by
  multiplying the fixed-point bit representation with pow2(-shift).

  The itpp::Fix::unfix and itpp::CFix::unfix methods can always be used:
  \code
  Fix a(3.14, 5);

  cout << a.unfix() << endl;
  \endcode
  The resulting output is 3.125.

  \note If you write templated code, you are better off if you use the
  functions unfix or to<T> described in sections \ref fix_fcn_unfix and
  \ref fix_fcn_to instead of the unfix method.

  Equivalently, the double(Fix&) and complex<double>(CFix&) operators can be
  used, unless you define NO_IMPLICIT_FIX_CONVERSION before you include IT++
  in your program.
  \code
  Fix a(3.14, 5);

  cout << double(a) << endl;
  \endcode
  The resulting output is 3.125.

  Finally, Fix/Fixed can be converted to CFix/CFixed using the appropriate
  CFix/CFixed constructors.

  \subsection fix_ops_get Get data members

  This example shows how to get the data members of fixed-point variables:
  \code
  Fix a;

  int64_t the_bit_representation = a.get_re();
  int the_shift_factor = a.get_shift();
  int the_word_length = a.get_wordlen();
  e_mode the_sign_encoding_mode = a.get_e_mode();
  o_mode the_overflow_mode = a.get_o_mode();
  q_mode the_quantization_mode = a.get_q_mode();

  int64_t max_bit_representation = a.get_max();
  int64_t min_bit_representation = a.get_min();
  \endcode

  \note For CFix and CFixed, you get the bit representation for the imaginary
  part with the method get_im().

  \subsection fix_ops_io Input and output

  The print() method outputs the entire state; see section \ref fix_ops_get.
  \code
  Fix a;
  a.print();
  \endcode

  This code example shows how to input and output fixed-point numbers:
  \code
  CFix a(FIX8);
  a.set(0.0, 0.0, 4);

  cout << "Old a: " << a << endl;
  cout << "New a? ";
  cin >> a;
  cout << "New a: " << a << endl;
  \endcode

  Complex numbers can be input on both the C++ form and the IT++ form:
  \code
  Old a: 0+0i<4>
  New a? 1+2i
  New a: 1+2i<4>
  \endcode
  \code
  Old a: 0+0i<4>
  New a? (1,2)
  New a: 1+2i<4>
  \endcode

  \note The output_mode is OUTPUT_FIX_SHIFT in these examples; see section
  \ref fix_base_outputmode.

  In the above examples, only the data representation was changed, while the
  shift (4) was kept. It is however possible to enter another (positive or
  negative) shift factor as well:
  \code
  Old a: 0+0i<4>
  New a? 1+2i<5>
  New a: 1+2i<5>
  \endcode

  It is also possible to enter a floating-point value and a (positive or
  negative) shift, rather than the data representation and a shift, if a
  slightly different format is used:
  \code
  Old a: 0+0i<4>
  New a? 1+2i<<5
  New a: 32+64i<5>
  \endcode

  The resulting data representation is the entered floating-point value 1+2i
  multiplied by 2^5.

  \note In order to enter a negative shift, write 1+2i<<-5. It is not possible
  to write 1+2i>>5.

  Vectors and matrices support the same formats for fixed-point numbers.
  However, all old elements are assumed to have the same shift factor, and the
  shift factor of the first old element becomes the default shift factor for all
  elements of the vector or matrix. The same holds in e.g. an Array of vectors,
  i.e. the shift factor of the first old element in each vector becomes the
  default shift factor for all elements of that vector. The default shift factor
  in an empty vector or matrix is zero.

  \section fix_fcn Functions

  \subsection fix_fcn_isfix Function is_fix

  The function itpp::is_fix returns true only if the argument is of type Fix or
  CFix or an Array/Vec/Mat of Fix or CFix.
  \code
  Array<Array<Vec<Fix> > > aavf(FIX20);
  bool will_be_true = is_fix(aavf);
  \endcode

  \warning Unfortunately, the function is_fix returns false if the argument is
  of type Fixed or CFixed.

  \subsection fix_fcn_setfix Function set_fix

  The function itpp::set_fix sets <tt>y = x * pow2(n)</tt> if the first argument
  is of type Fix or CFix or an Array/Vec/Mat of Fix or CFix. If the first
  argument is of type double or complex<double> or an Array/Vec/Mat of double or
  complex<double>, the function just sets <tt>y = x</tt>.
  \code
  Fix fix_var(FIX20);
  set_fix(fix_var, 3.14, 10);
  // fix_var will equal the integer portion of 3.14 * pow2(10)

  double double_var(FIX20);
  set_fix(double_var, 3.14, 10);
  // double_var will just equal 3.14
  \endcode

  When the floating-point value is quantized, the quantization mode of the first
  argument (TRN in the example above, since fix_var has this quantization mode)
  will be used, unless some other quantization mode (e.g. RND) is specified as a
  fourth argument to set_fix:
  \code
  set_fix(fix_var, 3.14, 10, RND);
  \endcode

  \subsection fix_fcn_shiftfix Functions lshift_fix and rshift_fix

  The functions itpp::lshift_fix and itpp::rshift_fix perform left and right
  bit-shifts, respectively, if the first argument is of type Fix or CFix or an
  Array/Vec/Mat of Fix or CFix. If the first argument is of type double or
  complex<double> or an Array/Vec/Mat of double or complex<double>, the
  functions do not do anything at all.
  \code
  // This will right-shift fix_var 10 bits
  rshift_fix(fix_var, 10);

  // This will not affect double_var
  rshift_fix(double_var, 10);
  \endcode

  When a fixed-point value is right-shifted using rshift_fix, the quantization
  mode of the first argument (TRN in the example above, since fix_var has this
  quantization mode) will be used, unless some other quantization mode (e.g.
  RND) is specified as a third argument to rshift_fix:
  \code
  rshift_fix(fix_var, 10, RND);
  \endcode

  When a fixed-point value is left-shifted using lshift_fix, on the other hand,
  the overflow mode of the first argument is always used.

  \subsection fix_fcn_assert Function assert_fixshift

  The itpp::assert_fixshift function can be used to verify that the shift factor
  has the expected value:
  \code
  Fix a(3.14, 5);

  // We will pass this check since 5 = 5
  assert_fixshift(a, 5);

  // The program will terminate (if ASSERT_LEVEL > 0) since 5 != 6
  assert_fixshift(a, 6);
  \endcode

  If the first argument is of type double or complex<double> instead, no test
  will be performed (since they have no shift factors).

  \subsection fix_fcn_unfix Function unfix

  The itpp::unfix function converts a fixed-point variable to a floating-point
  variable (or an Array/Vec/Mat of fixed-point variables to an Array/Vec/Mat of
  floating-point variables) by multiplying the fixed-point bit representation
  with pow2(-shift), using the unfix method mentioned above.
  \code
  Fix a(3.14, 5);
  double b = unfix(a);

  Array<Mat<CFix> > c(FIX40);
  cin >> c;
  Array<cmat> d = unfix(c);
  \endcode

  If the argument is a floating-point variable (or an Array/Vec/Mat of
  floating-point variables) instead, the function just returns the argument.

  \subsection fix_fcn_to Function to<T>

  The itpp::to<T> function is a very general conversion function. It converts a
  floating- or fixed-point variable to a floating- or fixed-point variable (or
  a floating- or fixed-point Array/Vec/Mat to a floating- or fixed-point
  Array/Vec/Mat) depending on the types of the argument and of the specified
  template parameter.
  \code
  // Convert a Vec<double> to a Vec<Fix> and assign it to f
  Vec<double> e = "1.0 2.0 3.0";
  Vec<Fix> f;
  f = to<Fix>(e);  // convert e "to Fix"

  // Convert an Array<Array<Mat<Fix> > > called g to
  // an Array<Array<Mat<CFix> > > and assign it to h
  Array<Array<Mat<CFix> > > h;
  h = to<CFix>(g);  // convert g "to CFix"
  \endcode

  \note The variants to<double>(x) and to<complex<double> >(x) provide the same
  functionality as unfix(x).

  \warning Be aware that the variants to<Fix>(x) and to<CFix>(x) will return
  a fixed-point variable with the shift factor(s) set to zero if x is a
  floating-point variable. However, if x is a fixed-point variable, the shift
  will be copied to the return variable.

  \subsection fix_fcn_other Other functions

  <ul>
  <li> itpp::abs
  <li> itpp::real
  <li> itpp::imag
  <li> itpp::conj
  </ul>
*/
//!@{

//! Representation for fixed-point data types
typedef int64_t fixrep;
//! Max word length
const int MAX_WORDLEN = 64;

//! Table for fast multiplication or division by 2^n
const uint64_t UINT64_POW2[64] = {
  uint64_t(1),     uint64_t(1) << 1,  uint64_t(1) << 2,  uint64_t(1) << 3,  uint64_t(1) << 4,
  uint64_t(1) << 5,  uint64_t(1) << 6,  uint64_t(1) << 7,  uint64_t(1) << 8,  uint64_t(1) << 9,
  uint64_t(1) << 10, uint64_t(1) << 11, uint64_t(1) << 12, uint64_t(1) << 13, uint64_t(1) << 14,
  uint64_t(1) << 15, uint64_t(1) << 16, uint64_t(1) << 17, uint64_t(1) << 18, uint64_t(1) << 19,
  uint64_t(1) << 20, uint64_t(1) << 21, uint64_t(1) << 22, uint64_t(1) << 23, uint64_t(1) << 24,
  uint64_t(1) << 25, uint64_t(1) << 26, uint64_t(1) << 27, uint64_t(1) << 28, uint64_t(1) << 29,
  uint64_t(1) << 30, uint64_t(1) << 31, uint64_t(1) << 32, uint64_t(1) << 33, uint64_t(1) << 34,
  uint64_t(1) << 35, uint64_t(1) << 36, uint64_t(1) << 37, uint64_t(1) << 38, uint64_t(1) << 39,
  uint64_t(1) << 40, uint64_t(1) << 41, uint64_t(1) << 42, uint64_t(1) << 43, uint64_t(1) << 44,
  uint64_t(1) << 45, uint64_t(1) << 46, uint64_t(1) << 47, uint64_t(1) << 48, uint64_t(1) << 49,
  uint64_t(1) << 50, uint64_t(1) << 51, uint64_t(1) << 52, uint64_t(1) << 53, uint64_t(1) << 54,
  uint64_t(1) << 55, uint64_t(1) << 56, uint64_t(1) << 57, uint64_t(1) << 58, uint64_t(1) << 59,
  uint64_t(1) << 60, uint64_t(1) << 61, uint64_t(1) << 62, uint64_t(1) << 63
};

//! Table for fast multiplication by 2^(n-64)
const double DOUBLE_POW2[128] = {
  0.5 / UINT64_POW2[63], 1.0 / UINT64_POW2[63], 1.0 / UINT64_POW2[62], 1.0 / UINT64_POW2[61],
  1.0 / UINT64_POW2[60], 1.0 / UINT64_POW2[59], 1.0 / UINT64_POW2[58], 1.0 / UINT64_POW2[57],
  1.0 / UINT64_POW2[56], 1.0 / UINT64_POW2[55], 1.0 / UINT64_POW2[54], 1.0 / UINT64_POW2[53],
  1.0 / UINT64_POW2[52], 1.0 / UINT64_POW2[51], 1.0 / UINT64_POW2[50], 1.0 / UINT64_POW2[49],
  1.0 / UINT64_POW2[48], 1.0 / UINT64_POW2[47], 1.0 / UINT64_POW2[46], 1.0 / UINT64_POW2[45],
  1.0 / UINT64_POW2[44], 1.0 / UINT64_POW2[43], 1.0 / UINT64_POW2[42], 1.0 / UINT64_POW2[41],
  1.0 / UINT64_POW2[40], 1.0 / UINT64_POW2[39], 1.0 / UINT64_POW2[38], 1.0 / UINT64_POW2[37],
  1.0 / UINT64_POW2[36], 1.0 / UINT64_POW2[35], 1.0 / UINT64_POW2[34], 1.0 / UINT64_POW2[33],
  1.0 / UINT64_POW2[32], 1.0 / UINT64_POW2[31], 1.0 / UINT64_POW2[30], 1.0 / UINT64_POW2[29],
  1.0 / UINT64_POW2[28], 1.0 / UINT64_POW2[27], 1.0 / UINT64_POW2[26], 1.0 / UINT64_POW2[25],
  1.0 / UINT64_POW2[24], 1.0 / UINT64_POW2[23], 1.0 / UINT64_POW2[22], 1.0 / UINT64_POW2[21],
  1.0 / UINT64_POW2[20], 1.0 / UINT64_POW2[19], 1.0 / UINT64_POW2[18], 1.0 / UINT64_POW2[17],
  1.0 / UINT64_POW2[16], 1.0 / UINT64_POW2[15], 1.0 / UINT64_POW2[14], 1.0 / UINT64_POW2[13],
  1.0 / UINT64_POW2[12], 1.0 / UINT64_POW2[11], 1.0 / UINT64_POW2[10], 1.0 / UINT64_POW2[9],
  1.0 / UINT64_POW2[8],  1.0 / UINT64_POW2[7],  1.0 / UINT64_POW2[6],  1.0 / UINT64_POW2[5],
  1.0 / UINT64_POW2[4],  1.0 / UINT64_POW2[3],  1.0 / UINT64_POW2[2],  1.0 / UINT64_POW2[1],
  1.0,                 1.0*UINT64_POW2[1],  1.0*UINT64_POW2[2],  1.0*UINT64_POW2[3],
  1.0*UINT64_POW2[4],  1.0*UINT64_POW2[5],  1.0*UINT64_POW2[6],  1.0*UINT64_POW2[7],
  1.0*UINT64_POW2[8],  1.0*UINT64_POW2[9],  1.0*UINT64_POW2[10], 1.0*UINT64_POW2[11],
  1.0*UINT64_POW2[12], 1.0*UINT64_POW2[13], 1.0*UINT64_POW2[14], 1.0*UINT64_POW2[15],
  1.0*UINT64_POW2[16], 1.0*UINT64_POW2[17], 1.0*UINT64_POW2[18], 1.0*UINT64_POW2[19],
  1.0*UINT64_POW2[20], 1.0*UINT64_POW2[21], 1.0*UINT64_POW2[22], 1.0*UINT64_POW2[23],
  1.0*UINT64_POW2[24], 1.0*UINT64_POW2[25], 1.0*UINT64_POW2[26], 1.0*UINT64_POW2[27],
  1.0*UINT64_POW2[28], 1.0*UINT64_POW2[29], 1.0*UINT64_POW2[30], 1.0*UINT64_POW2[31],
  1.0*UINT64_POW2[32], 1.0*UINT64_POW2[33], 1.0*UINT64_POW2[34], 1.0*UINT64_POW2[35],
  1.0*UINT64_POW2[36], 1.0*UINT64_POW2[37], 1.0*UINT64_POW2[38], 1.0*UINT64_POW2[39],
  1.0*UINT64_POW2[40], 1.0*UINT64_POW2[41], 1.0*UINT64_POW2[42], 1.0*UINT64_POW2[43],
  1.0*UINT64_POW2[44], 1.0*UINT64_POW2[45], 1.0*UINT64_POW2[46], 1.0*UINT64_POW2[47],
  1.0*UINT64_POW2[48], 1.0*UINT64_POW2[49], 1.0*UINT64_POW2[50], 1.0*UINT64_POW2[51],
  1.0*UINT64_POW2[52], 1.0*UINT64_POW2[53], 1.0*UINT64_POW2[54], 1.0*UINT64_POW2[55],
  1.0*UINT64_POW2[56], 1.0*UINT64_POW2[57], 1.0*UINT64_POW2[58], 1.0*UINT64_POW2[59],
  1.0*UINT64_POW2[60], 1.0*UINT64_POW2[61], 1.0*UINT64_POW2[62], 1.0*UINT64_POW2[63]
};

//! Sign encoding modes (aligned with SystemC)
enum e_mode {
  TC,                 //!< Two's complement
  US                  //!< Unsigned
};

//! Overflow modes (aligned with SystemC)
enum o_mode {
  SAT,                //!< Saturation
  SAT_ZERO,           //!< Saturation to zero (Not implemented)
  SAT_SYM,            //!< Symmetrical saturation (Not implemented)
  WRAP,               //!< Wrap-around
  WRAP_SM             //!< Sign magnitued wrap-around (Not implemented)
};

//! Quantization modes (aligned with SystemC)
enum q_mode {
  RND,                //!< Rounding to plus infinity
  RND_ZERO,           //!< Rounding to zero
  RND_MIN_INF,        //!< Rounding to minus infinity
  RND_INF,            //!< Rounding to infinity
  RND_CONV,           //!< Convergent rounding with half-way value rounded to even value
  RND_CONV_ODD,       //!< Convergent rounding with half-way value rounded to odd value (not defined in SystemC)
  TRN,                //!< Truncation
  TRN_ZERO            //!< Truncation to zero
};

//! Output modes
enum output_mode {
  OUTPUT_FIX,         //!< Output fixed-point representation only
  OUTPUT_FIX_SHIFT,   //!< Output fixed-point representation followed by \<shift\> (default)
  OUTPUT_FLOAT,       //!< Output floating-point value
  OUTPUT_FLOAT_SHIFT  //!< Output floating-point value followed by \<\<shift
};

/*!
  \brief Base class for fixed-point data types

  See the Detailed Description in the \ref fixed module.
*/
class ITPP_EXPORT Fix_Base
{
public:
  //! Default constructor
  explicit Fix_Base(int s = 0, int w = MAX_WORDLEN, e_mode e = TC, o_mode o = WRAP, q_mode q = TRN, Stat *ptr = 0)
      : shift(s), wordlen(w), emode(e), omode(o), qmode(q), stat_ptr(ptr) {init();}
  //! Copy constructor
  Fix_Base(const Fix_Base &x)
      : shift(x.shift), wordlen(MAX_WORDLEN), emode(TC), omode(WRAP), qmode(TRN), stat_ptr(0) {init();}
  //! Destructor
  virtual ~Fix_Base() {}

  //! Set shift (without shifting)
  void set_shift(int s) {shift = s;}
  //! Set output mode to OUTPUT_FIX, OUTPUT_FIX_SHIFT, OUTPUT_FLOAT or OUTPUT_FLOAT_SHIFT. Static member function.
  static void set_output_mode(output_mode o) {outputmode = o;}
  //! Set output mode to "OUTPUT_FIX", "OUTPUT_FIX_SHIFT", "OUTPUT_FLOAT" or "OUTPUT_FLOAT_SHIFT". Static member function.
  static void set_output_mode(std::string o);

  //! Get shift
  int get_shift() const {return shift;}
  //! Get word length
  int get_wordlen() const {return wordlen;}
  //! Get sign encoding mode
  e_mode get_e_mode() const {return emode;}
  //! Get overflow mode
  o_mode get_o_mode() const {return omode;}
  //! Get quantization mode
  q_mode get_q_mode() const {return qmode;}
  //! Get output mode
  output_mode get_output_mode() const {return outputmode;}
  //! Get maximum value of data representation
  fixrep get_max() const {return max;}
  //! Get minimum value of data representation
  fixrep get_min() const {return min;}
  //! Print restrictions
  virtual void print() const;

protected:
  //! Accumulated bitshift (positive means left-shifted, negative means right-shifted)
  int shift;
  //! Word length
  int wordlen;
  //! Sign encoding mode
  e_mode emode;
  //! Overflow mode
  o_mode omode;
  //! Quantization mode
  q_mode qmode;
  //! Pointer to statistics object
  Stat *stat_ptr;
  //! Minimum allowed value (help variable to speed up calculations)
  fixrep min;
  //! Maximum allowed value (help variable to speed up calculations)
  fixrep max;
  //! Number of unused (MSB) bits (help variable to speed up calculations)
  int n_unused_bits;

  //! Calculate help variables min, max and n_unused_bits
  void init();
  //! Handle overflows using overflow mode \c omode and make call to statistics object (if any)
  fixrep apply_o_mode(fixrep x) const;
  //! Convert from \c double to \c fixrep using \c shift and quantization mode \c qmode, then call <tt>limit()</tt>
  fixrep scale_and_apply_modes(double x) const {return scale_and_apply_modes(x, qmode);}
  //! Convert from \c double to \c fixrep using \c shift and quantization mode \c q, then call <tt>limit()</tt>
  fixrep scale_and_apply_modes(double x, q_mode q) const;
  //! Right shift \c n bits using quantization mode \c qmode and make call to statistics object (if any)
  fixrep rshift_and_apply_q_mode(fixrep x, int n) const {return rshift_and_apply_q_mode(x, n, qmode);}
  //! Right shift \c n bits using quantization mode \c q and make call to statistics object (if any)
  fixrep rshift_and_apply_q_mode(fixrep x, int n, q_mode q) const;

private:
  //! Output mode
  static output_mode outputmode;
};

//! Set output mode
inline std::ostream &operator<<(std::ostream &os, const output_mode &o)
{
  Fix_Base::set_output_mode(o);
  return os;
}

//!@}

} // namespace itpp

#endif // #ifndef FIX_BASE_H
