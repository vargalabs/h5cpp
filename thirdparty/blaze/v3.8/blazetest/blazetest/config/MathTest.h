//=================================================================================================
/*!
//  \file blazetest/config/MathTest.h
//  \brief General configuration file for the math tests of the blaze test suite
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================


//*************************************************************************************************
/*!\brief First data type for the math tests.
//
// The math tests are always conducted with two data types to test mixed type arithmetic
// expressions. This type definition defines the first of these two data types.
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TYPE_A=int ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TYPE_A int
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TYPE_A
#define BLAZETEST_MATHTEST_TYPE_A int
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Second data type for the math tests.
//
// The math tests are always conducted with two data types to test mixed type arithmetic
// expressions. This type definition defines the second of these two data types.
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TYPE_B=double ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TYPE_B double
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TYPE_B
#define BLAZETEST_MATHTEST_TYPE_B double
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the addition tests.
//
// This compilation switch triggers the addition test scenarios. In case the addition tests are
// activated, addition tests in combination with all activated operations are included in the
// tests. In case the addition tests are disabled, all kinds of addition tests are excluded.
//
// The following settings are possible:
//
//   - 0: The addition tests are not included in the compilation process and not executed
//   - 1: The addition tests are included in the compilation process, but not executed
//   - 2: The addition tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_ADDITION=2 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_ADDITION 2
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_ADDITION
#define BLAZETEST_MATHTEST_TEST_ADDITION 2
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the subtraction tests.
//
// This compilation switch triggers the subtraction test scenarios. In case the subtraction tests
// are activated, subtraction tests in combination with all activated operations are included in
// the tests. In case the subtraction tests are disabled, all kinds of subtraction tests are
// excluded.
//
// The following settings are possible:
//
//   - 0: The subtraction tests are not included in the compilation process and not executed
//   - 1: The subtraction tests are included in the compilation process, but not executed
//   - 2: The subtraction tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_SUBTRACTION=2 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_SUBTRACTION 2
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_SUBTRACTION
#define BLAZETEST_MATHTEST_TEST_SUBTRACTION 2
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the multiplication tests.
//
// This compilation switch triggers the multiplication test scenarios. In case the multiplication
// tests are activated, multiplication tests in combination with all activated operations are
// included in the tests. In case the multiplication tests are disabled, all kinds of subtraction
// tests are excluded.
//
// The following settings are possible:
//
//   - 0: The multiplication tests are not included in the compilation process and not executed
//   - 1: The multiplication tests are included in the compilation process, but not executed
//   - 2: The multiplication tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_MULTIPLICATION=2 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_MULTIPLICATION 2
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_MULTIPLICATION
#define BLAZETEST_MATHTEST_TEST_MULTIPLICATION 2
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the division tests.
//
// This compilation switch triggers the division test scenarios. In case the division tests are
// activated, division tests in combination with all activated operations are included in the
// tests. In case the division tests are disabled, all kinds of division tests are excluded.
//
// The following settings are possible:
//
//   - 0: The division tests are not included in the compilation process and not executed
//   - 1: The division tests are included in the compilation process, but not executed
//   - 2: The division tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_DIVISION=2 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_DIVISION 2
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_DIVISION
#define BLAZETEST_MATHTEST_TEST_DIVISION 2
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the minimum tests.
//
// This compilation switch triggers the minimum test scenarios. In case the minimum tests are
// activated, minimum tests in combination with all activated operations are included in the
// tests. In case the minimum tests are disabled, all kinds of minimum tests are excluded.
//
// The following settings are possible:
//
//   - 0: The minimum tests are not included in the compilation process and not executed
//   - 1: The minimum tests are included in the compilation process, but not executed
//   - 2: The minimum tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_MINIMUM=2 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_MINIMUM 2
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_MINIMUM
#define BLAZETEST_MATHTEST_TEST_MINIMUM 2
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the maximum tests.
//
// This compilation switch triggers the maximum test scenarios. In case the maximum tests are
// activated, maximum tests in combination with all activated operations are included in the
// tests. In case the maximum tests are disabled, all kinds of maximum tests are excluded.
//
// The following settings are possible:
//
//   - 0: The maximum tests are not included in the compilation process and not executed
//   - 1: The maximum tests are included in the compilation process, but not executed
//   - 2: The maximum tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_MAXIMUM=2 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_MAXIMUM 2
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_MAXIMUM
#define BLAZETEST_MATHTEST_TEST_MAXIMUM 2
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the basic tests.
//
// This compilation switch triggers the basic tests for all test scenarios. In case the
// basic tests are activated, each operation is tested separately and without any other
// combined operation. The following example demonstrates this by means of the vector
// addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = a + b;  // Basic vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The basic tests are not included in the compilation process and not executed
//   - 1: The basic tests are included in the compilation process, but not executed
//   - 2: The basic tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_BASIC_OPERATION=2 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_BASIC_OPERATION 2
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_BASIC_OPERATION
#define BLAZETEST_MATHTEST_TEST_BASIC_OPERATION 2
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the negation tests.
//
// This compilation switch triggers the negation tests for all test scenarios. In case the
// negation tests are activated, each operation is tested in combination with a negation.
// The following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = -( a + b );  // Negated vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The negation tests are not included in the compilation process and not executed
//   - 1: The negation tests are included in the compilation process, but not executed
//   - 2: The negation tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_NEGATED_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION
#define BLAZETEST_MATHTEST_TEST_NEGATED_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the scaling tests.
//
// This compilation switch triggers the scaling tests for all test scenarios. In case the
// scaling tests are activated, each operation is tested in combination with a scalar
// multiplication or division. The following example demonstrates this by means of the
// vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = 1.1 * ( a + b );  // Left-multiplied vector addition
   c = ( a + b ) * 1.1;  // Right-multiplied vector addition
   c = ( a + b ) / 1.1;  // Scalar-divided vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The scaling tests are not included in the compilation process and not executed
//   - 1: The scaling tests are included in the compilation process, but not executed
//   - 2: The scaling tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_SCALED_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_SCALED_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_SCALED_OPERATION
#define BLAZETEST_MATHTEST_TEST_SCALED_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the transpose tests.
//
// This compilation switch triggers the transpose tests for all test scenarios. In case the
// transpose tests are activated, each operation is tested in combination with a transpose
// operation. The following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<double,false> a, b;
   blaze::DynamicVector<double,true> c;
   c = trans( a + b );  // Transpose vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The transpose tests are not included in the compilation process and not executed
//   - 1: The transpose tests are included in the compilation process, but not executed
//   - 2: The transpose tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_TRANS_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_TRANS_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_TRANS_OPERATION
#define BLAZETEST_MATHTEST_TEST_TRANS_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the conjugate transpose tests.
//
// This compilation switch triggers the conjugate transpose tests for all test scenarios. In case
// the conjugate transpose tests are activated, each operation is tested in combination with a
// conjugate transpose operation. The following example demonstrates this by means of the vector
// addition:

   \code
   blaze::DynamicVector< complex<double>, false > a, b;
   blaze::DynamicVector< complex<double>, true > c;
   c = ctrans( a + b );  // Conjugate transpose vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The conjugate transpose tests are not included in the compilation process and not executed
//   - 1: The conjugate transpose tests are included in the compilation process, but not executed
//   - 2: The conjugate transpose tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_CTRANS_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION
#define BLAZETEST_MATHTEST_TEST_CTRANS_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the \a abs tests.
//
// This compilation switch triggers the \a abs tests for all test scenarios. In case the \a abs
// tests are activated, each operation is tested in combination with an \a abs operation. The
// following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<int> a, b, c;
   c = abs( a + b );  // Absolute value vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The \a abs tests are not included in the compilation process and not executed
//   - 1: The \a abs tests are included in the compilation process, but not executed
//   - 2: The \a abs tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_ABS_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_ABS_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_ABS_OPERATION
#define BLAZETEST_MATHTEST_TEST_ABS_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the complex conjugate tests.
//
// This compilation switch triggers the complex conjugate tests for all test scenarios. In case
// the complex conjugate tests are activated, each operation is tested in combination with a
// complex conjugate operation. The following example demonstrates this by means of the vector
// addition:

   \code
   blaze::DynamicVector< complex<double> > a, b, c;
   c = conj( a + b );  // Complex conjugate vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The complex conjugate tests are not included in the compilation process and not executed
//   - 1: The complex conjugate tests are included in the compilation process, but not executed
//   - 2: The complex conjugate tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_CONJ_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_CONJ_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_CONJ_OPERATION
#define BLAZETEST_MATHTEST_TEST_CONJ_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the \a real tests.
//
// This compilation switch triggers the \a real tests for all test scenarios. In case the \a real
// tests are activated, each operation is tested in combination with an \a real operation. The
// following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<int> a, b, c;
   c = real( a + b );  // Real part vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The \a real tests are not included in the compilation process and not executed
//   - 1: The \a real tests are included in the compilation process, but not executed
//   - 2: The \a real tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_REAL_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_REAL_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_REAL_OPERATION
#define BLAZETEST_MATHTEST_TEST_REAL_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the \a imag tests.
//
// This compilation switch triggers the \a imag tests for all test scenarios. In case the \a imag
// tests are activated, each operation is tested in combination with an \a imag operation. The
// following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<int> a, b, c;
   c = imag( a + b );  // Imaginary part vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The \a imag tests are not included in the compilation process and not executed
//   - 1: The \a imag tests are included in the compilation process, but not executed
//   - 2: The \a imag tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_IMAG_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_IMAG_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_IMAG_OPERATION
#define BLAZETEST_MATHTEST_TEST_IMAG_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the inversion tests.
//
// This compilation switch triggers the inversion tests for all test scenarios. In case the
// inversion tests are activated, each operation is tested in combination with an inversion
// operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double,false> A, B, C;
   C = inv( A + B );  // Inverted matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The inversion tests are not included in the compilation process and not executed
//   - 1: The inversion tests are included in the compilation process, but not executed
//   - 2: The inversion tests are included in the compilation process and executed
//
// \note In case the inversion tests are activated both \a TypeA and \a TypeB must be set to
// BLAS compatible data types (i.e. float, double, complex<float>, or complex<double>).
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_INV_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_INV_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_INV_OPERATION
#define BLAZETEST_MATHTEST_TEST_INV_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the \a eval tests.
//
// This compilation switch triggers the \a eval tests for all test scenarios. In case the \a eval
// tests are activated, each operation is tested in combination with an \a eval operation. The
// following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = eval( a + b );  // Explicit evaluation of the vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The \a eval tests are not included in the compilation process and not executed
//   - 1: The \a eval tests are included in the compilation process, but not executed
//   - 2: The \a eval tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_EVAL_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_EVAL_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_EVAL_OPERATION
#define BLAZETEST_MATHTEST_TEST_EVAL_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the serial tests.
//
// This compilation switch triggers the serial tests for all test scenarios. In case the serial
// tests are activated, each operation is tested in combination with an serial operation. The
// following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = serial( a + b );  // Explicit serialization of the vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The serial tests are not included in the compilation process and not executed
//   - 1: The serial tests are included in the compilation process, but not executed
//   - 2: The serial tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_SERIAL_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_SERIAL_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_SERIAL_OPERATION
#define BLAZETEST_MATHTEST_TEST_SERIAL_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the \a noalias tests.
//
// This compilation switch triggers the \a noalias tests for all test scenarios. In case the
// \a noalias tests are activated, each operation is tested in combination with an \a noalias
// operation. The following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = noalias( a + b );  // Non-aliased evaluation of the vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The \a noalias tests are not included in the compilation process and not executed
//   - 1: The \a noalias tests are included in the compilation process, but not executed
//   - 2: The \a noalias tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_NOALIAS_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_NOALIAS_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_NOALIAS_OPERATION
#define BLAZETEST_MATHTEST_TEST_NOALIAS_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the \a nosimd tests.
//
// This compilation switch triggers the \a nosimd tests for all test scenarios. In case the
// \a nosimd tests are activated, each operation is tested in combination with an \a nosimd
// operation. The following example demonstrates this by means of the vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   c = nosimd( a + b );  // Non-SIMD evaluation of the vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The \a nosimd tests are not included in the compilation process and not executed
//   - 1: The \a nosimd tests are included in the compilation process, but not executed
//   - 2: The \a nosimd tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_NOSIMD_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_NOSIMD_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_NOSIMD_OPERATION
#define BLAZETEST_MATHTEST_TEST_NOSIMD_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the declsym tests.
//
// This compilation switch triggers the declsym tests for all test scenarios. In case the declsym
// tests are activated, each operation is tested in combination with a declsym operation. The
// following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   C = declsym( A + B );  // Explicitly declaring the matrix addition as symmetric
   \endcode

// The following settings are possible:
//
//   - 0: The declsym tests are not included in the compilation process and not executed
//   - 1: The declsym tests are included in the compilation process, but not executed
//   - 2: The declsym tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION
#define BLAZETEST_MATHTEST_TEST_DECLSYM_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the declherm tests.
//
// This compilation switch triggers the declherm tests for all test scenarios. In case the
// declherm tests are activated, each operation is tested in combination with a declherm
// operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   C = declherm( A + B );  // Explicitly declaring the matrix addition as Hermitian
   \endcode

// The following settings are possible:
//
//   - 0: The declherm tests are not included in the compilation process and not executed
//   - 1: The declherm tests are included in the compilation process, but not executed
//   - 2: The declherm tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION
#define BLAZETEST_MATHTEST_TEST_DECLHERM_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the decllow tests.
//
// This compilation switch triggers the decllow tests for all test scenarios. In case the
// decllow tests are activated, each operation is tested in combination with a decllow
// operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   C = decllow( A + B );  // Explicitly declaring the matrix addition as lower
   \endcode

// The following settings are possible:
//
//   - 0: The decllow tests are not included in the compilation process and not executed
//   - 1: The decllow tests are included in the compilation process, but not executed
//   - 2: The decllow tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION
#define BLAZETEST_MATHTEST_TEST_DECLLOW_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the declupp tests.
//
// This compilation switch triggers the declupp tests for all test scenarios. In case the
// declupp tests are activated, each operation is tested in combination with a declupp
// operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   C = declupp( A + B );  // Explicitly declaring the matrix addition as upper
   \endcode

// The following settings are possible:
//
//   - 0: The declupp tests are not included in the compilation process and not executed
//   - 1: The declupp tests are included in the compilation process, but not executed
//   - 2: The declupp tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION
#define BLAZETEST_MATHTEST_TEST_DECLUPP_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the decldiag tests.
//
// This compilation switch triggers the decldiag tests for all test scenarios. In case the
// decldiag tests are activated, each operation is tested in combination with a decldiag
// operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   C = decldiag( A + B );  // Explicitly declaring the matrix addition as diagonal
   \endcode

// The following settings are possible:
//
//   - 0: The decldiag tests are not included in the compilation process and not executed
//   - 1: The decldiag tests are included in the compilation process, but not executed
//   - 2: The decldiag tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION
#define BLAZETEST_MATHTEST_TEST_DECLDIAG_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the subvector tests.
//
// This compilation switch triggers the subvector tests for all test scenarios. In case the
// subvector tests are activated, all operations resulting in vectors are tested in combination
// with a \a subvector() operation. The following example gives an impression by means of the
// vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   for( ... )
      subvector( c, ... ) = subvector( a + b, ... );  // Subvector-wise vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The subvector tests are not included in the compilation process and not executed
//   - 1: The subvector tests are included in the compilation process, but not executed
//   - 2: The subvector tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION
#define BLAZETEST_MATHTEST_TEST_SUBVECTOR_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the element selection tests.
//
// This compilation switch triggers the element selection tests for all test scenarios. In case
// the element selection tests are activated, all operations resulting in vectors are tested in
// combination with an \a elements() operation. The following example demonstrates this by means
// of the vector addition:

   \code
   blaze::DynamicVector<double> a, b, c;
   for( ... )
      elements( c, ... ) = elements( a + b, ... );  // Element-wise vector addition
   \endcode

// The following settings are possible:
//
//   - 0: The element selection tests are not included in the compilation process and not executed
//   - 1: The element selection tests are included in the compilation process, but not executed
//   - 2: The element selection tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION
#define BLAZETEST_MATHTEST_TEST_ELEMENTS_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the submatrix tests.
//
// This compilation switch triggers the submatrix tests for all test scenarios. In case the
// submatrix tests are activated, all operations resulting in matrices are tested in combination
// with a \a submatrix() operation. The following example gives an impression by means of the
// matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   for( ... )
      for( ... )
         submatrix( C, ... ) = submatrix( A + B, ... );  // Submatrix-wise matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The submatrix tests are not included in the compilation process and not executed
//   - 1: The submatrix tests are included in the compilation process, but not executed
//   - 2: The submatrix tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION
#define BLAZETEST_MATHTEST_TEST_SUBMATRIX_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the row tests.
//
// This compilation switch triggers the row tests for all test scenarios. In case the row tests
// are activated, all operations resulting in matrices are tested in combination with a \a row()
// operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   for( size_t i=0UL; i<A.rows(); ++i )
      row( C, i ) = row( A + B, i );  // Row-wise matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The row tests are not included in the compilation process and not executed
//   - 1: The row tests are included in the compilation process, but not executed
//   - 2: The row tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_ROW_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_ROW_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_ROW_OPERATION
#define BLAZETEST_MATHTEST_TEST_ROW_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the row selection tests.
//
// This compilation switch triggers the row selection tests for all test scenarios. In case
// the row selection tests are activated, all operations resulting in matrices are tested in
// combination with a \a rows() operation. The following example demonstrates this by means
// of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   for( ... )
      rows( C, ... ) = rows( A + B, ... );  // Rows-wise matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The row selection tests are not included in the compilation process and not executed
//   - 1: The row selection tests are included in the compilation process, but not executed
//   - 2: The row selection tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_ROWS_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_ROWS_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_ROWS_OPERATION
#define BLAZETEST_MATHTEST_TEST_ROWS_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the column tests.
//
// This compilation switch triggers the column tests for all test scenarios. In case the column
// tests are activated, all operations resulting in matrices are tested in combination with a
// \a column() operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   for( size_t i=0UL; i<A.column(); ++i )
      column( C, i ) = column( A + B, i );  // Column-wise matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The column tests are not included in the compilation process and not executed
//   - 1: The column tests are included in the compilation process, but not executed
//   - 2: The column tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_COLUMN_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION
#define BLAZETEST_MATHTEST_TEST_COLUMN_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the column selection tests.
//
// This compilation switch triggers the column selection tests for all test scenarios. In case
// the column selection tests are activated, all operations resulting in matrices are tested in
// combination with a \a columns() operation. The following example demonstrates this by means
// of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   for( ... )
      columns( C, ... ) = columns( A + B, ... );  // Columns-wise matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The column selection tests are not included in the compilation process and not executed
//   - 1: The column selection tests are included in the compilation process, but not executed
//   - 2: The column selection tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION
#define BLAZETEST_MATHTEST_TEST_COLUMNS_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compilation switch for the band tests.
//
// This compilation switch triggers the band tests for all test scenarios. In case the band tests
// are activated, all operations resulting in matrices are tested in combination with a \a band()
// operation. The following example demonstrates this by means of the matrix addition:

   \code
   blaze::DynamicMatrix<double> A, B, C;
   const ptrdiff_t ibegin( 1UL - A.rows() );
   const ptrdiff_t iend  ( A.column() );
   for( ptrdiff_t i=ibegin; i<iend; ++i )
      band( C, i ) = band( A + B, i );  // Band-wise matrix addition
   \endcode

// The following settings are possible:
//
//   - 0: The band tests are not included in the compilation process and not executed
//   - 1: The band tests are included in the compilation process, but not executed
//   - 2: The band tests are included in the compilation process and executed
//
// \note It is possible to specify this compilation switch via command line or by defining this
// symbol manually before including any Blazetest header file:

   \code
   g++ ... -DBLAZETEST_MATHTEST_TEST_BAND_OPERATION=0 ...
   \endcode

   \code
   #define BLAZETEST_MATHTEST_TEST_BAND_OPERATION 0
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_MATHTEST_TEST_BAND_OPERATION
#define BLAZETEST_MATHTEST_TEST_BAND_OPERATION 0
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Number of repetitions for a single test case.
//
// The \a BLAZETEST_REPETITIONS value specifies the number of repetitions for each single test
// case. In each repetition the test case is run with random input values. Therefore a higher
// number of repetitions increases the likelihood of detecting implementationn errors. On the
// downside the execution of all tests takes longer.
//
// \note It is possible to specify the number of repetitions via command line or by defining this
// symbol manually before including any Blaze header file:

   \code
   g++ ... -DBLAZETEST_REPETITIONS=3 ...
   \endcode

   \code
   #define BLAZETEST_REPETITIONS 3UL
   #include <blazetest/config/MathTest.h>
   \endcode
*/
#ifndef BLAZETEST_REPETITIONS
#define BLAZETEST_REPETITIONS 3UL
#endif
//*************************************************************************************************
