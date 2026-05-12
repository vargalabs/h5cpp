//=================================================================================================
/*!
//  \file src/mathtest/matrices/initializermatrix/ClassTest.cpp
//  \brief Source file for the InitializerMatrix class test
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
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <blazetest/mathtest/matrices/initializermatrix/ClassTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


namespace blazetest {

namespace mathtest {

namespace matrices {

namespace initializermatrix {

//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the InitializerMatrix class test.
//
// \exception std::runtime_error Operation error detected.
*/
ClassTest::ClassTest()
{
   testConstructors();
   testFunctionCall();
   testAt();
   testIterator();
   testNonZeros();
   testSwap();
}
//*************************************************************************************************




//=================================================================================================
//
//  TEST FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Test of the InitializerMatrix constructors.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of all constructors of the InitializerMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testConstructors()
{
   using blaze::initializer_list;


   //=====================================================================================
   // Single argument constructor
   //=====================================================================================

   {
      test_ = "InitializerMatrix single argument constructor (0x0)";

      initializer_list< initializer_list<int> > list = {};

      blaze::InitializerMatrix<int> mat( list );

      checkRows    ( mat, 0UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "InitializerMatrix single argument constructor (3x4)";

      initializer_list< initializer_list<int> > list = { { 1, 0, 3, 4 },
                                                         { 0 },
                                                         { 2, 0, 5 } };

      blaze::InitializerMatrix<int> mat( list );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 5UL );
   }


   //=====================================================================================
   // Two argument constructor
   //=====================================================================================

   {
      test_ = "InitializerMatrix two argument constructor (3x0)";

      initializer_list< initializer_list<int> > list = { {},
                                                         {},
                                                         {} };

      blaze::InitializerMatrix<int> mat( list, 0UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 0UL );
      checkNonZeros( mat, 0UL );
   }

   {
      test_ = "InitializerMatrix two argument constructor (3x4)";

      initializer_list< initializer_list<int> > list = { { 1, 0, 3, 4 },
                                                         { 0 },
                                                         { 2, 0, 5 } };

      blaze::InitializerMatrix<int> mat( list, 4UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 4UL );
      checkNonZeros( mat, 5UL );
   }

   {
      test_ = "InitializerMatrix two argument constructor (3x6)";

      initializer_list< initializer_list<int> > list = { { 1, 0, 3, 4 },
                                                         { 0 },
                                                         { 2, 0, 5 } };

      blaze::InitializerMatrix<int> mat( list, 6UL );

      checkRows    ( mat, 3UL );
      checkColumns ( mat, 6UL );
      checkNonZeros( mat, 5UL );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the InitializerMatrix function call operator.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the function call operator
// of the InitializerMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testFunctionCall()
{
   using blaze::initializer_list;


   test_ = "InitializerMatrix::operator()";

   initializer_list< initializer_list<int> > list = { { 1, 0, 3, 4 },
                                                      { 0 },
                                                      { 2, 0, 5 } };

   blaze::InitializerMatrix<int> mat( list, 6UL );

   // Access to the element (0,2)
   if( mat(0,2) != 3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Function call operator failed\n"
          << " Details:\n"
          << "   Result:\n" << mat << "\n"
          << "   Expected result:\n( 1 0 3 4 0 0 )\n( 0 0 0 0 0 0 )\n( 2 0 5 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Access to the element (1,2)
   if( mat(1,2) != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Function call operator failed\n"
          << " Details:\n"
          << "   Result:\n" << mat << "\n"
          << "   Expected result:\n( 1 0 3 4 0 0 )\n( 0 0 0 0 0 0 )\n( 2 0 5 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c at() member function of the InitializerMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of adding and accessing elements via the \c at() member function
// of the InitializerMatrix class template. In case an error is detected, a \a std::runtime_error
// exception is thrown.
*/
void ClassTest::testAt()
{
   using blaze::initializer_list;


   test_ = "InitializerMatrix::operator()";

   initializer_list< initializer_list<int> > list = { { 1, 0, 3, 4 },
                                                      { 0 },
                                                      { 2, 0, 5 } };

   blaze::InitializerMatrix<int> mat( list, 6UL );

   // Access to the element (0,2)
   if( mat.at(0,2) != 3 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Function call operator failed\n"
          << " Details:\n"
          << "   Result:\n" << mat << "\n"
          << "   Expected result:\n( 1 0 3 4 0 0 )\n( 0 0 0 0 0 0 )\n( 2 0 5 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Access to the element (1,2)
   if( mat.at(1,2) != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Function call operator failed\n"
          << " Details:\n"
          << "   Result:\n" << mat << "\n"
          << "   Expected result:\n( 1 0 3 4 0 0 )\n( 0 0 0 0 0 0 )\n( 2 0 5 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }

   // Attempt to access the element (3,0)
   try {
      mat.at(3,0);

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Result:\n" << mat << "\n"
          << "   Expected result:\n( 1 0 3 4 0 0 )\n( 0 0 0 0 0 0 )\n( 2 0 5 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}

   // Attempt to access the element (2,6)
   try {
      mat.at(2,6);

      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Out-of-bound access succeeded\n"
          << " Details:\n"
          << "   Result:\n" << mat << "\n"
          << "   Expected result:\n( 1 0 3 4 0 0 )\n( 0 0 0 0 0 0 )\n( 2 0 5 0 0 0 )\n";
      throw std::runtime_error( oss.str() );
   }
   catch( std::out_of_range& ) {}
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the InitializerMatrix iterator implementation.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the iterator implementation of the InitializerMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testIterator()
{
   using blaze::initializer_list;


   using MatrixType    = blaze::InitializerMatrix<int>;
   using Iterator      = MatrixType::Iterator;
   using ConstIterator = MatrixType::ConstIterator;

   initializer_list< initializer_list<int> > list = { {  0,  1 },
                                                      { -2,  0, -3 },
                                                      {  0,  4,  5 } };

   MatrixType mat( list, 4UL );

   // Testing the Iterator default constructor
   {
      test_ = "Iterator default constructor";

      Iterator it{};

      if( it != Iterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing the ConstIterator default constructor
   {
      test_ = "ConstIterator default constructor";

      ConstIterator it{};

      if( it != ConstIterator() ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator default constructor\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing conversion from Iterator to ConstIterator
   {
      test_ = "Iterator/ConstIterator conversion";

      ConstIterator it( begin( mat, 1UL ) );

      if( it == end( mat, 1UL ) || *it != -2 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Failed iterator conversion detected\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in 0th row via Iterator (end-begin)
   {
      test_ = "Iterator subtraction (end-begin)";

      const ptrdiff_t number( end( mat, 0UL ) - begin( mat, 0UL ) );

      if( number != 4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in 0th row via Iterator (begin-end)
   {
      test_ = "Iterator subtraction (begin-end)";

      const ptrdiff_t number( begin( mat, 0UL ) - end( mat, 0UL ) );

      if( number != -4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: -4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in 1st row via ConstIterator (end-begin)
   {
      test_ = "ConstIterator subtraction (end-begin)";

      const ptrdiff_t number( cend( mat, 1UL ) - cbegin( mat, 1UL ) );

      if( number != 4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: 4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Counting the number of elements in 1st row via ConstIterator (begin-end)
   {
      test_ = "ConstIterator subtraction (begin-end)";

      const ptrdiff_t number( cbegin( mat, 1UL ) - cend( mat, 1UL ) );

      if( number != -4L ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid number of elements detected\n"
             << " Details:\n"
             << "   Number of elements         : " << number << "\n"
             << "   Expected number of elements: -4\n";
         throw std::runtime_error( oss.str() );
      }
   }

   // Testing read-only access via ConstIterator
   {
      test_ = "Read-only access via ConstIterator";

      ConstIterator it ( cbegin( mat, 2UL ) );
      ConstIterator end( cend( mat, 2UL ) );

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Invalid initial iterator detected\n";
         throw std::runtime_error( oss.str() );
      }

      ++it;

      if( it == end || *it != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      --it;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator pre-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it++;

      if( it == end || *it != 4 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-increment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it--;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator post-decrement failed\n";
         throw std::runtime_error( oss.str() );
      }

      it += 2UL;

      if( it == end || *it != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator addition assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it -= 2UL;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator subtraction assignment failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it + 2UL;

      if( it == end || *it != 5 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar addition failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = it - 2UL;

      if( it == end || *it != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Iterator/scalar subtraction failed\n";
         throw std::runtime_error( oss.str() );
      }

      it = 4UL + it;

      if( it != end ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Scalar/iterator addition failed\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c nonZeros() member function of the InitializerMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c nonZeros() member function of the InitializerMatrix class
// template. In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testNonZeros()
{
   using blaze::initializer_list;


   test_ = "InitializerMatrix::nonZeros()";

   {
      initializer_list< initializer_list<int> > list = { { 0, 0, 0 },
                                                         { 0, 0, 0 } };

      blaze::InitializerMatrix<int> mat( list );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 0UL );
      checkNonZeros( mat, 0UL, 0UL );
      checkNonZeros( mat, 1UL, 0UL );

      if( mat(0,0) != 0 || mat(0,1) != 0 || mat(0,2) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 0 || mat(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 0 0 )\n( 0 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      initializer_list< initializer_list<int> > list = { { 0, 1, 2 },
                                                         { 0, 3, 0 } };

      blaze::InitializerMatrix<int> mat( list );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 3UL );
      checkCapacity( mat, 6UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );

      if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 2 ||
          mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 1 2 )\n( 0 3 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }

   {
      initializer_list< initializer_list<int> > list = { { 0, 1, 2 },
                                                         { 0, 3, 0 } };

      blaze::InitializerMatrix<int> mat( list, 4UL );

      checkRows    ( mat, 2UL );
      checkColumns ( mat, 4UL );
      checkCapacity( mat, 8UL );
      checkNonZeros( mat, 3UL );
      checkNonZeros( mat, 0UL, 2UL );
      checkNonZeros( mat, 1UL, 1UL );

      if( mat(0,0) != 0 || mat(0,1) != 1 || mat(0,2) != 2 || mat(0,3) != 0 ||
          mat(1,0) != 0 || mat(1,1) != 3 || mat(1,2) != 0 || mat(1,3) != 0 ) {
         std::ostringstream oss;
         oss << " Test: " << test_ << "\n"
             << " Error: Initialization failed\n"
             << " Details:\n"
             << "   Result:\n" << mat << "\n"
             << "   Expected result:\n( 0 1 2 0 )\n( 0 3 0 0 )\n";
         throw std::runtime_error( oss.str() );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Test of the \c swap() functionality of the InitializerMatrix class template.
//
// \return void
// \exception std::runtime_error Error detected.
//
// This function performs a test of the \c swap() function of the InitializerMatrix class template.
// In case an error is detected, a \a std::runtime_error exception is thrown.
*/
void ClassTest::testSwap()
{
   using blaze::initializer_list;


   test_ = "InitializerMatrix swap";

   initializer_list< initializer_list<int> > list1{ { 1, 2 },
                                                    { 0, 3 },
                                                    { 4 } };
   initializer_list< initializer_list<int> > list2{ { 6, 5, 4 },
                                                    { 3, 2, 1 }  };

   blaze::InitializerMatrix<int> mat1( list1 );
   blaze::InitializerMatrix<int> mat2( list2, 4UL );

   swap( mat1, mat2 );

   checkRows    ( mat1, 2UL );
   checkColumns ( mat1, 4UL );
   checkCapacity( mat1, 8UL );
   checkNonZeros( mat1, 6UL );
   checkNonZeros( mat1, 0UL, 3UL );
   checkNonZeros( mat1, 1UL, 3UL );

   if( mat1(0,0) != 6 || mat1(0,1) != 5 || mat1(0,2) != 4 ||
       mat1(1,0) != 3 || mat1(1,1) != 2 || mat1(1,2) != 1 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the first matrix failed\n"
          << " Details:\n"
          << "   Result:\n" << mat1 << "\n"
          << "   Expected result:\n( 6 5 4 )\n( 3 2 1 )\n";
      throw std::runtime_error( oss.str() );
   }

   checkRows    ( mat2, 3UL );
   checkColumns ( mat2, 2UL );
   checkCapacity( mat2, 6UL );
   checkNonZeros( mat2, 4UL );
   checkNonZeros( mat2, 0UL, 2UL );
   checkNonZeros( mat2, 1UL, 1UL );
   checkNonZeros( mat2, 2UL, 1UL );

   if( mat2(0,0) != 1 || mat2(0,1) != 2 ||
       mat2(1,0) != 0 || mat2(1,1) != 3 ||
       mat2(2,0) != 4 || mat2(2,1) != 0 ) {
      std::ostringstream oss;
      oss << " Test: " << test_ << "\n"
          << " Error: Swapping the second matrix failed\n"
          << " Details:\n"
          << "   Result:\n" << mat2 << "\n"
          << "   Expected result:\n( 1 2 )\n( 0 3 )\n( 4, 0 )\n";
      throw std::runtime_error( oss.str() );
   }
}
//*************************************************************************************************

} // namespace initializermatrix

} // namespace matrices

} // namespace mathtest

} // namespace blazetest




//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running InitializerMatrix class test..." << std::endl;

   try
   {
      RUN_INITIALIZERMATRIX_CLASS_TEST;
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during InitializerMatrix class test:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
