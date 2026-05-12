//=================================================================================================
/*!
//  \file blazemark/util/Parser.h
//  \brief Header file for the Parser class template
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

#ifndef _BLAZEMARK_UTIL_PARSER_H_
#define _BLAZEMARK_UTIL_PARSER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <blazemark/system/Types.h>


namespace blazemark {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Benchmark run extractor.
//
// The Parser class provides the functionality to extract the parameters for benchmark runs
// from a given parameter file. Via the template parameter it is possible to configure the
// parser to benchmark runs for either dense or sparse vectors and/or matrices.
//
//
// \section dense Benchmark runs for dense vectors/matrices
//
// For benchmarks for dense vectors and/or matrices, the Parser has to be instantiated with
// blazemark::DenseRun:

   \code
   blazemark::Parser< blazemark::DenseRun > parser;
   std::vector< blazemark::DenseRun > runs;
   parser.parse( "dvecdvecadd.prm", runs );
   \endcode

// In the example, 'dvecdvecadd.prm' is the name of the parameter file (in this case for
// the dense vector/dense vector addition), which contains the parameters for the individual
// benchmark runs. The following example shows a parameter file for benchmarks for dense
// vectors/matrices:

   \code
   (   100, 50000 )
   (  1000,  5000 )
   ( 10000        )
   \endcode

// In this example, three separate benchmark runs are defined. The first value specifies
// the target size of the dense vectors or the number of rows and columns of the matrices,
// respectively. The second value, which is optional, specifies the number of steps that
// should be performed within the benchmark. In case the second value is omitted, the
// necessary number of steps will automatically be evaluated.
//
//
// \section sparse Benchmark runs for sparse vectors/matrices
//
// For benchmarks involving sparse vectors and/or matrices, the Parser has to be instantiated
// with blazemark::SparseRun:

   \code
   blazemark::Parser< blazemark::SparseRun > parser;
   std::vector< blazemark::SparseRun > runs;
   parser.parse( "svecsvecadd.prm", runs );
   \endcode

// In the example, 'svecsvecadd.prm' is the name of the parameter file (in this case for
// the sparse vector/sparse vector addition), which contains the parameters for the individual
// benchmark runs. The following example shows a parameter file for benchmarks for sparse
// vectors/matrices:

   \code
   (  100, 10, 50000 )
   (  100, 40, 10000 )
   ( 1000, 10        )
   ( 1000, 40        )
   \endcode

// In this particular example, a total of four benchmark runs is defined. The first value
// specifies the target size of the sparse vectors or the number of rows and columns of the
// matrices, respectively. The second value specifies the number of non-zero elements in
// the sparse vectors or in case of matrices the number of non-zero elements per row. The
// third value, which is optional, specifies the number of steps that should be performed
// within the benchmark. In case the third value is omitted, the necessary number of steps
// will automatically be evaluated.
*/
template< typename RunType >  // Type of the benchmark run
class Parser
{
 private:
   //**Type definitions****************************************************************************
   using StreamPos  = std::stringstream::pos_type;  //!< Stream position.
   using Pair       = std::pair<StreamPos,size_t>;  //!< Pair of a stream position and a line.
   using LineVector = std::vector<Pair>;            //!< Vector of position/line pairs.
   //**********************************************************************************************

 public:
   //**Constructors********************************************************************************
   // No explicitly declared constructors.
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   // No explicitly declared destructor.
   //**********************************************************************************************

   //**Copy assignment operator********************************************************************
   // No explicitly declared copy assignment operator.
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   void parse( const char* const filename, std::vector<RunType>& runs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   size_t getLineNumber();
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   LineVector lineNumbers_;     //!< Vector for the line numbers of the input file.
   std::stringstream input_;    //!< Buffer for the preprocessed benchmark run parameters.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Extracting all benchmark run parameters from a given parameter file.
//
// \param filename The name of the parameter file.
// \param runs Vector to be filled with the extracted benchmark runs.
// \return void
// \exception std::invalid_argument Could not open parameter file.
//
// This function extracts all benchmark runs of the specified type from the given parameter
// file. All previously extracted benchmark runs are deleted.
*/
template< typename RunType >  // Type of the benchmark run
void Parser<RunType>::parse( const char* const filename, std::vector<RunType>& runs )
{
   ///////////////////////////////////////////////////////
   // Resetting the input stream and line number vector

   runs.clear();
   lineNumbers_.clear();
   input_.str( "" );
   input_.clear();


   ////////////////////////////////////////////////////////
   // Preprocessing the rigid body parameter file stream

   std::ifstream in( filename, std::ifstream::in );
   if( !in.is_open() ) {
      throw std::invalid_argument( "Could not open parameter file" );
   }

   bool comment( false );
   size_t lineCounter( 0 );
   std::string::size_type pos1, pos2;
   std::string line;

   line.reserve( 100 );

   while( std::getline( in, line ) )
   {
      ++lineCounter;

      if( comment ) {
         if( ( pos1 = line.find( "*/", 0 ) ) != std::string::npos ) {
            line.erase( line.begin(), line.begin()+pos1+2 );
            comment = false;
         }
         else continue;
      }

      if( ( pos1 = line.find( "//", 0 ) ) != std::string::npos ) {
         line.erase( line.begin()+pos1, line.end() );
      }
      if( ( pos1 = line.find( "/*", 0 ) ) != std::string::npos ) {
         if( ( pos2 = line.find( "*/", pos1+2 ) ) != std::string::npos ) {
            line.replace( line.begin()+pos1, line.begin()+pos2+2, " " );
         }
         else {
            line.erase( line.begin()+pos1, line.end() );
            comment = true;
         }
      }

      // Adding the line to the input string
      lineNumbers_.push_back( Pair( input_.tellp(), lineCounter ) );
      input_ << line << "\n";
   }

   in.close();


   ///////////////////////////////////
   // Extracting the benchmark runs

   RunType run;

   // Extracting all benchmark run parameters
   while( (input_ >> std::ws) && input_.peek() != std::stringstream::traits_type::eof() )
   {
      // Extracting the parameters for a single benchmark run
      if( !(input_ >> run) ) {
         input_.clear();
         std::ostringstream oss;
         oss << "Input error in line " << getLineNumber() << " of parameter file '" << filename << "'";
         throw std::logic_error( oss.str() );
      }

      // Storing the benchmark run
      runs.push_back( run );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Estimating the input file line depending on the stream position of the preprocessed stream.
//
// \return The line number of the input file.
*/
template< typename RunType >  // Type of the benchmark run
size_t Parser<RunType>::getLineNumber()
{
   // Estimating the stream position
   const StreamPos pos( input_.tellg() );

   // Calculating the line number in the input file
   LineVector::const_iterator it=lineNumbers_.begin();
   for( ; it!=lineNumbers_.end()-1; ++it ) {
      if( (it+1)->first > pos )
         return it->second;
   }
   return it->second;
}
//*************************************************************************************************

} // namespace blazemark

#endif
