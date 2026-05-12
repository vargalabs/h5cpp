// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2015 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2015 National ICT Australia (NICTA)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------

#include <armadillo>
#include "catch.hpp"

using namespace arma;


TEST_CASE("fn_solve_1")
  {
  // square-sized A
  
  mat A =
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768;\
    ";
  
  mat B = 
    "\
     0.906104;\
     0.372084;\
     0.826970;\
     0.911125;\
     0.764937;\
    ";
  
  mat X1;
  mat X2;
  
  bool status1 = solve(X1, A, B);
  bool status2 = solve(X2, A, B, solve_opts::fast + solve_opts::no_approx);
  
  REQUIRE( status1 == true );
  REQUIRE( status2 == true );
  
  REQUIRE( norm(vectorise(X1-X2),1) == Approx(0.0) );
  
  REQUIRE( X1(0) == Approx(-1.486404584527645e+00) );
  REQUIRE( X1(1) == Approx(-1.298396730449162e+01) );
  REQUIRE( X1(2) == Approx( 9.470697159029561e+00) );
  REQUIRE( X1(3) == Approx(-9.356796987754391e+00) );
  REQUIRE( X1(4) == Approx( 9.375696966670487e+00) );
  }



TEST_CASE("fn_solve_2")
  {
  // square-sized A; rank-deficient
  
  mat A =
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768;\
    ";
  
  A.col(0) *= datum::eps;
  
  mat B = 
    "\
     0.906104;\
     0.372084;\
     0.826970;\
     0.911125;\
     0.764937;\
    ";
  
  mat X1;
  mat X2;
  
  bool status1 = solve(X1, A, B);
  bool status2 = solve(X2, A, B, solve_opts::no_approx);
  
  REQUIRE( status1 == true  );
  REQUIRE( status2 == false );
  
  REQUIRE( X1(0) == Approx( 9.929054488765161e-15) );
  REQUIRE( X1(1) == Approx(-1.088301076771930e+01) );
  REQUIRE( X1(2) == Approx( 8.435894619970558e+00) );
  REQUIRE( X1(3) == Approx(-6.820665349768506e+00) );
  REQUIRE( X1(4) == Approx( 6.785813315231441e+00) );
  }



TEST_CASE("fn_solve_3")
  {
  // non-square-sized A
  
  mat A =
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745   0.051408;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153   0.035437;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317  -0.454499;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040   0.373833;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768   0.258704;\
    ";
  
  mat B = 
    "\
     0.906104;\
     0.372084;\
     0.826970;\
     0.911125;\
     0.764937;\
    ";
  
  mat X1;
  mat X2;
  
  bool status1 = solve(X1, A, B);
  bool status2 = solve(X2, A, B, solve_opts::fast + solve_opts::no_approx);
  
  REQUIRE( status1 == true );
  REQUIRE( status2 == true );
  
  REQUIRE( norm(vectorise(X1-X2),1) == Approx(0.0) );
  
  REQUIRE( X1(0) == Approx( 4.115847103194906) );
  REQUIRE( X1(1) == Approx(-3.878220128389104) );
  REQUIRE( X1(2) == Approx( 4.328088888881608) );
  REQUIRE( X1(3) == Approx(-2.949916817226818) );
  REQUIRE( X1(4) == Approx(-1.602040603621000) );
  REQUIRE( X1(5) == Approx(-5.985543296434588) );
  }
