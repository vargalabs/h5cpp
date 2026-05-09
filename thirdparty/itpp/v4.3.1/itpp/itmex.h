/*!
 * \file
 * \brief Conversion routines between IT++ and Matlab
 * \author Tony Ottosson and Pal Frenger
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

#ifndef ITMEX_H
#define ITMEX_H

#include <itpp/itbase.h>
#include <mex.h>


namespace itpp
{

//--------------------------------------------------------
// mex -> it++
//--------------------------------------------------------

/*!
  \addtogroup mexfiles
  \brief Conversions between IT++ and Matlab for writing mex-files.
  \author Tony Ottosson and Pal Frenger

  These routines are used to help writng mex-files for speeding up
  matlab simulations.
  A simple mex-file that performs QPSK modulation is given below.

  \code
  #include <itpp/itcomm.h>
  #include <itpp/itmex.h>

  using namespace itpp;

  void mexFunction(int n_output, mxArray *output[], int n_input, const mxArray *input[])
  {
  // Check the number of inputs and output arguments
  if(n_output!=1) mexErrMsgTxt("Wrong number of output variables!");
  if(n_input!=1) mexErrMsgTxt("Wrong number of input variables!");

  // Convert input variables to IT++ format
  bvec input_bits = mxArray2bvec(input[0]);

  // ------------------ Start of routine ---------------------------
  cvec output_symbols;
  QPSK qpsk;

  output_symbols = qpsk.modulate_bits(input_bits);
  // ------------------ End of routine -----------------------------

  // Create output vectors
  output[0] = mxCreateDoubleMatrix(1,output_symbols.size(), mxCOMPLEX);

  // Convert the IT++ format to Matlab format for output
  cvec2mxArray(output_symbols, output[0]);
  }
  \endcode
*/

/*!
 * \addtogroup mexfiles
 * @{
 */

// --------------------------------------------------------
// mex -> IT++
// --------------------------------------------------------

//! Convert the matlab-format mxArray to bin
bin mxArray2bin(const mxArray *in);
//! Convert the matlab-format mxArray to short
short mxArray2short(const mxArray *in);
//! Convert the matlab-format mxArray to int
int mxArray2int(const mxArray *in);
//! Convert the matlab-format mxArray to double
double mxArray2double(const mxArray *in);
//! Convert the matlab-format mxArray to complex<double>
std::complex<double> mxArray2double_complex(const mxArray *in);
//! Convert the matlab-format mxArray to string
std::string mxArray2string(const mxArray *in);

//! Convert the matlab-format mxArray to bvec
bvec mxArray2bvec(const mxArray *in);
//! Convert the matlab-format mxArray to svec
svec mxArray2svec(const mxArray *in);
//! Convert the matlab-format mxArray to ivec
ivec mxArray2ivec(const mxArray *in);
//! Convert the matlab-format mxArray to vec
vec mxArray2vec(const mxArray *in);
//! Convert the matlab-format mxArray to cvec
cvec mxArray2cvec(const mxArray *in);

//! Convert the matlab-format mxArray to bmat
bmat mxArray2bmat(const mxArray *in);
//! Convert the matlab-format mxArray to smat
smat mxArray2smat(const mxArray *in);
//! Convert the matlab-format mxArray to imat
imat mxArray2imat(const mxArray *in);
//! Convert the matlab-format mxArray to mat
mat mxArray2mat(const mxArray *in);
//! Convert the matlab-format mxArray to cmat
cmat mxArray2cmat(const mxArray *in);

// --------------------------------------------------------
// IT++ -> mex
// --------------------------------------------------------

//! Convert bin to the matlab-format mxArray
void bin2mxArray(const bin &in, mxArray *out);
//! Convert short to the matlab-format mxArray
void short2mxArray(const short &in, mxArray *out);
//! Convert int to the matlab-format mxArray
void int2mxArray(const int &in, mxArray *out);
//! Convert double to the matlab-format mxArray
void double2mxArray(const double &in, mxArray *out);
//! Convert complex<double> to the matlab-format mxArray
void double_complex2mxArray(const std::complex<double> &in, mxArray *out);
//! Convert string to the matlab-format mxArray
void string2mxArray(const std::string &in, mxArray* &out);

//! Convert bvec to the matlab-format mxArray
void bvec2mxArray(const bvec &in, mxArray *out);
//! Convert svec to the matlab-format mxArray
void svec2mxArray(const svec &in, mxArray *out);
//! Convert ivec to the matlab-format mxArray
void ivec2mxArray(const ivec &in, mxArray *out);
//! Convert vec to the matlab-format mxArray
void vec2mxArray(const vec &in, mxArray *out);
//! Convert cvec to the matlab-format mxArray
void cvec2mxArray(const cvec &in, mxArray *out);

//! Convert bmat to the matlab-format mxArray
void bmat2mxArray(const bmat &in, mxArray *out);
//! Convert smat to the matlab-format mxArray
void smat2mxArray(const smat &in, mxArray *out);
//! Convert imat to the matlab-format mxArray
void imat2mxArray(const imat &in, mxArray *out);
//! Convert mat to the matlab-format mxArray
void mat2mxArray(const mat &in, mxArray *out);
//! Convert cmat to the matlab-format mxArray
void cmat2mxArray(const cmat &in, mxArray *out);

// --------------------------------------------------------
// mex -> C
// --------------------------------------------------------

//! Convert the matlab-format mxArray to C-format pointer to short
void mxArray2Csvec(const mxArray *in, short *out);
//! Convert the matlab-format mxArray to C-format pointer to int
void mxArray2Civec(const mxArray *in, int *out);
//! Convert the matlab-format mxArray to C-format pointer to double
void mxArray2Cvec(const mxArray *in, double *out);
//! Convert the matlab-format mxArray to C-format pointers to double (real and imaginary parts)
void mxArray2Ccvec(const mxArray *in, double *out_real, double *out_imag);

//! Convert the matlab-format mxArray to C-format pointer to pointer to short
void mxArray2Csmat(const mxArray *in, short **out);
//! Convert the matlab-format mxArray to C-format pointer to pointer to int
void mxArray2Cimat(const mxArray *in, int **out);
//! Convert the matlab-format mxArray to C-format pointer to pointer to double
void mxArray2Cmat(const mxArray *in, double **out);
//! Convert the matlab-format mxArray to C-format pointer to pointer to double (real and imaginary parts)
void mxArray2Ccmat(const mxArray *in, double **out_real, double **out_imag);

// --------------------------------------------------------
// C -> mex
// --------------------------------------------------------

//! Convert C-format pointer to short to matlab-format mxArray
void Csvec2mxArray(short *in, mxArray *out);
//! Convert C-format pointer to int to matlab-format mxArray
void Civec2mxArray(int *in, mxArray *out);
//! Convert C-format pointer to double to matlab-format mxArray
void Cvec2mxArray(double *in, mxArray *out);
//! Convert C-format pointers to double (real and imaginary parts) to matlab-format mxArray
void Ccvec2mxArray(double *in_real, double *in_imag, mxArray *out);

//! Convert C-format pointer to pointer to short to matlab-format mxArray
void Csmat2mxArray(short **in, mxArray *out);
//! Convert C-format pointer to pointer to int to matlab-format mxArray
void Cimat2mxArray(int **in, mxArray *out);
//! Convert C-format pointer to pointer to double to matlab-format mxArray
void Cmat2mxArray(double **in, mxArray *out);
//! Convert C-format pointer to pointer to double (real and imaginary parts) to matlab-format mxArray
void Ccmat2mxArray(double **in_real, double **in_imag, mxArray *out);

/*!
 * @}
 */


bin mxArray2bin(const mxArray *in)
{
  int size;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2bin: Pointer to data is NULL");
  size = mxGetNumberOfElements(in);
  if (size != 1) mexErrMsgTxt("mxArray2bin: Size of data is not equal to one");

  return (((*temp) > 0.0) ? bin(1) : bin(0));
}

short mxArray2short(const mxArray *in)
{
  int size;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2short: Pointer to data is NULL");
  size = mxGetNumberOfElements(in);
  if (size != 1) mexErrMsgTxt("mxArray2short: Size of data is not equal to one");

  return (short)(*temp);
}

int mxArray2int(const mxArray *in)
{
  int size;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2int: Pointer to data is NULL");
  size = mxGetNumberOfElements(in);
  if (size != 1) mexErrMsgTxt("mxArray2int: Size of data is not equal to one");

  return (int)(*temp);
}

double mxArray2double(const mxArray *in)
{
  int size;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2double: Pointer to data is NULL");
  size = mxGetNumberOfElements(in);
  if (size != 1) mexErrMsgTxt("mxArray2double: Size of data is not equal to one");

  return (*temp);
}

std::complex<double> mxArray2double_complex(const mxArray *in)
{
  int size;
  double* tempR = (double*) mxGetPr(in);
  double* tempI = (double*) mxGetPi(in);

  if ((tempR == 0) && (tempI == 0)) mexErrMsgTxt("mxArray2double_complex: Pointer to data is NULL");

  size = mxGetNumberOfElements(in);
  if (size != 1) mexErrMsgTxt("mxArray2double_complex: Size of data is not equal to one");

  if (tempR == 0) {
    return std::complex<double>(0.0 , (*tempI));
  }
  else if (tempI == 0) {
    return std::complex<double>((*tempR), 0.0);
  }
  else {
    return std::complex<double>((*tempR), (*tempI));
  }

}

std::string mxArray2string(const mxArray *in)
{
  if (in == 0)
    mexErrMsgTxt("mxArray2string: Pointer to data is NULL");
  std::string str = mxArrayToString(in);
  if (str.data() == 0)
    mexErrMsgTxt("mxArray2string: Could not convert mxArray to string");
  return str;
}

bvec mxArray2bvec(const mxArray *in)
{
  bvec out;
  int i, size;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2bvec: Pointer to data is NULL");

  size = mxGetNumberOfElements(in);
  if (size == 0) mexErrMsgTxt("mxArray2bvec: Size of data is zero");

  out.set_size(size, false);

  for (i = 0; i < size; i++) {
    out(i) = (((*temp++) > 1e-5) ? bin(1) : bin(0));
  }

  return out;

}

svec mxArray2svec(const mxArray *in)
{
  svec out;
  int i, size;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2svec: Pointer to data is NULL");

  size = mxGetNumberOfElements(in);
  if (size == 0) mexErrMsgTxt("mxArray2svec: Size of data is zero");

  out.set_size(size, false);

  for (i = 0; i < size; i++) {
    out(i) = (short)(*temp++);
  }

  return out;

}

ivec mxArray2ivec(const mxArray *in)
{
  ivec out;
  int i, size;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2ivec: Pointer to data is NULL");

  size = mxGetNumberOfElements(in);
  if (size == 0) mexErrMsgTxt("mxArray2ivec: Size of data is zero");

  out.set_size(size, false);

  for (i = 0; i < size; i++) {
    out(i) = (int)(*temp++);
  }

  return out;

}

vec mxArray2vec(const mxArray *in)
{
  vec out;
  int i, size;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2vec: Pointer to data is NULL");

  size = mxGetNumberOfElements(in);
  if (size == 0) mexErrMsgTxt("mxArray2vec: Size of data is zero");

  out.set_size(size, false);

  for (i = 0; i < size; i++) {
    out(i) = (*temp++);
  }

  return out;

}

cvec mxArray2cvec(const mxArray *in)
{
  cvec out;
  int i, size;
  double* tempR = (double*) mxGetPr(in);
  double* tempI = (double*) mxGetPi(in);

  if ((tempR == 0) && (tempI == 0)) mexErrMsgTxt("mxArray2cvec: Pointer data is NULL");

  size = mxGetNumberOfElements(in);
  if (size == 0) mexErrMsgTxt("mxArray2cvec: Size of data is zero");

  out.set_size(size, false);

  if (tempR == 0) {
    for (i = 0; i < size; i++) { out(i) = std::complex<double>(0.0, (*tempI++)); }
  }
  else if (tempI == 0) {
    for (i = 0; i < size; i++) { out(i) = std::complex<double>((*tempR++), 0.0); }
  }
  else {
    for (i = 0; i < size; i++) { out(i) = std::complex<double>((*tempR++), (*tempI++)); }
  }

  return out;

}

bmat mxArray2bmat(const mxArray *in)
{
  bmat out;
  int r, c, rows, cols;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2bmat: Pointer to data is NULL");

  rows = mxGetM(in);
  if (rows == 0) mexErrMsgTxt("mxArray2bmat: Data has zero rows");
  cols = mxGetN(in);
  if (cols == 0) mexErrMsgTxt("mxArray2bmat: Data has zero columns");

  out.set_size(rows, cols, false);

  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      out(r, c) = (((*temp++) > 0.0) ? bin(1) : bin(0));
    }
  }

  return out;

}

smat mxArray2smat(const mxArray *in)
{
  smat out;
  int r, c, rows, cols;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2smat: Pointer to data is NULL");

  rows = mxGetM(in);
  if (rows == 0) mexErrMsgTxt("mxArray2smat: Data has zero rows");
  cols = mxGetN(in);
  if (cols == 0) mexErrMsgTxt("mxArray2smat: Data has zero columns");

  out.set_size(rows, cols, false);

  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      out(r, c) = (short)(*temp++);
    }
  }

  return out;

}

imat mxArray2imat(const mxArray *in)
{
  imat out;
  int r, c, rows, cols;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2imat: Pointer to data is NULL");

  rows = mxGetM(in);
  if (rows == 0) mexErrMsgTxt("mxArray2imat: Data has zero rows");
  cols = mxGetN(in);
  if (cols == 0) mexErrMsgTxt("mxArray2imat: Data has zero columns");
  out.set_size(rows, cols, false);

  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      out(r, c) = (int)(*temp++);
    }
  }

  return out;

}

mat mxArray2mat(const mxArray *in)
{
  mat out;
  int r, c, rows, cols;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2mat: Pointer to data is NULL");

  rows = mxGetM(in);
  if (rows == 0) mexErrMsgTxt("mxArray2mat: Data has zero rows");
  cols = mxGetN(in);
  if (cols == 0) mexErrMsgTxt("mxArray2mat: Data has zero columns");
  out.set_size(rows, cols, false);

  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      out(r, c) = (*temp++);
    }
  }

  return out;

}

cmat mxArray2cmat(const mxArray *in)
{
  cmat out;
  int r, c, rows, cols;
  double* tempR = (double*) mxGetPr(in);
  double* tempI = (double*) mxGetPi(in);

  if ((tempR == 0) && (tempI == 0)) mexErrMsgTxt("mxArray2cmat: Pointer to data is NULL");

  rows = mxGetM(in);
  if (rows == 0) mexErrMsgTxt("mxArray2cmat: Data has zero rows");
  cols = mxGetN(in);
  if (cols == 0) mexErrMsgTxt("mxArray2cmat: Data has zero columns");
  out.set_size(rows, cols, false);

  if (tempR == 0) {
    for (c = 0; c < cols; c++) { for (r = 0; r < rows; r++) { out(r, c) = std::complex<double>(0.0 , (*tempI++)); } }
  }
  else if (tempI == 0) {
    for (c = 0; c < cols; c++) { for (r = 0; r < rows; r++) { out(r, c) = std::complex<double>((*tempR++), 0.0); } }
  }
  else {
    for (c = 0; c < cols; c++) { for (r = 0; r < rows; r++) { out(r, c) = std::complex<double>((*tempR++), (*tempI++)); } }
  }

  return out;

}

void double2mxArray(const double &in, mxArray *out)
{
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("double2mxArray: Pointer to data is NULL");

  *temp = (double) in;
}

void double_complex2mxArray(const std::complex<double> &in, mxArray *out)
{
  double* tempR = (double *) mxGetPr(out);
  double* tempI = (double *) mxGetPi(out);
  if (tempR == 0) mexErrMsgTxt("double_complex2mxArray: Pointer to real valued part is NULL");
  if (tempI == 0) mexErrMsgTxt("double_complex2mxArray: Pointer to imaginary valued part is NULL");

  *tempR = (double) in.real();
  *tempI = (double) in.imag();
}

void string2mxArray(const std::string &in, mxArray* &out)
{
  if (in.data() == 0)
    mexErrMsgTxt("string2mxArray: Pointer to string is NULL");
  out = mxCreateString(in.data());
  if (out == 0)
    mexErrMsgTxt("string2mxArray: Could not convert string to mxArray");
}

void bvec2mxArray(const bvec &in, mxArray *out)
{
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("bvec2mxArray: Pointer to data is NULL");
  if (in.size() == 0) mexErrMsgTxt("bvec2mxArray: Size of data is zero");
  for (int i = 0; i < in.size(); i++) {
    if (in(i))
      *temp++ = 1.0;
    else
      *temp++ = 0.0;
  }
}

void ivec2mxArray(const ivec &in, mxArray *out)
{
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("ivec2mxArray: Pointer to data is NULL");
  if (in.size() == 0) mexErrMsgTxt("ivec2mxArray: Size of data is zero");

  for (int i = 0; i < in.size(); i++) {
    *temp++ = (double) in(i);
  }
}

void vec2mxArray(const vec &in, mxArray *out)
{
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("vec2mxArray: Pointer to data is NULL");
  if (in.size() == 0) mexErrMsgTxt("vec2mxArray: Size of data is zero");

  for (int i = 0; i < in.size(); i++) {
    *temp++ = (double) in(i);
  }
}

void cvec2mxArray(const cvec &in, mxArray *out)
{
  double* tempR = (double *) mxGetPr(out);
  double* tempI = (double *) mxGetPi(out);
  if (tempR == 0) mexErrMsgTxt("cvec2mxArray: Pointer to real valued part is NULL");
  if (tempI == 0) mexErrMsgTxt("cvec2mxArray: Pointer to imaginary valued part is NULL");
  if (in.size() == 0) mexErrMsgTxt("cvec2mxArray: Size of data is zero");

  for (int i = 0; i < in.size(); i++) {
    *tempR++ = (double) in(i).real();
    *tempI++ = (double) in(i).imag();
  }
}

void bmat2mxArray(const bmat &in, mxArray *out)
{
  int rows, cols, r, c;

  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("bmat2mxArray: Pointer to data is NULL");

  rows = in.rows();
  cols = in.cols();
  if (rows == 0) mexErrMsgTxt("bmat2mxArray: Data has zero rows");
  if (cols == 0) mexErrMsgTxt("bmat2mxArray: Data has zero columns");

  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      if (in(r, c))
        *temp++ = 1.0;
      else
        *temp++ = 0.0;
    }
  }

}

void smat2mxArray(const smat &in, mxArray *out)
{
  int rows, cols, r, c;

  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("smat2mxArray: Pointer to data is NULL");

  rows = in.rows();
  cols = in.cols();
  if (rows == 0) mexErrMsgTxt("smat2mxArray: Data has zero rows");
  if (cols == 0) mexErrMsgTxt("smat2mxArray: Data has zero columns");

  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      *temp++ = (double) in(r, c);
    }
  }

}

void imat2mxArray(const imat &in, mxArray *out)
{
  int rows, cols, r, c;

  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("imat2mxArray: Pointer to data is NULL");

  rows = in.rows();
  cols = in.cols();
  if (rows == 0) mexErrMsgTxt("imat2mxArray: Data has zero rows");
  if (cols == 0) mexErrMsgTxt("imat2mxArray: Data has zero columns");

  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      *temp++ = (double) in(r, c);
    }
  }

}

void mat2mxArray(const mat &in, mxArray *out)
{
  int rows, cols, r, c;

  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("mat2mxArray: Pointer to data is NULL");

  rows = in.rows();
  cols = in.cols();
  if (rows == 0) mexErrMsgTxt("mat2mxArray: Data has zero rows");
  if (cols == 0) mexErrMsgTxt("mat2mxArray: Data has zero columns");

  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      *temp++ = in(r, c);
    }
  }

}

void cmat2mxArray(const cmat &in, mxArray *out)
{
  int rows, cols, r, c;

  double* tempR = (double *) mxGetPr(out);
  double* tempI = (double *) mxGetPi(out);
  if (tempR == 0) mexErrMsgTxt("cvec2mxArray: Pointer to real valued part is NULL");
  if (tempI == 0) mexErrMsgTxt("cvec2mxArray: Pointer to imaginary valued part is NULL");

  rows = in.rows();
  cols = in.cols();
  if (rows == 0) mexErrMsgTxt("cvec2mxArray: Data has zero rows");
  if (cols == 0) mexErrMsgTxt("cvec2mxArray: Data has zero columns");

  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      *tempR++ = (double) in(r, c).real();
      *tempI++ = (double) in(r, c).imag();
    }
  }

}

void mxArray2Csvec(const mxArray *in, short *out)
{
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2Csvec: Pointer to data is NULL");
  int size = mxGetNumberOfElements(in);
  if (size == 0) mexErrMsgTxt("mxArray2Csvec: Size of data is zero");
  for (int i = 0; i < size; i++) { out[i] = (short)(*temp++); }
}

void mxArray2Civec(const mxArray *in, int *out)
{
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2Civec: Pointer to data is NULL");
  int size = mxGetNumberOfElements(in);
  if (size == 0) mexErrMsgTxt("mxArray2Civec: Size of data is zero");
  for (int i = 0; i < size; i++) { out[i] = (int)(*temp++); }
}

void mxArray2Cvec(const mxArray *in, double *out)
{
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2Cvec: Pointer to data is NULL");
  int size = mxGetNumberOfElements(in);
  if (size == 0) mexErrMsgTxt("mxArray2Cvec: Size of data is zero");
  for (int i = 0; i < size; i++) { out[i] = (*temp++); }
}

void mxArray2Ccvec(const mxArray *in, double *out_real, double *out_imag)
{
  double* tempR = (double*) mxGetPr(in);
  double* tempI = (double*) mxGetPi(in);
  if (tempR == 0) mexErrMsgTxt("mxArray2Ccvec: Pointer to real valued part is NULL");
  if (tempI == 0) mexErrMsgTxt("mxArray2Ccvec: Pointer to imaginary valued part is NULL");
  int size = mxGetNumberOfElements(in);
  if (size == 0) mexErrMsgTxt("mxArray2Ccvec: Size of data is zero");
  for (int i = 0; i < size; i++) { out_real[i] = (*tempR++); out_imag[i] = (*tempI++); }
}

void mxArray2Csmat(const mxArray *in, short **out)
{
  int r, c;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2Csmat: Pointer to data is NULL");
  int rows = mxGetM(in);
  if (rows == 0) mexErrMsgTxt("mxArray2Csmat: Data has zero rows");
  int cols = mxGetN(in);
  if (cols == 0) mexErrMsgTxt("mxArray2Csmat: Data has zero columns");
  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      out[r][c] = (short)(*temp++);
    }
  }
}

void mxArray2Cimat(const mxArray *in, int **out)
{
  int r, c;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2Cimat: Pointer to data is NULL");
  int rows = mxGetM(in);
  if (rows == 0) mexErrMsgTxt("mxArray2Cimat: Data has zero rows");
  int cols = mxGetN(in);
  if (cols == 0) mexErrMsgTxt("mxArray2Cimat: Data has zero columns");
  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      out[r][c] = (int)(*temp++);
    }
  }
}

void mxArray2Cmat(const mxArray *in, double **out)
{
  int r, c;
  double* temp = (double*) mxGetPr(in);
  if (temp == 0) mexErrMsgTxt("mxArray2Cmat: Pointer to data is NULL");
  int rows = mxGetM(in);
  if (rows == 0) mexErrMsgTxt("mxArray2Cmat: Data has zero rows");
  int cols = mxGetN(in);
  if (cols == 0) mexErrMsgTxt("mxArray2Cmat: Data has zero columns");
  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      out[r][c] = (*temp++);
    }
  }
}

void mxArray2Ccmat(const mxArray *in, double **out_real, double **out_imag)
{
  int r, c;
  double* tempR = (double*) mxGetPr(in);
  double* tempI = (double*) mxGetPi(in);
  if (tempR == 0) mexErrMsgTxt("mxArray2Cmat: Pointer to real valued part is NULL");
  if (tempI == 0) mexErrMsgTxt("mxArray2Cmat: Pointer to imaginary valued part is NULL");
  int rows = mxGetM(in);
  if (rows == 0) mexErrMsgTxt("mxArray2Cmat: Data has zero rows");
  int cols = mxGetN(in);
  if (cols == 0) mexErrMsgTxt("mxArray2Cmat: Data has zero columns");
  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      out_real[r][c] = (*tempR++);
      out_imag[r][c] = (*tempI++);
    }
  }
}

void Csvec2mxArray(short *in, mxArray *out)
{
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("Csvec2mxArray: Pointer to data is NULL");
  int size = mxGetNumberOfElements(out);
  if (size == 0) mexErrMsgTxt("Csvec2mxArray: Size of data is zero");
  for (int i = 0; i < size; i++) { *temp++ = (double) in[i]; }
}

void Civec2mxArray(int *in, mxArray *out)
{
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("Civec2mxArray: Pointer to data is NULL");
  int size = mxGetNumberOfElements(out);
  if (size == 0) mexErrMsgTxt("Civec2mxArray: Size of data is zero");
  for (int i = 0; i < size; i++) { *temp++ = (double) in[i]; }
}

void Cvec2mxArray(double *in, mxArray *out)
{
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("Cvec2mxArray: Pointer to data is NULL");
  int size = mxGetNumberOfElements(out);
  if (size == 0) mexErrMsgTxt("Cvec2mxArray: Size of data is zero");
  for (int i = 0; i < size; i++) { *temp++ = in[i]; }
}

void Ccvec2mxArray(double *in_real, double *in_imag, mxArray *out)
{
  double* tempR = (double *) mxGetPr(out);
  double* tempI = (double *) mxGetPi(out);
  if (tempR == 0) mexErrMsgTxt("Ccvec2mxArray: Pointer to real valued part is NULL");
  if (tempI == 0) mexErrMsgTxt("Ccvec2mxArray: Pointer to imaginary valued part is NULL");
  int size = mxGetNumberOfElements(out);
  if (size == 0) mexErrMsgTxt("Ccvec2mxArray: Size of data is zero");
  for (int i = 0; i < size; i++) { *tempR++ = in_real[i]; *tempI++ = in_imag[i]; }
}

void Csmat2mxArray(short **in, mxArray *out)
{
  int r, c;
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("Csmat2mxArray: Pointer to data is NULL");
  int rows = mxGetM(out);
  if (rows == 0) mexErrMsgTxt("Csmat2mxArray: Data has zero rows");
  int cols = mxGetN(out);
  if (cols == 0) mexErrMsgTxt("Csmat2mxArray: Data has zero columns");
  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      *temp++ = (short) in[r][c];
    }
  }
}

void Cimat2mxArray(int **in, mxArray *out)
{
  int r, c;
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("Cimat2mxArray: Pointer to data is NULL");
  int rows = mxGetM(out);
  if (rows == 0) mexErrMsgTxt("Cimat2mxArray: Data has zero rows");
  int cols = mxGetN(out);
  if (cols == 0) mexErrMsgTxt("Cimat2mxArray: Data has zero columns");
  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      *temp++ = (int) in[r][c];
    }
  }
}

void Cmat2mxArray(double **in, mxArray *out)
{
  int r, c;
  double* temp = (double *) mxGetPr(out);
  if (temp == 0) mexErrMsgTxt("Cmat2mxArray: Pointer to data is NULL");
  int rows = mxGetM(out);
  if (rows == 0) mexErrMsgTxt("Cmat2mxArray: Data has zero rows");
  int cols = mxGetN(out);
  if (cols == 0) mexErrMsgTxt("Cmat2mxArray: Data has zero columns");
  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      *temp++ = in[r][c];
    }
  }
}

void Ccmat2mxArray(double **in_real, double **in_imag, mxArray *out)
{
  int r, c;
  double* tempR = (double *) mxGetPr(out);
  double* tempI = (double *) mxGetPi(out);
  if (tempR == 0) mexErrMsgTxt("Ccmat2mxArray: Pointer to real valued part is NULL");
  if (tempI == 0) mexErrMsgTxt("Ccmat2mxArray: Pointer to imaginary valued part is NULL");
  int rows = mxGetM(out);
  if (rows == 0) mexErrMsgTxt("Ccmat2mxArray: Data has zero rows");
  int cols = mxGetN(out);
  if (cols == 0) mexErrMsgTxt("Ccmat2mxArray: Data has zero columns");
  for (c = 0; c < cols; c++) {
    for (r = 0; r < rows; r++) {
      *tempR++ = in_real[r][c];
      *tempI++ = in_imag[r][c];
    }
  }
}

} // namespace itpp

#endif // #ifndef ITMEX_H
