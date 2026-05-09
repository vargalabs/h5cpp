/*!
 * \file
 * \brief Implementation of PNM graphics format I/O function
 * \author
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

#include <itpp/srccode/pnm.h>
#include <itpp/base/itassert.h>
#include <fstream>

//! \cond

using std::istream;
using std::ostream;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::ios;
using std::ios_base;
using std::streampos;


namespace itpp
{


// Suppress the additional white characters and return the comments
static void pnm_read_comments(istream & i, string & comments);

// Write comment in the image file
static void pnm_write_comments(ostream & o, const string & comments);

// Read/Write the header for the pnm file format
static bool pnm_read_header(ifstream & file, char & pnm_type,
                            int & width, int & height, int & max_val,
                            string & comments, char pnm_type_required = '0');

static bool pnm_write_header(ofstream & file, char type,
                             int width, int height, int max_val,
                             const string & comments);


//--------------------------------------------------------------
// General PNM functions
//--------------------------------------------------------------
char pnm_type(const string & filename)
{
  ifstream file;
  char pnm_type;

  file.open(filename.c_str(), ifstream::in | ifstream::binary);

  string comments;
  int width, height, max_val;
  pnm_read_header(file, pnm_type, width, height, max_val, comments);

  return pnm_type;
}


//--------------------------------------------------------------
bool pnm_info(const string & filename, char & pnm_type,
              int & width, int & height, int & max_val,
              string & comments)
{
  ifstream file;

  file.open(filename.c_str(), ifstream::in | ifstream::binary);

  pnm_read_header(file, pnm_type, width, height, max_val, comments);

  return true;
}


//--------------------------------------------------------------
// PGM related functions (gray images)
//--------------------------------------------------------------

bool pgm_read(const string & filename,
              imat & m, string & comments)
{
  ifstream file;
  int width, height, max_val, i, j;
  comments = "";

  file.open(filename.c_str(), ifstream::in | ifstream::binary);

  // The format code is 'P5' for pgm files
  char pnm_type;
  if (!pnm_read_header(file, pnm_type, width, height, max_val, comments, '5'))
    return false;

  // Format the returned matrix
  m.set_size(height, width, false);

  // Retrieve the integer value from the file
  for (i = 0 ; i < height; i++)
    for (j = 0; j < width; j++)
      m(i, j) = file.get();

  return true;
}


//--------------------------------------------------------------
// Simplified version of read_pgm
imat pgm_read(const string & filename)
{
  imat I;
  string comments;
  if (!pgm_read(filename, I, comments)) {
    it_warning("pgm_read (PGM file->imat) failed ");
  }
  return I;
}


//--------------------------------------------------------------
bool pgm_read(const string & filename, imat &m,
              int r1, int r2, int c1, int c2)
{
  ifstream file;
  int width, height, max_val, i, j;

  // This is a dummy variable.
  // Its purpose is the call of function pnm_read_header.
  string comments;

  file.open(filename.c_str(), ifstream::in | ifstream::binary);

  char pnm_type;
  if (!pnm_read_header(file, pnm_type, width, height, max_val, comments, '5'))
    return false;

  // Inversion of the column/row numbers may be required
  if (r1 > r2) {
    int rtmp = r2;
    r2 = r1;
    r1 = rtmp;
  }

  if (c1 > c2) {
    int ctmp = c2;
    c2 = c1;
    c1 = ctmp;
  }

  it_error_if((r1 < 0) || (c1 < 0),
              "Bad parameter value: row and column number must be >=0");
  it_error_if((r2 >= height) || (c1 >= width), "Bad parameter value: "
              "row or column number exceeds the image heigth");

  m.set_size(r2 - r1 + 1, c2 - c1 + 1, false);
  file.seekg(r1 * width + c1, ios::cur);

  for (i = 0 ; i < m.rows() ; i++) {
    for (j = 0 ; j < m.cols() ; j++)
      m(i, j) = file.get();
    file.seekg(width - (c2 - c1 + 1), ios::cur);
  }

  return true;
}


//--------------------------------------------------------------
bool pgm_write(const string & filename,
               const imat &m, const string & comments)
{

  ofstream file;
  int i, j;

  file.open(filename.c_str(), ofstream::out | ofstream::binary);

  if (!pnm_write_header(file, '5', m.cols(), m.rows(), 255, comments))
    return false;

  for (i = 0; i < m.rows(); i++)
    for (j = 0; j < m.cols(); j++)
      file.put(static_cast<char>(m(i, j)));

  if (!file)
    return false;

  return true;
}


//--------------------------------------------------------------
// PPM related functions (color images)
//--------------------------------------------------------------

bool ppm_read(const string & filename,
              imat &r, imat &g, imat &b,
              string & comments)
{
  ifstream file;
  int width, height, max_val, i, j;

  file.open(filename.c_str(), ifstream::in | ifstream::binary);

  char pnm_type;
  if (!pnm_read_header(file, pnm_type, width, height, max_val, comments, '6'))
    return false;

  r.set_size(height, width, false);
  g.set_size(height, width, false);
  b.set_size(height, width, false);
  for (i = 0; i < height; i++)
    for (j = 0; j < width; j++) {
      r(i, j) = file.get();
      g(i, j) = file.get();
      b(i, j) = file.get();
    }

  return true;
}


//--------------------------------------------------------------
// Same function but suppress the comments
bool ppm_read(const string & filename,
              imat &r, imat &g, imat &b)
{
  string comments; // This is a dummy variable

  return ppm_read(filename, r, g, b, comments);
}

//--------------------------------------------------------------
bool ppm_read(const string & filename,
              imat &r, imat &g, imat &b,
              int r1, int r2, int c1, int c2)
{
  ifstream file;
  int width, height, max_val, i, j;

  // This is a dummy variable. Its purpose is the call of function pnm_read_header.
  string comments;

  file.open(filename.c_str(), ifstream::in | ifstream::binary);

  char pnm_type;
  if (!pnm_read_header(file, pnm_type, width, height, max_val, comments, '6'))
    return false;

  // Inversion of the column/row numbers may be required
  if (r1 > r2) {
    // Funny way to do it... (without using any temporary variable)
    r1 += r2;
    r2 = r1 - r2;
    r1 -= r2;
  }

  if (c1 > c2) {
    // Conventionnal way to do it
    int ctmp = c2;
    c2 = c1;
    c1 = ctmp;
  }

  it_error_if((r1 < 0) || (c1 < 0),
              "Bad parameter value: row and column number must be >=0");
  it_error_if((r2 >= height) || (c1 >= width), "Bad parameter value: "
              "row or column number exceeds the image heigth");

  r.set_size(r2 - r1 + 1, c2 - c1 + 1, false);
  g.set_size(r2 - r1 + 1, c2 - c1 + 1, false);
  b.set_size(r2 - r1 + 1, c2 - c1 + 1, false);
  file.seekg(3 *(r1 * width + c1), ios::cur);

  for (i = 0; i < r.rows(); i++) {
    for (j = 0; j < r.cols(); j++) {
      r(i, j) = file.get();
      g(i, j) = file.get();
      b(i, j) = file.get();
    }
    file.seekg(3 * (width - (c2 - c1 + 1)), ios::cur);
  }

  return true;
}


//--------------------------------------------------------------
bool ppm_write(const string & filename,
               const imat &r, const imat &g, const imat &b,
               const string & comments,
               int max_val)
{
  ofstream file;
  int i, j;

  it_assert_debug(r.cols() == g.cols() && g.cols() == b.cols() &&
                  r.rows() == g.rows() && g.rows() == b.rows(),
                  "Matrices r, g and b must have the same size in ppm_write()");

  file.open(filename.c_str(), ofstream::out | ofstream::binary);

  if (max_val < 0 || max_val > 65535) {
    it_warning("Proposed maximal value is incorrect");
    return false;
  }

  if (!pnm_write_header(file, '6', r.cols(), r.rows(), max_val, comments))
    return false;

  for (i = 0; i < r.rows(); i++)
    for (j = 0; j < r.cols(); j++) {
      file.put(static_cast<char>(r(i, j)));
      file.put(static_cast<char>(g(i, j)));
      file.put(static_cast<char>(b(i, j)));
    }

  if (!file)
    return false;

  return true;
}


//--------------------------------------------------------------
imat img_double2int(const mat & m,
                    int max_val,
                    double double_min,
                    double double_max)
{
  int i, j;
  imat M(m.rows(), m.cols());

  for (i = 0 ; i < m.rows() ; i++)
    for (j = 0 ; j < m.cols() ; j++)
      if (m(i, j) <= double_min)
        M(i, j) = 0;

      else if (m(i, j) >= double_max)
        M(i, j) = max_val;

      else
        M(i, j) = (int)(max_val * (m(i, j) - double_min)
                        / (double_max - double_min) + 0.5);

  return M;
}

//--------------------------------------------------------------
mat img_int2double(const imat & m,
                   int max_val,
                   double double_min,
                   double double_max)
{
  int i, j;
  mat M(m.rows(), m.cols());

  for (i = 0 ; i < m.rows() ; i++)
    for (j = 0 ; j < m.cols() ; j++)
      if (m(i, j) <= 0)
        M(i, j) = double_min;

      else if (m(i, j) >= max_val)
        M(i, j) = double_max;

      else
        // This rounding works well when m(i,j) is positive
        M(i, j) = double_min + (double_max - double_min)
                  * m(i, j) / (double) max_val;

  return M;
}


//--------------------------------------------------------------
// Static functions: Used in this file only
//--------------------------------------------------------------

//--------------------------------------------------------------
static void pnm_read_comments(istream & i, string & comments)
{
  while (isspace(i.peek())) {
    while (isspace(i.peek()))
      i.get();

    if (i.peek() == '#')
      while (i.peek() != '\r' && i.peek() != '\n')
        comments += static_cast<char>(i.get());
  }
}


//--------------------------------------------------------------
static void pnm_write_comments(ostream & o, const string & comments)
{
  istringstream comments_stream(comments);
  char comment_line[ 256 ];

  // Put header and comment
  while (!comments_stream.eof()) {
    o << "#";
    comments_stream.get(comment_line, 256);
    o << comment_line << endl;
  }
}


//--------------------------------------------------------------
// Read the header of a pnm file
static bool pnm_read_header(ifstream & file, char & pnm_type,
                            int & width, int & height, int & max_val,
                            string & comments, char pnm_type_required)
{
  bool return_code = true;

  if (file.get() != 'P')
    return_code = false;
  it_error_if(!return_code, "Invalid format file: code of file format has "
              "not been found");

  // Read the type of the pnm file
  file.get(pnm_type);
  it_error_if((pnm_type < '1') || (pnm_type > '6'),
              "Bad file code P" << pnm_type);

  // If a type has been specified
  if (pnm_type_required != '0')
    if (pnm_type_required != pnm_type) {
      string err_msg("Found file code P");
      err_msg += pnm_type + " instead of P" + pnm_type_required;
      it_error(err_msg);
    }

  // Retrieve the image format and the comments
  pnm_read_comments(file, comments);
  file >> width;
  pnm_read_comments(file, comments);
  file >> height;
  pnm_read_comments(file, comments);

  it_error_if((height < 0) || (width < 0), "Bad image size");

  // Maximal values is not present in PBM files
  if (pnm_type == '2' || pnm_type == '3' || pnm_type == '5' || pnm_type == '6')
    file >> max_val;

  file.get(); // Eat the last whitespace

  // According to the pnm specification, the maximal value should not
  // be greater than 65536 and lower than 0
  it_error_if((max_val >= 65536) || (max_val < 0),
              "Invalid maximum number in pnm header");

  // For type P5 and P6, the value have to be lower than 255
  it_error_if((pnm_type == '5' || pnm_type == '6') && (max_val > 255),
              "Invalid maximum number in pnm header");

  return file.good();
}


//--------------------------------------------------------------
static bool pnm_write_header(ofstream &file, char pnm_type,
                             int width, int height, int max_val,
                             const string & comments)
{
  file << 'P' << pnm_type << endl;
  pnm_write_comments(file, comments);
  file << width << ' ' << height << endl;

  // Maximal values is not present in PBM files
  if (pnm_type == '2' || pnm_type == '3' || pnm_type == '5' || pnm_type == '6')
    file << max_val << endl;

  return file.good();
}

} // namespace itpp

//! \endcond
