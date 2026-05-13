/*!
 * \file
 * \brief Implementation of an argument parser class
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

#include <itpp/base/parser.h>
#include <fstream>
#include <sstream>


using std::cout;
using std::endl;

namespace itpp
{

Parser::Parser()
{
  VERBOSE = true;
}

Parser::Parser(const std::string &filename)
{
  VERBOSE = true;
  init(filename);
}

Parser::Parser(int argc, char *argv[])
{
  VERBOSE = true;
  init(argc, argv);
}

Parser::Parser(const std::string &filename, int argc, char *argv[])
{
  VERBOSE = true;
  init(filename, argc, argv);
}

Parser::Parser(const Array<std::string> &setup)
{
  VERBOSE = true;
  init(setup);
}

void Parser::pre_parsing(void)
{
  int i, j, n, k;
  int count = 0;
  int size = SetupStrings.size();
  bool cont_line;
  std::string Line, NewLine;
  Array<std::string> TempSetupStrings;

  // Remove lines starting with '%' or have zero length:
  for (i = 0; i < size; i++) {
    Line = SetupStrings(i);
    if ((Line.length() != 0) && (Line[0] != '%')) {
      SetupStrings(count) = Line;
      count++;
    }
  }
  SetupStrings.set_size(count, true);
  size = SetupStrings.size();

  //Check for '...' continuation as in Matlab:
  count = 0;
  NewLine = "";
  for (i = 0; i < SetupStrings.size(); i++) {
    Line = SetupStrings(i);
    n = int(Line.size());
    cont_line = false;
    for (k = 0; k < (n - 2); k++) {
      if ((Line[k] == '.') && (Line[k+1] == '.') && (Line[k+2] == '.')) {
        cont_line = true;
        break;
      }
    }
    if (cont_line) {
      for (j = 0; j < k; j++) { NewLine += Line[j]; }
    }
    else {
      NewLine += Line;
      SetupStrings(count) = NewLine;
      count++;
      NewLine = "";
    }
  }
  SetupStrings.set_size(count, true);

  //Strip unneccesary blanks, tabs, and newlines:
  size = SetupStrings.size();
  for (i = 0; i < size; i++) {
    NewLine = "";
    Line = SetupStrings(i);
    n = int(Line.length());
    j = 0; //counter in Line
    while (j < n) {
      switch (Line[j]) {
      case '\n':
        //Remove newlines
        j++;
        break;
      case ' ':
        //Remove blanks
        j++;
        break;
      case '\t':
        //Remove tabs
        j++;
        break;
      case '[':
        //Don't remove blanks between '[' and ']'
        while ((j < n) && (Line[j] != ']')) { NewLine += Line[j]; j++; }
        if (j < n) { NewLine += Line[j]; j++; }
        break;
      case '{':
        //Don't remove blanks between '{' and '}'
        while ((j < n) && (Line[j] != '}')) { NewLine += Line[j]; j++; }
        if (j < n) { NewLine += Line[j]; j++; }
        break;
      case '"':
        //Don't remove blanks between '"' and '"'
        NewLine += Line[j];
        j++; //Read in the first '"'
        while ((j < n) && (Line[j] != '"')) { NewLine += Line[j]; j++; }
        NewLine += Line[j];
        j++;
        break;
      case '\'':
        //Don't remove blanks between '\'' and '\''
        NewLine += Line[j];
        j++; //Read in the first '\''
        while ((j < n) && (Line[j] != '\'')) { NewLine += Line[j]; j++; }
        NewLine += Line[j];
        j++;
        break;
      case '%':
        //Remove trailing comments
        j = n;
        break;
      default:
        //Keep the current character:
        NewLine += Line[j];
        j++;
        break;
      }
    }
    SetupStrings(i) = NewLine;
  }

  // Split lines with several expressions (i.e. a=3;b=[1 2 3],c="Hello World") on the same line
  // (separated by comma or semicolon)
  TempSetupStrings.set_size(size, false);
  count = 0; //Counter in TempSetupStrings
  for (i = 0; i < size; i++) {

    NewLine = "";
    Line = SetupStrings(i);
    n = int(Line.length());
    j = 0;

    while (j < n) {

      switch (Line[j]) {

      case '[':
        //A vector or a matrix
        while ((j < n) && (Line[j] != ']')) { NewLine += Line[j]; j++; }
        if (Line[j] == ']') { NewLine += Line[j]; j++; }
        if (j == n) {
          if (count >= TempSetupStrings.size()) { TempSetupStrings.set_size(2*count, true); }
          TempSetupStrings(count) = NewLine;
          NewLine = "";
          count++;
        }
        break;

      case '"':
        //A string
        NewLine += Line[j];
        j++; //Read in the first '"'
        while ((j < n) && (Line[j] != '"')) { NewLine += Line[j]; j++; }
        if (Line[j] == '"') { NewLine += Line[j]; j++; }
        if (j == n) {
          if (count >= TempSetupStrings.size()) { TempSetupStrings.set_size(2*count, true); }
          TempSetupStrings(count) = NewLine;
          NewLine = "";
          count++;
        }
        break;


      case '\'':
        //A string
        NewLine += Line[j];
        j++; //Read in the first '\''
        while ((j < n) && (Line[j] != '\'')) { NewLine += Line[j]; j++; }
        if (Line[j] == '\'') { NewLine += Line[j]; j++; }
        if (j == n) {
          if (count >= TempSetupStrings.size()) { TempSetupStrings.set_size(2*count, true); }
          TempSetupStrings(count) = NewLine;
          NewLine = "";
          count++;
        }
        break;

      case ',':
        //Expression ends here
        if (count >= TempSetupStrings.size()) { TempSetupStrings.set_size(2*count, true); }
        TempSetupStrings(count) = NewLine;
        NewLine = "";
        count++;
        j++;
        break;

      case ';':
        //Expression ends here
        if (count >= TempSetupStrings.size()) { TempSetupStrings.set_size(2*count, true); }
        TempSetupStrings(count) = NewLine + ';';
        NewLine = "";
        count++;
        j++;
        break;

      default:
        //Copy the current character:
        NewLine += Line[j];
        j++;
        if (j == n) {
          if (count >= TempSetupStrings.size()) { TempSetupStrings.set_size(2*count, true); }
          TempSetupStrings(count) = NewLine;
          NewLine = "";
          count++;
        }
        break;

      }

    }

  }
  TempSetupStrings.set_size(count, true);
  SetupStrings = TempSetupStrings;
}

void Parser::init(const std::string &filename)
{
  std::string Line;
  SetupStrings.set_size(0, false);
  std::ifstream SetupFile(filename.c_str());
  it_assert(SetupFile.is_open(),
            "Parser::init(): Could not open `" + filename + "' file");

  while (getline(SetupFile, Line, '\n')) {
    SetupStrings.set_size(SetupStrings.size() + 1, true);
    SetupStrings(SetupStrings.size() - 1) = Line;
  }

  SetupFile.close();
  pre_parsing();
}

void Parser::init(int argc, char *argv[])
{
  SetupStrings.set_size(argc);
  int i;

  for (i = 0; i < argc; i++) {
    SetupStrings(i) = argv[i];
  }
  pre_parsing();
}

void Parser::init(const std::string &filename, int argc, char *argv[])
{
  std::string Line;
  int i;
  std::ifstream SetupFile(filename.c_str());
  it_assert(SetupFile.is_open(),
            "Parser::init(): Could not open `" + filename + "' file");

  //Read the command line parameters:
  SetupStrings.set_size(argc);
  for (i = 0; i < argc; i++) {
    SetupStrings(i) = argv[i];
  }

  //Read the file parameters:
  while (getline(SetupFile, Line, '\n')) {
    SetupStrings.set_size(SetupStrings.size() + 1, true);
    SetupStrings(SetupStrings.size() - 1) = Line;
  }

  SetupFile.close();
  pre_parsing();
}

void Parser::init(const Array<std::string> &setup)
{
  SetupStrings = setup;
  pre_parsing();
}

void Parser::set_silentmode(bool v)
{
  VERBOSE = !v;
}

bool Parser::exist(const std::string &name)
{
  bool error_flag, print_flag;
  std::string temp = findname(name, error_flag, print_flag);
  if (error_flag) {
    return false;
  }
  else {
    return true;
  }
}

template<>
bool Parser::get(std::string &var, const std::string &name, int num)
{
  bool error_flag, print_flag;
  std::string str = findname(name, error_flag, print_flag, num);
  if (error_flag) {
    if (VERBOSE) {
      cout << name << " = '" << var << "';" << endl;
    }
  }
  else {
    var = str;
    if (print_flag) {
      cout << name << " = '" << var << "'" << endl;
    }
    else if (VERBOSE) {
      cout << name << " = '" << var << "';" << endl;
    }
  }
  return !error_flag;
}

template<>
bool Parser::get(int &var, const std::string &name, int num)
{
  ivec out;
  bool error_flag, print_flag;
  out = ivec(findname(name, error_flag, print_flag, num));
  if (error_flag) {
    if (VERBOSE) {
      cout << name << " = " << var << ";" << endl;
    }
  }
  else {
    it_assert(out.size() == 1, "Parser::get(int): Improper variable string: "
              + name);
    var = out(0);
    if (print_flag) {
      cout << name << " = " << var << endl;
    }
    else if (VERBOSE) {
      cout << name << " = " << var << ";" << endl;
    }
  }
  return !error_flag;
}

template<>
bool Parser::get(bool &var, const std::string &name, int num)
{
  std::string ss;
  bool error_flag, print_flag;
  ss = findname(name, error_flag, print_flag, num);
  if (error_flag) {
    if (VERBOSE) {
      cout << name << " = " << var << ";" << endl;
    }
  }
  else {
    if ((ss == "true") || (ss == "1")) {
      var = true;
    }
    else if ((ss == "false") || (ss == "0")) {
      var = false;
    }
    else {
      it_error("Parser::get(bool): Improper variable string: " + name);
    }
    if (print_flag) {
      cout << name << " = " << var << endl;
    }
    else if (VERBOSE) {
      cout << name << " = " << var << ";" << endl;
    }
  }
  return !error_flag;
}

bool Parser::get_bool(const std::string &name, int num)
{
  std::string ss;
  bool out = false;
  bool error_flag, print_flag;
  ss = findname(name, error_flag, print_flag, num);
  it_assert(!error_flag, "Parser::get_bool(): Can not find variable: " + name);
  if ((ss == "true") || (ss == "1")) {
    out = true;
  }
  else if ((ss == "false") || (ss == "0")) {
    out = false;
  }
  else {
    it_error("Parser::get_bool(): Improper variable string: " + name);
  }
  if (print_flag) { cout << "Parsing bool   : " << name << " = " << out << endl; }
  return out;
}

int Parser::get_int(const std::string &name, int num)
{
  ivec out;
  bool error_flag, print_flag;
  out = ivec(findname(name, error_flag, print_flag, num));
  it_assert(!error_flag, "Parser::get_int(): Can not find variable: " + name);
  it_assert(out.size() == 1, "Parser::get_int(): Improper variable string: "
            + name);
  if (print_flag) {
    cout << "Parsing int   : " << name << " = " << out(0) << endl;
  }
  return out(0);
}

double Parser::get_double(const std::string &name, int num)
{
  double out;
  bool error_flag, print_flag;
  std::istringstream ss(findname(name, error_flag, print_flag, num));
  ss >> out;
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing double: " << name << " = " << out << endl; }
  return out;
}

std::string Parser::get_string(const std::string &name, int num)
{
  std::string out;
  bool error_flag, print_flag;
  out = findname(name, error_flag, print_flag, num);
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing string: " << name << " = " << out << endl; }
  return out;
}

vec Parser::get_vec(const std::string &name, int num)
{
  vec out;
  bool error_flag, print_flag;
  out = vec(findname(name, error_flag, print_flag, num));
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing vec   : " << name << " = " << out << endl; }
  return out;
}

ivec Parser::get_ivec(const std::string &name, int num)
{
  ivec out;
  bool error_flag, print_flag;
  out = ivec(findname(name, error_flag, print_flag, num));
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing ivec  : " << name << " = " << out << endl; }
  return out;
}

svec Parser::get_svec(const std::string &name, int num)
{
  svec out;
  bool error_flag, print_flag;
  out = svec(findname(name, error_flag, print_flag, num));
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing svec  : " << name << " = " << out << endl; }
  return out;
}

bvec Parser::get_bvec(const std::string &name, int num)
{
  bvec out;
  bool error_flag, print_flag;
  out = bvec(findname(name, error_flag, print_flag, num));
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing bvec  : " << name << " = " << out << endl; }
  return out;
}

mat Parser::get_mat(const std::string &name, int num)
{
  mat out;
  bool error_flag, print_flag;
  out = mat(findname(name, error_flag, print_flag, num));
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing mat   : " << name << " = " << out << endl; }
  return out;
}

imat Parser::get_imat(const std::string &name, int num)
{
  imat out;
  bool error_flag, print_flag;
  out = imat(findname(name, error_flag, print_flag, num));
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing imat  : " << name << " = " << out << endl; }
  return out;
}

smat Parser::get_smat(const std::string &name, int num)
{
  smat out;
  bool error_flag, print_flag;
  out = smat(findname(name, error_flag, print_flag, num));
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing smat  : " << name << " = " << out << endl; }
  return out;
}

bmat Parser::get_bmat(const std::string &name, int num)
{
  bmat out;
  bool error_flag, print_flag;
  out = bmat(findname(name, error_flag, print_flag, num));
  if (error_flag) { it_error("Parser: Can not find variable: " + name); }
  if (print_flag) { cout << "Parsing bmat  : " << name << " = " << out << endl; }
  return out;
}

std::string Parser::findname(const std::string &name, bool &error_flag, bool &print_flag, int num, bool keep_brackets)
{
  std::string Name, Out, Line, Temp;
  int n, j = 0, i = 0, index = -1;
  bool found = false, vec_mat = false;

  error_flag = false;
  print_flag = false;
  if (num < 0) { num = 0; }
  Name = "";

  // Go through all words in input string to find a match
  while (i < SetupStrings.size()) {
    Line = SetupStrings(i);
    i++;

    // Find equal sign "=", and take the left part to be the Name
    if (Line.find_first_of("=") != std::string::npos) {
      Name = Line.substr(0, Line.find_first_of("="));
    }
    else {
      Name = "";
    }

    if (Name == name) {
      if (found) {
        //cout << "Parser Warning: Duplicate Entry of variable " << name << endl;
      } else {
        found = true;
        index = i - 1;
      }
    }
  }

  // if we have a match
  if ((found) && (index + num <= SetupStrings.size())) {
    Line = SetupStrings(index + num);
    if (num != 0) {
      Temp = Line;
    }
    else {
      Temp = Line.substr(Line.find_first_of("=") + 1);
    }
  }
  else {
    error_flag = true;
    return "";
  }

  //Remove [, ],",' and ending ;. Set the print_flag:
  n = int(Temp.size());
  Out = "";
  for (i = 0; i < n; i++) {
    switch (Temp[i]) {
    case '[':
      vec_mat = true;
      if (keep_brackets) { Out += Temp[i]; }
      if (i == (n - 1)) { print_flag = true; }
      break;
    case ']':
      if (keep_brackets) { Out += Temp[i]; }
      if (i == (n - 1)) { print_flag = true; }
      break;
    case '"':
      if (i == (n - 1)) { print_flag = true; }
      break;
    case '\'':
      if (i == (n - 1)) { print_flag = true; }
      break;
    case ';':
      if (i == (n - 1)) { print_flag = false; }
      else { Out += Temp[i]; }
      break;
    default:
      Out += Temp[i];
      if (i == (n - 1)) { print_flag = true; }
      break;
    }
  }

  //Final parsing of vectors and matrices:
  if (vec_mat) {

    Temp = Out;
    Out  = "";
    n    = int(Temp.size());
    j    = 0;

    while ((Temp[j]   == ' ') || (Temp[j]   == '\t') || (Temp[j]   == '\n')) { j++; } //Remove spaces/tabs/newline in beginning
    while ((Temp[n-1] == ' ') || (Temp[n-1] == '\t') || (Temp[n-1] == '\n')) { n--; } //Remove spaces/tabs/newline at the end

    while (j < n) {

      //Read in data element:
      while ((j < n) && (Temp[j] != ' ') && (Temp[j] != '\t') && (Temp[j] != '\n') && (Temp[j] != ',') && (Temp[j] != ';')) {
        Out += Temp[j];
        j++;
      }

      //Read separator:
      if (j < (n - 1)) {
        //Remove spaces before separator
        while ((j < n) && ((Temp[j] == ' ') || (Temp[j] == '\t') || (Temp[j] == '\n'))) { j++; }
        //Insert Separator:
        if (j < n) {
          if (Temp[j] == ';') { Out += ';'; j++; }
          else if (Temp[j] == ',') { Out += ' '; j++; }
          else if (Temp[j] == '+') { Out += ' '; Out += Temp[j]; j++; }
          else if (Temp[j] == '-') { Out += ' '; Out += Temp[j]; j++; }
          else { Out += ' '; }
        }
        //Remove spaces after separator:
        while ((j < n) && ((Temp[j] == ' ') || (Temp[j] == '\t') || (Temp[j] == '\n'))) { j++; }
      }

    }

  }
  if (!VERBOSE) print_flag = false;
  return Out;
}

} // namespace itpp
