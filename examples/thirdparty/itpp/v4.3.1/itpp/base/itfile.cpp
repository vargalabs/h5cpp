/*!
 * \file
 * \brief Implementation of classes for the IT++ file format
 * \author Tony Ottosson, Tobias Ringstrom and Adam Piatyszek
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


#include <itpp/base/itfile.h>


namespace itpp
{

char it_file_base::file_magic[4] = { 'I', 'T', '+', '+' };
char it_file_base::file_version = 3;

// ----------------------------------------------------------------------
// it_ifile class
// ----------------------------------------------------------------------

it_ifile::it_ifile() {}

it_ifile::it_ifile(const std::string &name)
{
  open(name);
}

void it_ifile::open(const std::string &name)
{
  it_assert(exist(name), "it_ifile::open(): File does not exist");
  s.open_readonly(name, bfstream_base::l_endian);
  if (!read_check_file_header()) {
    s.close();
    it_error("it_ifile::open(): Corrupt file (not an it_file)");
  }
}

void it_ifile::close()
{
  s.close();
}

bool it_ifile::seek(const std::string &name)
{
  data_header h;
  std::streampos p;

  s.clear();
  s.seekg(sizeof(file_header));

  while (true) {
    p = s.tellg();
    read_data_header(h);
    if (s.eof()) {
      s.clear();
      return false;
    }
    if (h.type != "" && h.name == name) {
      s.seekg(p);
      break;
    }
    s.seekg(p + static_cast<std::streamoff>(h.block_bytes));
  }

  return true;
}

bool it_ifile::seek(int n)
{
  data_header h;
  std::streampos p;

  s.clear();
  s.seekg(sizeof(file_header));
  for (int i = 0; i <= n; i++) {
    p = s.tellg();
    read_data_header(h);
    if (s.eof()) {
      s.clear();
      return false;
    }
    if (h.type == "")
      i--;
    s.seekg((i == n) ? p : p + static_cast<std::streamoff>(h.block_bytes));
  }
  return true;
}

void it_ifile::info(std::string &name, std::string &type,
                    std::string &desc, uint64_t &bytes)
{
  data_header h;
  std::streampos p;

  p = s.tellg();
  read_data_header(h);
  s.seekg(p);
  name = h.name;
  type = h.type;
  desc = h.desc;
  bytes = h.data_bytes;
}

bool it_ifile::read_check_file_header()
{
  file_header h;
  s.read(reinterpret_cast<char *>(&h), sizeof(h));
  return (memcmp(h.magic, file_magic, 4) == 0
          && (h.version == file_version));
}

void it_ifile::read_data_header(data_header &h)
{
  s.clear();
  s >> h.hdr_bytes;
  s >> h.data_bytes;
  s >> h.block_bytes;
  s >> h.name;
  s >> h.type;
  s >> h.desc;
}

void it_ifile::low_level_read(char &x)
{
  s >> x;
}

void it_ifile::low_level_read(uint64_t &x)
{
  s >> x;
}

void it_ifile::low_level_read(bool &x)
{
  char tmp;
  s >> tmp;
  x = (tmp == 0) ? false : true;
}


void it_ifile::low_level_read(bin &x)
{
  char tmp;
  s >> tmp;
  x = tmp;
}

void it_ifile::low_level_read(short &x)
{
  int16_t tmp;
  s >> tmp;
  x = tmp;
}

void it_ifile::low_level_read(int &x)
{
  int32_t tmp;
  s >> tmp;
  x = tmp;
}

void it_ifile::low_level_read(float &x)
{
  s >> x;
}

void it_ifile::low_level_read(double &x)
{
  s >> x;
}

void it_ifile::low_level_read(std::complex<float> &x)
{
  float x_real, x_imag;
  s >> x_real;
  s >> x_imag;
  x = std::complex<float>(x_real, x_imag);
}

void it_ifile::low_level_read(std::complex<double> &x)
{
  double x_real, x_imag;
  s >> x_real;
  s >> x_imag;
  x = std::complex<double>(x_real, x_imag);
}

void it_ifile::low_level_read(bvec &v)
{
  uint64_t size;
  char tmp;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> tmp;
    v(i) = tmp;
  }
}

void it_ifile::low_level_read(svec &v)
{
  uint64_t size;
  int16_t val;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val;
    v(i) = val;
  }
}

void it_ifile::low_level_read(ivec &v)
{
  uint64_t size;
  int32_t val;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val;
    v(i) = val;
  }
}

void it_ifile::low_level_read_lo(vec &v)
{
  uint64_t size;
  float val;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val;
    v(i) = static_cast<double>(val);
  }
}

void it_ifile::low_level_read_hi(vec &v)
{
  uint64_t size;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i)
    s >> v(i);
}

void it_ifile::low_level_read_lo(cvec &v)
{
  uint64_t size;
  float val_real, val_imag;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<double>(val_real, val_imag);
  }
}

void it_ifile::low_level_read_hi(cvec &v)
{
  uint64_t size;
  double val_real, val_imag;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<double>(val_real, val_imag);
  }
}

void it_ifile::low_level_read(std::string &str)
{
  uint64_t size;
  s >> size;
  std::string::size_type size2 = static_cast<std::string::size_type>(size);
  str.resize(size2);
  for (std::string::size_type i = 0; i < size2; ++i)
    s >> str[i];
}

void it_ifile::low_level_read(bmat &m)
{
  uint64_t i, j;
  char tmp;
  s >> i >> j;
  m.set_size(static_cast<int>(i), static_cast<int>(j), false);
  for (int j = 0; j < m.cols(); ++j) {
    for (int i = 0; i < m.rows(); ++i) {
      s >> tmp;
      m(i, j) = tmp;
    }
  }
}

void it_ifile::low_level_read(smat &m)
{
  uint64_t i, j;
  int16_t val;
  s >> i >> j;
  m.set_size(static_cast<int>(i), static_cast<int>(j), false);
  for (int j = 0; j < m.cols(); ++j)
    for (int i = 0; i < m.rows(); ++i) {
      s >> val;
      m(i, j) = val;
    }
}

void it_ifile::low_level_read(imat &m)
{
  uint64_t i, j;
  int32_t val;
  s >> i >> j;
  m.set_size(static_cast<int>(i), static_cast<int>(j), false);
  for (int j = 0; j < m.cols(); ++j)
    for (int i = 0; i < m.rows(); ++i) {
      s >> val;
      m(i, j) = val;
    }
}

void it_ifile::low_level_read_lo(mat &m)
{
  uint64_t i, j;
  float val;
  s >> i >> j;
  m.set_size(static_cast<int>(i), static_cast<int>(j), false);
  for (int j = 0; j < m.cols(); ++j)
    for (int i = 0; i < m.rows(); ++i) {
      s >> val;
      m(i, j) = static_cast<double>(val);
    }
}

void it_ifile::low_level_read_hi(mat &m)
{
  uint64_t i, j;
  s >> i >> j;
  m.set_size(static_cast<int>(i), static_cast<int>(j), false);
  for (int j = 0; j < m.cols(); ++j)
    for (int i = 0; i < m.rows(); ++i)
      s >> m(i, j);
}

void it_ifile::low_level_read_lo(cmat &m)
{
  uint64_t i, j;
  float val_real, val_imag;
  s >> i >> j;
  m.set_size(static_cast<int>(i), static_cast<int>(j), false);
  for (int j = 0; j < m.cols(); ++j)
    for (int i = 0; i < m.rows(); ++i) {
      s >> val_real;
      s >> val_imag;
      m(i, j) = std::complex<double>(val_real, val_imag);
    }
}

void it_ifile::low_level_read_hi(cmat &m)
{
  uint64_t i, j;
  double val_real, val_imag;
  s >> i >> j;
  m.set_size(static_cast<int>(i), static_cast<int>(j), false);
  for (int j = 0; j < m.cols(); ++j)
    for (int i = 0; i < m.rows(); ++i) {
      s >> val_real;
      s >> val_imag;
      m(i, j) = std::complex<double>(val_real, val_imag);
    }
}

void it_ifile::low_level_read(Array<bin> &v)
{
  uint64_t size;
  char tmp;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> tmp;
    v(i) = tmp;
  }
}

void it_ifile::low_level_read(Array<short> &v)
{
  uint64_t size;
  int16_t val;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val;
    v(i) = val;
  }
}

void it_ifile::low_level_read(Array<int> &v)
{
  uint64_t size;
  int32_t val;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val;
    v(i) = val;
  }
}

void it_ifile::low_level_read(Array<float> &v)
{
  uint64_t size;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i)
    s >> v(i);
}

void it_ifile::low_level_read_lo(Array<double> &v)
{
  uint64_t size;
  float val;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val;
    v(i) = static_cast<double>(val);
  }
}

void it_ifile::low_level_read_hi(Array<double> &v)
{
  uint64_t size;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i)
    s >> v(i);
}

void it_ifile::low_level_read(Array<std::complex<float> > &v)
{
  uint64_t size;
  float val_real, val_imag;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<float>(val_real, val_imag);
  }
}

void it_ifile::low_level_read_lo(Array<std::complex<double> > &v)
{
  uint64_t size;
  float val_real, val_imag;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<double>(val_real, val_imag);
  }
}

void it_ifile::low_level_read_hi(Array<std::complex<double> > &v)
{
  uint64_t size;
  double val_real, val_imag;
  s >> size;
  v.set_size(static_cast<int>(size), false);
  for (int i = 0; i < v.size(); ++i) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<double>(val_real, val_imag);
  }
}


// ----------------------------------------------------------------------
// it_file class
// ----------------------------------------------------------------------

it_file::it_file(): low_prec(false), _strings(new Strings_Holder) {}

it_file::it_file(const std::string &name, bool trunc):
    low_prec(false), _strings(new Strings_Holder)
{
  open(name, trunc);
}

void it_file::open(const std::string &name, bool trunc)
{
  if (!exist(name))
    trunc = true;

  s.open(name, trunc, bfstream_base::l_endian);
  it_assert(s.is_open(), "it_file::open(): Could not open file for writing");

  if (trunc)
    write_file_header();
  else if (!read_check_file_header()) {
    s.close();
    it_error("it_file::open(): Corrupt file (not an it_file)");
  }

  fname() = name;
}

void it_file::close()
{
  s.close();
}

void it_file::flush()
{
  s.flush();
}

void it_file::write_file_header()
{
  s.write(file_magic, 4);
  s.put(file_version);
}

void it_file::write_data_header(const std::string &type, uint64_t size)
{
  it_error_if(next_name() == "", "it_file::write_data_header(): Can not "
              "write without a name");
  write_data_header(type, next_name(), size, next_desc());
  next_name() = "";
  next_desc() = "";
}

void it_file::write_data_header(const std::string &type,
                                const std::string &name, uint64_t size,
                                const std::string &desc)
{
  data_header h1, h2;

  // Prepare a new data header
  h1.hdr_bytes = 3 * sizeof(uint64_t) + type.size() + 1 + name.size() + 1
                 + desc.size() + 1;
  h1.data_bytes = size;
  h1.block_bytes = h1.hdr_bytes + h1.data_bytes;
  h1.name = name;
  h1.type = type;
  h1.desc = desc;

  // If variable exists, remove it first
  if (exists(name))
    remove();

  // Try to find an empty space
  s.clear();
  s.seekg(sizeof(file_header)); // skip file header
  while (true) {
    // save the current position
    std::streampos p = s.tellp();
    // read block at the current position
    read_data_header(h2);
    // if empty file, stop the search and set write pointer to the end of
    // file
    if (s.eof()) {
      s.clear();
      s.seekp(0, std::ios::end);
      break;
    }
    // save the size of the current read block
    std::streamoff skip = static_cast<std::streamoff>(h2.block_bytes);
    // check if we have enough empty space from previously deleted data
    if ((h2.type == "") && (h2.block_bytes >= h1.block_bytes)) {
      h1.block_bytes = h2.block_bytes;
      s.seekp(p);
      break;
    }
    // if not, maybe we can squeeze the current block to find space
    else if ((h2.block_bytes - h2.hdr_bytes - h2.data_bytes)
             >= h1.block_bytes) {
      h1.block_bytes = h2.block_bytes - h2.hdr_bytes - h2.data_bytes;
      h2.block_bytes = h2.hdr_bytes + h2.data_bytes;
      s.seekp(p);
      // rewrite squeezed data block
      write_data_header_here(h2);
      s.seekp(p + static_cast<std::streamoff>(h2.block_bytes));
      break;
    }
    // otherwise, skip the current block and try again
    s.seekg(p + skip);
  } // while(true)

  write_data_header_here(h1);
}

void it_file::write_data_header_here(const data_header &h)
{
  s << h.hdr_bytes << h.data_bytes << h.block_bytes
  << h.name << h.type << h.desc;
}

void it_file::remove(const std::string &name)
{
  seek(name);
  remove();
}

void it_file::remove()
{
  data_header h;
  std::streampos p;

  p = s.tellp();
  read_data_header(h);
  h.type = "";
  h.name = "";
  h.desc = "";
  h.hdr_bytes = 3 * sizeof(uint64_t) + 1 + 1 + 1;
  h.data_bytes = 0;
  s.seekp(p);
  write_data_header_here(h);
  s.seekp(p + static_cast<std::streamoff>(h.block_bytes));
}

bool it_file::exists(const std::string &name)
{
  return seek(name);
}

void it_file::pack()
{
  it_assert(s.is_open(), "it_file::pack(): File has to be open");

  // check total file size
  s.seekg(0, std::ios::end);
  std::streampos p = s.tellg();
  s.seekg(0, std::ios::beg);
  s.clear();

  // allocate buffer of size equal to file size
  char* buffer = new char[int(p)];
  char* b_ptr = buffer;

  // copy file header and start counting the size of compacted file
  uint64_t size;
  for (size = 0; size < sizeof(file_header); ++size)
    s.get(*b_ptr++);

  // remove empty space between data blocks
  data_header h;
  while (true) {
    p = s.tellg();
    read_data_header(h);
    if (s.eof()) {
      s.clear();
      break;
    }
    if (h.type != "") {
      s.seekg(p);
      for (uint64_t i = 0; i < h.hdr_bytes + h.data_bytes; ++i)
        s.get(*b_ptr++);
      size += h.hdr_bytes + h.data_bytes;
    }
    s.seekg(p + static_cast<std::streamoff>(h.block_bytes));
  }

  // close and reopen file truncating it
  s.close();
  s.open(fname(), true, bfstream_base::l_endian);
  // write compacted data to the reopend empty file
  for (uint64_t i = 0; i < size; ++i)
    s.put(buffer[i]);

  // free buffer memory
  delete buffer;

  // go back to the first data block (skiping file header)
  s.seekg(sizeof(file_header));

  // update block_bytes in headers of compacted data blocks
  while (true) {
    p = s.tellg();
    read_data_header(h);
    if (s.eof()) {
      s.clear();
      break;
    }
    if (h.hdr_bytes + h.data_bytes < h.block_bytes) {
      h.block_bytes = h.hdr_bytes + h.data_bytes;
      s.seekp(p);
      write_data_header_here(h);
    }
    s.seekg(p + static_cast<std::streamoff>(h.block_bytes));
  }
}

void it_file::low_level_write(char x)
{
  s << x;
}

void it_file::low_level_write(uint64_t x)
{
  s << x;
}

void it_file::low_level_write(bool x)
{
  s << static_cast<char>(x);
}

void it_file::low_level_write(bin x)
{
  s << x.value();
}

void it_file::low_level_write(short x)
{
  s << static_cast<int16_t>(x);
}

void it_file::low_level_write(int x)
{
  s << static_cast<int32_t>(x);
}

void it_file::low_level_write(float x)
{
  s << x;
}

void it_file::low_level_write(double x)
{
  s << x;
}

void it_file::low_level_write(const std::complex<float> &x)
{
  s << x.real();
  s << x.imag();
}

void it_file::low_level_write(const std::complex<double> &x)
{
  s << x.real();
  s << x.imag();
}

void it_file::low_level_write(const bvec &v)
{
  s << static_cast<uint64_t>(v.size());
  for (int i = 0; i < v.size(); ++i)
    s << v(i).value();
}

void it_file::low_level_write(const svec &v)
{
  s << static_cast<uint64_t>(v.size());
  for (int i = 0; i < v.size(); ++i)
    s << static_cast<int16_t>(v(i));
}

void it_file::low_level_write(const ivec &v)
{
  s << static_cast<uint64_t>(v.size());
  for (int i = 0; i < v.size(); ++i)
    s << static_cast<int32_t>(v(i));
}

void it_file::low_level_write(const vec &v)
{
  s << static_cast<uint64_t>(v.size());
  if (get_low_precision()) {
    for (int i = 0; i < v.size(); ++i)
      s << static_cast<float>(v(i));
  }
  else {
    for (int i = 0; i < v.size(); ++i)
      s << v(i);
  }
}

void it_file::low_level_write(const cvec &v)
{
  s << static_cast<uint64_t>(v.size());
  if (get_low_precision()) {
    for (int i = 0; i < v.size(); ++i) {
      s << static_cast<float>(v(i).real());
      s << static_cast<float>(v(i).imag());
    }
  }
  else {
    for (int i = 0; i < v.size(); ++i) {
      s << v(i).real();
      s << v(i).imag();
    }
  }
}

void it_file::low_level_write(const std::string &str)
{
  s << static_cast<uint64_t>(str.size());
  for (std::string::size_type i = 0; i < str.size(); ++i)
    s << str[i];
}

void it_file::low_level_write(const bmat &m)
{
  s << static_cast<uint64_t>(m.rows())
  << static_cast<uint64_t>(m.cols());
  for (int j = 0; j < m.cols(); ++j)
    for (int i = 0; i < m.rows(); ++i)
      s << m(i, j).value();
}

void it_file::low_level_write(const smat &m)
{
  s << static_cast<uint64_t>(m.rows())
  << static_cast<uint64_t>(m.cols());
  for (int j = 0; j < m.cols(); ++j)
    for (int i = 0; i < m.rows(); ++i)
      s << static_cast<int16_t>(m(i, j));
}

void it_file::low_level_write(const imat &m)
{
  s << static_cast<uint64_t>(m.rows())
  << static_cast<uint64_t>(m.cols());
  for (int j = 0; j < m.cols(); ++j)
    for (int i = 0; i < m.rows(); ++i)
      s << static_cast<int32_t>(m(i, j));
}

void it_file::low_level_write(const mat &m)
{
  s << static_cast<uint64_t>(m.rows())
  << static_cast<uint64_t>(m.cols());
  if (get_low_precision()) {
    for (int j = 0; j < m.cols(); ++j)
      for (int i = 0; i < m.rows(); ++i)
        s << static_cast<float>(m(i, j));
  }
  else {
    for (int j = 0; j < m.cols(); ++j)
      for (int i = 0; i < m.rows(); ++i)
        s << m(i, j);
  }
}

void it_file::low_level_write(const cmat &m)
{
  s << static_cast<uint64_t>(m.rows())
  << static_cast<uint64_t>(m.cols());
  if (get_low_precision()) {
    for (int j = 0; j < m.cols(); ++j)
      for (int i = 0; i < m.rows(); ++i) {
        s << static_cast<float>(m(i, j).real());
        s << static_cast<float>(m(i, j).imag());
      }
  }
  else {
    for (int j = 0; j < m.cols(); ++j)
      for (int i = 0; i < m.rows(); ++i) {
        s << m(i, j).real();
        s << m(i, j).imag();
      }
  }
}

void it_file::low_level_write(const Array<bin> &v)
{
  s << static_cast<uint64_t>(v.size());
  for (int i = 0; i < v.size(); ++i)
    s << v(i).value();
}

void it_file::low_level_write(const Array<short> &v)
{
  s << static_cast<uint64_t>(v.size());
  for (int i = 0; i < v.size(); ++i)
    s << static_cast<int16_t>(v(i));
}

void it_file::low_level_write(const Array<int> &v)
{
  s << static_cast<uint64_t>(v.size());
  for (int i = 0; i < v.size(); ++i)
    s << static_cast<int32_t>(v(i));
}

void it_file::low_level_write(const Array<float> &v)
{
  s << static_cast<uint64_t>(v.size());
  for (int i = 0; i < v.size(); ++i)
    s << v(i);
}

void it_file::low_level_write(const Array<double> &v)
{
  s << static_cast<uint64_t>(v.size());
  if (get_low_precision()) {
    for (int i = 0; i < v.size(); ++i)
      s << static_cast<float>(v(i));
  }
  else {
    for (int i = 0; i < v.size(); ++i)
      s << static_cast<double>(v(i));
  }
}

void it_file::low_level_write(const Array<std::complex<float> > &v)
{
  s << static_cast<uint64_t>(v.size());
  for (int i = 0; i < v.size(); ++i) {
    s << v(i).real();
    s << v(i).imag();
  }
}

void it_file::low_level_write(const Array<std::complex<double> > &v)
{
  s << static_cast<uint64_t>(v.size());
  if (get_low_precision()) {
    for (int i = 0; i < v.size(); ++i) {
      s << static_cast<float>(v(i).real());
      s << static_cast<float>(v(i).imag());
    }
  }
  else {
    for (int i = 0; i < v.size(); ++i) {
      s << v(i).real();
      s << v(i).imag();
    }
  }
}


it_ifile &operator>>(it_ifile &f, char &x)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "int8", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(x);
  return f;
}

it_ifile &operator>>(it_ifile &f, bool &x)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "bool", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(x);
  return f;
}

it_ifile &operator>>(it_ifile &f, bin &x)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "bin", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(x);
  return f;
}

it_ifile &operator>>(it_ifile &f, short &x)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "int16", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(x);
  return f;
}

it_ifile &operator>>(it_ifile &f, int &x)
{
  it_file::data_header h;
  f.read_data_header(h);
  if (h.type == "int32")
    f.low_level_read(x);
  else if (h.type == "int16") {
    short x16;
    f.low_level_read(x16);
    x = static_cast<int>(x16);
  }
  else
    it_error("it_ifile::operator>>(): Wrong type");

  return f;
}

it_ifile &operator>>(it_ifile &f, float &x)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "float32", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(x);
  return f;
}

it_ifile &operator>>(it_ifile &f, double &x)
{
  it_file::data_header h;
  f.read_data_header(h);
  if (h.type == "float64")
    f.low_level_read(x);
  else if (h.type == "float32") {
    float f32;
    f.low_level_read(f32);
    x = static_cast<double>(f32);
  }
  else
    it_error("it_ifile::operator>>(): Wrong type");

  return f;
}

it_ifile &operator>>(it_ifile &f, std::complex<float> &x)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "cfloat32",
            "it_ifile::operator>>(): Wrong type");
  f.low_level_read(x);
  return f;
}

it_ifile &operator>>(it_ifile &f, std::complex<double> &x)
{
  it_file::data_header h;
  f.read_data_header(h);
  if (h.type == "cfloat64")
    f.low_level_read(x);
  else if (h.type == "cfloat32") {
    std::complex<float> f32_c;
    f.low_level_read(f32_c);
    x = static_cast<std::complex<double> >(f32_c);
  }
  else
    it_error("it_ifile::operator>>(): Wrong type");

  return f;
}

it_ifile &operator>>(it_ifile &f, bvec &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "bvec", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(v);
  return f;
}

it_ifile &operator>>(it_ifile &f, svec &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "svec", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(v);
  return f;
}

it_ifile &operator>>(it_ifile &f, ivec &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "ivec", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(v);
  return f;
}

it_ifile &operator>>(it_ifile &f, vec &v)
{
  it_ifile::data_header h;

  f.read_data_header(h);
  if (h.type == "fvec")
    f.low_level_read_lo(v);
  else if (h.type == "dvec")
    f.low_level_read_hi(v);
  else
    it_error("it_ifile::operator>>(): Wrong type");

  return f;
}

it_ifile &operator>>(it_ifile &f, cvec &v)
{
  it_file::data_header h;

  f.read_data_header(h);
  if (h.type == "fcvec")
    f.low_level_read_lo(v);
  else if (h.type == "dcvec")
    f.low_level_read_hi(v);
  else
    it_error("it_ifile::operator>>(): Wrong type");

  return f;
}

it_ifile &operator>>(it_ifile &f, std::string &str)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "string", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(str);
  return f;
}

it_ifile &operator>>(it_ifile &f, bmat &m)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "bmat", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(m);
  return f;
}

it_ifile &operator>>(it_ifile &f, smat &m)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "smat", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(m);
  return f;
}

it_ifile &operator>>(it_ifile &f, imat &m)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "imat", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(m);
  return f;
}

it_ifile &operator>>(it_ifile &f, mat &m)
{
  it_file::data_header h;

  f.read_data_header(h);
  if (h.type == "fmat")
    f.low_level_read_lo(m);
  else if (h.type == "dmat")
    f.low_level_read_hi(m);
  else
    it_error("it_ifile::operator>>(): Wrong type");

  return f;
}

it_ifile &operator>>(it_ifile &f, cmat &m)
{
  it_file::data_header h;
  f.read_data_header(h);
  if (h.type == "fcmat")
    f.low_level_read_lo(m);
  else if (h.type == "dcmat")
    f.low_level_read_hi(m);
  else
    it_error("it_ifile::operator>>(): Wrong type");

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<bin> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "bArray", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(v);
  return f;
}

it_ifile &operator>>(it_ifile &f, Array<short> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "sArray", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(v);
  return f;
}

it_ifile &operator>>(it_ifile &f, Array<int> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "iArray", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(v);
  return f;
}

it_ifile &operator>>(it_ifile &f, Array<float> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "fArray", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(v);
  return f;
}

it_ifile &operator>>(it_ifile &f, Array<double> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  if (h.type == "fArray")
    f.low_level_read_lo(v);
  else if (h.type == "dArray")
    f.low_level_read_hi(v);
  else
    it_error("it_ifile::operator>>(): Wrong type");

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<std::complex<float> > &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "fcArray", "it_ifile::operator>>(): Wrong type");
  f.low_level_read(v);
  return f;
}

it_ifile &operator>>(it_ifile &f, Array<std::complex<double> > &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  if (h.type == "fcArray")
    f.low_level_read_lo(v);
  else if (h.type == "dcArray")
    f.low_level_read_hi(v);
  else
    it_error("it_ifile::operator>>(): Wrong type");

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<bvec> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "bvecArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<svec> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "svecArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<ivec> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "ivecArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<vec> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "vecArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read_hi(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<cvec> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "cvecArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read_hi(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<std::string> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "stringArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<bmat> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "bmatArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<smat> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "smatArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<imat> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "imatArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<mat> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "matArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read_hi(v(i));

  return f;
}

it_ifile &operator>>(it_ifile &f, Array<cmat> &v)
{
  it_file::data_header h;
  f.read_data_header(h);
  it_assert(h.type == "cmatArray", "it_ifile::operator>>(): Wrong type");
  uint64_t n;
  f.low_level_read(n);
  int size = static_cast<int>(n);
  v.set_size(size, false);
  for (int i = 0; i < size; ++i)
    f.low_level_read_hi(v(i));

  return f;
}


it_file &operator<<(it_file &f, char x)
{
  f.write_data_header("int8", sizeof(char));
  f.low_level_write(x);
  return f;
}

it_file &operator<<(it_file &f, bool x)
{
  f.write_data_header("bool", sizeof(char));
  f.low_level_write(x);
  return f;
}

it_file &operator<<(it_file &f, bin x)
{
  f.write_data_header("bin", sizeof(char));
  f.low_level_write(x);
  return f;
}

it_file &operator<<(it_file &f, short x)
{
  f.write_data_header("int16", sizeof(int16_t));
  f.low_level_write(x);
  return f;
}

it_file &operator<<(it_file &f, int x)
{
  f.write_data_header("int32", sizeof(int32_t));
  f.low_level_write(x);
  return f;
}

it_file &operator<<(it_file &f, float x)
{
  f.write_data_header("float32", sizeof(float));
  f.low_level_write(x);
  return f;
}

it_file &operator<<(it_file &f, double x)
{
  f.write_data_header("float64", sizeof(double));
  f.low_level_write(x);
  return f;
}

it_file &operator<<(it_file &f, std::complex<float> x)
{
  f.write_data_header("cfloat32", 2 * sizeof(float));
  f.low_level_write(x);
  return f;
}

it_file &operator<<(it_file &f, std::complex<double> x)
{
  f.write_data_header("cfloat64", 2 * sizeof(double));
  f.low_level_write(x);
  return f;
}

it_file &operator<<(it_file &f, const bvec &v)
{
  f.write_data_header("bvec", sizeof(uint64_t) + v.size() * sizeof(char));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const svec &v)
{
  f.write_data_header("svec", sizeof(uint64_t) + v.size() * sizeof(int16_t));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const ivec &v)
{
  f.write_data_header("ivec", sizeof(uint64_t) + v.size() * sizeof(int32_t));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const vec &v)
{
  if (f.get_low_precision())
    f.write_data_header("fvec", sizeof(uint64_t)
                        + v.size() * sizeof(float));
  else
    f.write_data_header("dvec", sizeof(uint64_t)
                        + v.size() * sizeof(double));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const cvec &v)
{
  if (f.get_low_precision())
    f.write_data_header("fcvec", sizeof(uint64_t)
                        + v.size() * 2 * sizeof(float));
  else
    f.write_data_header("dcvec", sizeof(uint64_t)
                        + v.size() * 2 * sizeof(double));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const std::string &str)
{
  f.write_data_header("string", sizeof(uint64_t) + str.size() * sizeof(char));
  f.low_level_write(str);
  return f;
}

it_file &operator<<(it_file &f, const bmat &m)
{
  f.write_data_header("bmat", 2 * sizeof(uint64_t)
                      + m.rows() * m.cols() * sizeof(char));
  f.low_level_write(m);
  return f;
}

it_file &operator<<(it_file &f, const smat &m)
{
  f.write_data_header("smat", 2 * sizeof(uint64_t)
                      + m.rows() * m.cols() * sizeof(int16_t));
  f.low_level_write(m);
  return f;
}

it_file &operator<<(it_file &f, const imat &m)
{
  f.write_data_header("imat", 2 * sizeof(uint64_t)
                      + m.rows() * m.cols() * sizeof(int32_t));
  f.low_level_write(m);
  return f;
}

it_file &operator<<(it_file &f, const mat &m)
{
  if (f.get_low_precision())
    f.write_data_header("fmat", 2 * sizeof(uint64_t)
                        + m.rows() * m.cols() * sizeof(float));
  else
    f.write_data_header("dmat", 2 * sizeof(uint64_t)
                        + m.rows() * m.cols() * sizeof(double));
  f.low_level_write(m);
  return f;
}

it_file &operator<<(it_file &f, const cmat &m)
{
  if (f.get_low_precision())
    f.write_data_header("fcmat", 2 * sizeof(uint64_t)
                        + m.rows() * m.cols() * 2 * sizeof(float));
  else
    f.write_data_header("dcmat", 2 * sizeof(uint64_t)
                        + m.rows() * m.cols() * 2 * sizeof(double));
  f.low_level_write(m);
  return f;
}

it_file &operator<<(it_file &f, const Array<bin> &v)
{
  f.write_data_header("bArray", sizeof(uint64_t) + v.size() * sizeof(char));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const Array<short> &v)
{
  f.write_data_header("sArray", sizeof(uint64_t)
                      + v.size() * sizeof(int16_t));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const Array<int> &v)
{
  f.write_data_header("iArray", sizeof(uint64_t)
                      + v.size() * sizeof(int32_t));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const Array<float> &v)
{
  f.write_data_header("fArray", sizeof(uint64_t) + v.size() * sizeof(float));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const Array<double> &v)
{
  if (f.get_low_precision())
    f.write_data_header("fArray", sizeof(uint64_t)
                        + v.size() * sizeof(float));
  else
    f.write_data_header("dArray", sizeof(uint64_t)
                        + v.size() * sizeof(double));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const Array<std::complex<float> > &v)
{
  f.write_data_header("fcArray", sizeof(uint64_t)
                      + v.size() * 2 * sizeof(float));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const Array<std::complex<double> > &v)
{
  if (f.get_low_precision())
    f.write_data_header("fcArray", sizeof(uint64_t)
                        + v.size() * 2 * sizeof(float));
  else
    f.write_data_header("dcArray", sizeof(uint64_t)
                        + v.size() * 2 * sizeof(double));
  f.low_level_write(v);
  return f;
}

it_file &operator<<(it_file &f, const Array<bvec> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i).size();

  // write header
  f.write_data_header("bvecArray", sizeof(uint64_t) * (1 + v.size())
                      + sum_l * sizeof(char));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<svec> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i).size();

  // write header
  f.write_data_header("svecArray", sizeof(uint64_t) * (1 + v.size())
                      + sum_l * sizeof(int16_t));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<ivec> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i).size();

  // write header
  f.write_data_header("ivecArray", sizeof(uint64_t) * (1 + v.size())
                      + sum_l * sizeof(int32_t));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<vec> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i).size();

  // write header
  f.write_data_header("vecArray", sizeof(uint64_t) * (1 + v.size())
                      + sum_l * sizeof(double));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<cvec> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i).size();

  // write header
  f.write_data_header("cvecArray", sizeof(uint64_t) * (1 + v.size())
                      + sum_l * 2 * sizeof(double));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<std::string> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += int(v(i).size());

  // write header
  f.write_data_header("stringArray", sizeof(uint64_t) * (1 + v.size())
                      + sum_l * sizeof(char));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<bmat> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i)._datasize();

  // write header
  f.write_data_header("bmatArray", sizeof(uint64_t) * (1 + 2 * v.size())
                      + sum_l * sizeof(char));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<smat> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i)._datasize();

  // write header
  f.write_data_header("smatArray", sizeof(uint64_t) * (1 + 2 * v.size())
                      + sum_l * sizeof(int16_t));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<imat> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i)._datasize();

  // write header
  f.write_data_header("imatArray", sizeof(uint64_t) * (1 + 2 * v.size())
                      + sum_l * sizeof(int32_t));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<mat> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i)._datasize();

  // write header
  f.write_data_header("matArray", sizeof(uint64_t) * (1 + 2 * v.size())
                      + sum_l * sizeof(double));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}

it_file &operator<<(it_file &f, const Array<cmat> &v)
{
  // calculate total length of Array
  int sum_l = 0;
  for (int i = 0; i < v.size(); ++i)
    sum_l += v(i)._datasize();

  // write header
  f.write_data_header("cmatArray", sizeof(uint64_t) * (1 + 2 * v.size())
                      + sum_l * 2 * sizeof(double));
  // write the length of the array
  f.low_level_write(static_cast<uint64_t>(v.size()));

  // write one vector at a time (i.e. size and elements)
  for (int i = 0; i < v.size(); ++i)
    f.low_level_write(v(i));

  return f;
}



// ----------------------------------------------------------------------
// Deprecated implementation of IT++ file format version 2
// Will be removed in future versions
// ----------------------------------------------------------------------

char it_file_base_old::file_magic[4] = { 'I', 'T', '+', '+' };
char it_file_base_old::file_version = 2;

it_ifile_old::it_ifile_old()
{
}

it_ifile_old::it_ifile_old(const std::string &name)
{
  open(name);
}

void it_ifile_old::open(const std::string &name)
{
  it_assert(exist(name), "File does not exist");

  s.open_readonly(name);

  if (!read_check_file_header()) {
    s.close();
    it_error("Corrupt file (Not an it-file)");
  }

}

void it_ifile_old::close()
{
  s.close();
}

bool it_ifile_old::seek(const std::string &name)
{
  data_header h;
  std::streampos p;

  s.clear();
  s.seekg(sizeof(file_header));

  while (true) {
    p = s.tellg();
    read_data_header(h);
    if (s.eof()) {
      s.clear();
      return false;
    }
    if (h.type != "" && h.name == name) {
      s.seekg(p);
      break;
    }
    s.seekg(p + static_cast<std::streamoff>(h.block_bytes));
  }

  return true;
}

bool it_ifile_old::seek(int n)
{
  data_header h;
  std::streampos p;

  s.clear();
  s.seekg(sizeof(file_header));
  for (int i = 0; i <= n; i++) {
    p = s.tellg(); // changed from tellp() since probably an error
    read_data_header(h);
    if (s.eof()) {
      s.clear();
      return false;
    }
    if (h.type == "")
      i--;
    s.seekg(i == n ? p : p + static_cast<std::streamoff>(h.block_bytes));
  }
  return true;
}

void it_ifile_old::info(std::string &name, std::string &type, int &bytes)
{
  data_header h;
  std::streampos p;

  p = s.tellg(); // changed from tellp()
  read_data_header(h);
  s.seekg(p);
  name = h.name;
  type = h.type;
  bytes = h.data_bytes;
}

bool it_ifile_old::read_check_file_header()
{
  file_header h;

  memset(&h, 0, sizeof(h)); // Clear the struct
  s.read(reinterpret_cast<char *>(&h), sizeof(h));

  return (memcmp(h.magic, file_magic, 4) == 0 && (h.version <= file_version));
}

void it_ifile_old::read_data_header(data_header &h)
{
  std::streampos p = s.tellg();
  s.clear();
  s >> h.endianity;
  if (s.eof())
    return;
  s.set_endianity(static_cast<bfstream_base::endian>(h.endianity));
  uint32_t tmp;
  s >> tmp;
  h.hdr_bytes = tmp;
  s >> tmp;
  h.data_bytes = tmp;
  s >> tmp;
  h.block_bytes = tmp;
  s >> h.name;
  s >> h.type;
  s.seekg(p + static_cast<std::streamoff>(h.hdr_bytes));
}

void it_ifile_old::low_level_read(char &x)
{
  s >> x;
}

void it_ifile_old::low_level_read(bin &x)
{
  s >> x;
}

void it_ifile_old::low_level_read(short &x)
{
  s >> x;
}

void it_ifile_old::low_level_read(int &x)
{
  int32_t tmp;
  s >> tmp;
  x = tmp;
}

void it_ifile_old::low_level_read(float &x)
{
  s >> x;
}

void it_ifile_old::low_level_read(double &x)
{
  s >> x;
}

void it_ifile_old::low_level_read(std::complex<float> &x)
{
  float x_real, x_imag;
  s >> x_real;
  s >> x_imag;
  x = std::complex<float>(x_real, x_imag);
}

void it_ifile_old::low_level_read(std::complex<double> &x)
{
  double x_real, x_imag;
  s >> x_real;
  s >> x_imag;
  x = std::complex<double>(x_real, x_imag);
}

void it_ifile_old::low_level_read_lo(vec &v)
{
  int32_t i;
  float val;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val;
    v(i) = static_cast<double>(val);
  }
}

void it_ifile_old::low_level_read_hi(vec &v)
{
  int32_t i;
  double val;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val;
    v(i) = static_cast<double>(val);
  }
}

void it_ifile_old::low_level_read(ivec &v)
{
  int32_t i, val;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val;
    v(i) = val;
  }
}

void it_ifile_old::low_level_read(bvec &v)
{
  int32_t i;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++)
    s >> v(i);
}

void it_ifile_old::low_level_read_lo(cvec &v)
{
  int32_t i;
  float val_real, val_imag;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<double>(val_real, val_imag);
  }
}

void it_ifile_old::low_level_read_hi(cvec &v)
{
  int32_t i;
  double val_real, val_imag;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<double>(val_real, val_imag);
  }
}

void it_ifile_old::low_level_read(std::string &str)
{
  int32_t i, j;
  char val;
  str = "";

  s >> i;

  for (j = 0; j < i; j++) {
    s >> val;
    str += val;
  }
}

void it_ifile_old::low_level_read_lo(mat &m)
{
  int32_t i, j;
  float val;

  s >> i >> j;
  m.set_size(i, j, false);
  for (j = 0; j < m.cols(); j++)
    for (i = 0; i < m.rows(); i++) {
      s >> val;
      m(i, j) = static_cast<double>(val);
    }
}

void it_ifile_old::low_level_read_hi(mat &m)
{
  int32_t i, j;
  double val;

  s >> i >> j;
  m.set_size(i, j, false);
  for (j = 0; j < m.cols(); j++)
    for (i = 0; i < m.rows(); i++) {
      s >> val;
      m(i, j) = static_cast<double>(val);
    }
}

void it_ifile_old::low_level_read(imat &m)
{
  int32_t i, j, val;

  s >> i >> j;
  m.set_size(i, j, false);
  for (j = 0; j < m.cols(); j++)
    for (i = 0; i < m.rows(); i++) {
      s >> val;
      m(i, j) = val;
    }
}

void it_ifile_old::low_level_read(bmat &m)
{
  int32_t i, j;

  s >> i >> j;
  m.set_size(i, j, false);
  for (j = 0; j < m.cols(); j++)
    for (i = 0; i < m.rows(); i++)
      s >> m(i, j);
}

void it_ifile_old::low_level_read_lo(cmat &m)
{
  int32_t i, j;
  float val_real, val_imag;

  s >> i >> j;
  m.set_size(i, j, false);
  for (j = 0; j < m.cols(); j++)
    for (i = 0; i < m.rows(); i++) {
      s >> val_real;
      s >> val_imag;
      m(i, j) = std::complex<double>(val_real, val_imag);
    }
}

void it_ifile_old::low_level_read_hi(cmat &m)
{
  int32_t i, j;
  double val_real, val_imag;

  s >> i >> j;
  m.set_size(i, j, false);
  for (j = 0; j < m.cols(); j++)
    for (i = 0; i < m.rows(); i++) {
      s >> val_real;
      s >> val_imag;
      m(i, j) = std::complex<double>(val_real, val_imag);
    }
}


void it_ifile_old::low_level_read_lo(Array<float> &v)
{
  int32_t i;
  float val;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val;
    v(i) = val;
  }
}

void it_ifile_old::low_level_read_lo(Array<double> &v)
{
  int32_t i;
  float val;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val;
    v(i) = static_cast<double>(val);
  }
}

void it_ifile_old::low_level_read_hi(Array<double> &v)
{
  int32_t i;
  double val;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val;
    v(i) = static_cast<double>(val);
  }
}

void it_ifile_old::low_level_read(Array<int> &v)
{
  int32_t i, val;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val;
    v(i) = val;
  }
}

void it_ifile_old::low_level_read(Array<bin> &v)
{
  int32_t i;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++)
    s >> v(i);
}

void it_ifile_old::low_level_read_lo(Array<std::complex<float> > &v)
{
  int32_t i;
  float val_real, val_imag;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<float>(val_real, val_imag);
  }
}

void it_ifile_old::low_level_read_lo(Array<std::complex<double> > &v)
{
  int32_t i;
  float val_real, val_imag;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<double>(val_real, val_imag);
  }
}

void it_ifile_old::low_level_read_hi(Array<std::complex<double> > &v)
{
  int32_t i;
  double val_real, val_imag;

  s >> i;
  v.set_size(i, false);
  for (i = 0; i < v.size(); i++) {
    s >> val_real;
    s >> val_imag;
    v(i) = std::complex<double>(val_real, val_imag);
  }
}

it_file_old::it_file_old():_string(new String_Holder())
{
  low_prec = false;
}

it_file_old::it_file_old(const std::string &name, bool trunc):
_string(new String_Holder())
{
  low_prec = false;
  open(name, trunc);
}

void it_file_old::open(const std::string &name, bool trunc)
{
  if (!exist(name))
    trunc = true;

  s.open(name, trunc);
  it_error_if(!s.is_open(), "Could not open file for writing");

  if (trunc)
    write_file_header();
  else if (!read_check_file_header()) {
    s.close();
    it_error("Corrupt file (Not an it-file)");
  }
}

void it_file_old::close()
{
  s.close();
}

void it_file_old::flush()
{
  s.flush();
}

void it_file_old::write_file_header()
{
  s.write(file_magic, 4);
  s << file_version;
}

void it_file_old::write_data_header(const std::string &type, uint32_t size)
{
  it_error_if(next_name() == "", "Try to write without a name");
  write_data_header(type, next_name(), size);
  next_name() = "";
}

void it_file_old::write_data_header(const std::string &type,
                                    const std::string &name, uint32_t size)
{
  data_header h1, h2;
  std::streampos p;
  int availpos = 0;
  bool removed = false;
  int skip;

  h1.endianity = static_cast<char>(s.get_native_endianity());
  h1.hdr_bytes = 1 + 3 * 4 + int(type.size()) + 1 + int(name.size()) + 1;
  h1.data_bytes = size;
  h1.block_bytes = h1.hdr_bytes + h1.data_bytes;
  h1.name = name;
  h1.type = type;

  if (exists(name))
    remove();

  // Try to find an empty space
  s.clear();
  s.seekg(sizeof(file_header));
  while (true) {
    p = s.tellp();
    read_data_header(h2);
    if (s.eof()) {
      s.clear();
      break;
    }
    skip = h2.block_bytes;
    if (h2.type != "" && h2.name == name) {
      s.seekg(p);
      remove();
      s.seekg(p);
      read_data_header(h2);
      removed = true;
      if (availpos != 0)
        break;
    }
    if (availpos == 0) {
      if (h2.type == "" && h2.block_bytes >= h1.block_bytes) {
        h1.block_bytes = h2.block_bytes;
        availpos = int(p);
      }
      else if (h2.block_bytes - h2.hdr_bytes - h2.data_bytes >= h1.block_bytes) {
        h1.block_bytes = h2.block_bytes - h2.hdr_bytes - h2.data_bytes;
        h2.block_bytes = h2.hdr_bytes + h2.data_bytes;
        s.seekp(p);
        write_data_header_here(h2);
        availpos = static_cast<int>(p) + h2.block_bytes;
        if (removed)
          break;
      }
    }
    s.seekg(p + static_cast<std::streamoff>(skip));
  }
  if (availpos != 0)
    s.seekp(availpos);
  else
    s.seekp(0, std::ios::end);

  write_data_header_here(h1);
}

void it_file_old::write_data_header_here(const data_header &h)
{
  s.set_endianity(static_cast<bfstream_base::endian>(h.endianity));
  s << h.endianity << h.hdr_bytes << h.data_bytes << h.block_bytes << h.name << h.type;
}

void it_file_old::remove(const std::string &name)
{
  seek(name);
  remove();
}

void it_file_old::remove()
{
  data_header h;
  std::streampos p;

  p = s.tellp();
  read_data_header(h);
  h.type = "";
  h.name = "";
  h.hdr_bytes = 1 + 3 * 4 + 1 + 1;
  h.data_bytes = 0;
  s.seekp(p);
  write_data_header_here(h);
  s.seekp(p + static_cast<std::streamoff>(h.block_bytes));
}

bool it_file_old::exists(const std::string &name)
{
  if (seek(name))
    return true;
  else
    return false;
}

void it_file_old::pack()
{
  it_warning("pack() is not implemented!");
}

void it_file_old::low_level_write(char x)
{
  s << x;
}

void it_file_old::low_level_write(bin x)
{
  s << x.value();
}

void it_file_old::low_level_write(short x)
{
  s << x;
}

void it_file_old::low_level_write(int x)
{
  s << static_cast<int32_t>(x);
}

void it_file_old::low_level_write(float x)
{
  s << x;
}

void it_file_old::low_level_write(double x)
{
  s << x;
}

void it_file_old::low_level_write(const std::complex<float> &x)
{
  s << x.real();
  s << x.imag();
}

void it_file_old::low_level_write(const std::complex<double> &x)
{
  s << x.real();
  s << x.imag();
}

void it_file_old::low_level_write(const vec &v)
{
  if (get_low_precision()) {
    s << static_cast<int32_t>(v.size());
    for (int i = 0; i < v.size(); i++)
      s << static_cast<float>(v(i));
  }
  else {
    s << static_cast<int32_t>(v.size());
    for (int i = 0; i < v.size(); i++)
      s << static_cast<double>(v(i));
  }
}

void it_file_old::low_level_write(const ivec &v)
{
  s << static_cast<int32_t>(v.size());
  for (int i = 0; i < v.size(); i++)
    s << static_cast<int32_t>(v(i));
}

void it_file_old::low_level_write(const bvec &v)
{
  s << static_cast<int32_t>(v.size());
  for (int i = 0; i < v.size(); i++)
    s << v(i).value();
}

void it_file_old::low_level_write(const cvec &v)
{
  if (get_low_precision()) {
    s << static_cast<int32_t>(v.size());
    for (int i = 0; i < v.size(); i++) {
      s << static_cast<float>(v(i).real());
      s << static_cast<float>(v(i).imag());
    }
  }
  else {
    s << static_cast<int32_t>(v.size());
    for (int i = 0; i < v.size(); i++) {
      s << static_cast<double>(v(i).real());
      s << static_cast<double>(v(i).imag());
    }
  }
}

void it_file_old::low_level_write(const std::string &str)
{
  int size = int(str.size());
  s << static_cast<int32_t>(size);

  for (int i = 0; i < size; i++)
    s << str[i];
}

void it_file_old::low_level_write(const mat &m)
{
  int i, j;

  if (get_low_precision()) {
    s << static_cast<int32_t>(m.rows()) << static_cast<int32_t>(m.cols());
    for (j = 0; j < m.cols(); j++)
      for (i = 0; i < m.rows(); i++)
        s << static_cast<float>(m(i, j));
  }
  else {
    s << static_cast<int32_t>(m.rows()) << static_cast<int32_t>(m.cols());
    for (j = 0; j < m.cols(); j++)
      for (i = 0; i < m.rows(); i++)
        s << static_cast<double>(m(i, j));
  }
}

void it_file_old::low_level_write(const imat &m)
{
  int i, j;

  s << static_cast<int32_t>(m.rows()) << static_cast<int32_t>(m.cols());
  for (j = 0; j < m.cols(); j++)
    for (i = 0; i < m.rows(); i++)
      s << static_cast<int32_t>(m(i, j));
}

void it_file_old::low_level_write(const bmat &m)
{
  int i, j;

  s << static_cast<int32_t>(m.rows()) << static_cast<int32_t>(m.cols());
  for (j = 0; j < m.cols(); j++)
    for (i = 0; i < m.rows(); i++)
      s << m(i, j).value();
}

void it_file_old::low_level_write(const cmat &m)
{
  int i, j;

  if (get_low_precision()) {
    s << static_cast<int32_t>(m.rows()) << static_cast<int32_t>(m.cols());
    for (j = 0; j < m.cols(); j++)
      for (i = 0; i < m.rows(); i++) {
        s << static_cast<float>(m(i, j).real());
        s << static_cast<float>(m(i, j).imag());
      }

  }
  else {
    s << static_cast<int32_t>(m.rows()) << static_cast<int32_t>(m.cols());
    for (j = 0; j < m.cols(); j++)
      for (i = 0; i < m.rows(); i++) {
        s << static_cast<double>(m(i, j).real());
        s << static_cast<double>(m(i, j).imag());
      }
  }
}

void it_file_old::low_level_write(const Array<float> &v)
{
  s << static_cast<int32_t>(v.size());
  for (int i = 0; i < v.size(); i++)
    s << v(i);
}

void it_file_old::low_level_write(const Array<double> &v)
{
  if (get_low_precision()) {
    s << static_cast<int32_t>(v.size());
    for (int i = 0; i < v.size(); i++)
      s << static_cast<float>(v(i));
  }
  else {
    s << static_cast<int32_t>(v.size());
    for (int i = 0; i < v.size(); i++)
      s << static_cast<double>(v(i));
  }
}

void it_file_old::low_level_write(const Array<int> &v)
{
  s << static_cast<int32_t>(v.size());
  for (int i = 0; i < v.size(); i++)
    s << static_cast<int32_t>(v(i));
}

void it_file_old::low_level_write(const Array<bin> &v)
{
  s << static_cast<int32_t>(v.size());
  for (int i = 0; i < v.size(); i++)
    s << v(i).value();
}

void it_file_old::low_level_write(const Array<std::complex<float> > &v)
{
  s << static_cast<int32_t>(v.size());
  for (int i = 0; i < v.size(); i++) {
    s << v(i).real();
    s << v(i).imag();
  }
}

void it_file_old::low_level_write(const Array<std::complex<double> > &v)
{
  if (get_low_precision()) {
    s << static_cast<int32_t>(v.size());
    for (int i = 0; i < v.size(); i++) {
      s << static_cast<float>(v(i).real());
      s << static_cast<float>(v(i).imag());
    }
  }
  else {
    s << static_cast<int32_t>(v.size());
    for (int i = 0; i < v.size(); i++) {
      s << static_cast<double>(v(i).real());
      s << static_cast<double>(v(i).imag());
    }
  }
}

it_ifile_old &operator>>(it_ifile_old &f, char &x)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "int8")
    f.low_level_read(x);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, bin &x)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "bin")
    f.low_level_read(x);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, short &x)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "int16")
    f.low_level_read(x);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, int &x)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "int32")
    f.low_level_read(x);
  else if (h.type == "int16") {
    short x16;
    f.low_level_read(x16);
    x = x16;
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, double &x)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "float64")
    f.low_level_read(x);
  else if (h.type == "float32") {
    float f32;
    f.low_level_read(f32);
    x = f32;
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, float &x)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "float32")
    f.low_level_read(x);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, std::complex<float> &x)
{
  it_file_old::data_header h;

  f.read_data_header(h);

  if (h.type == "float32_complex") {
    std::complex<float> f32_c;
    f.low_level_read(f32_c);
    x = f32_c;
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, std::complex<double> &x)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "float64_complex")
    f.low_level_read(x);
  else if (h.type == "float32_complex") {
    std::complex<float> f32_c;
    f.low_level_read(f32_c);
    x = f32_c;
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, vec &v)
{
  it_ifile_old::data_header h;

  f.read_data_header(h);
  if (h.type == "fvec")
    f.low_level_read_lo(v);
  else if (h.type == "dvec")
    f.low_level_read_hi(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, ivec &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "ivec")
    f.low_level_read(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, bvec &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "bvec")
    f.low_level_read(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, cvec &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "fcvec")
    f.low_level_read_lo(v);
  else if (h.type == "dcvec")
    f.low_level_read_hi(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, std::string &str)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "string")
    f.low_level_read(str);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, mat &m)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "fmat")
    f.low_level_read_lo(m);
  else if (h.type == "dmat")
    f.low_level_read_hi(m);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, imat &m)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "imat")
    f.low_level_read(m);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, bmat &m)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "bmat")
    f.low_level_read(m);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, cmat &m)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "fcmat")
    f.low_level_read_lo(m);
  else if (h.type == "dcmat")
    f.low_level_read_hi(m);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<float> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "fArray")
    f.low_level_read_lo(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<double> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "fArray")
    f.low_level_read_lo(v);
  else if (h.type == "dArray")
    f.low_level_read_hi(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<int> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "iArray")
    f.low_level_read(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<bin> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "bArray")
    f.low_level_read(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<std::complex<float> > &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "fcArray")
    f.low_level_read_lo(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<std::complex<double> > &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "fcArray")
    f.low_level_read_lo(v);
  else if (h.type == "dcArray")
    f.low_level_read_hi(v);
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<vec> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "vecArray") {
    int n;
    f.low_level_read(n);
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      f.low_level_read_hi(v(i));
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<ivec> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "ivecArray") {
    int n;
    f.low_level_read(n);
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      f.low_level_read(v(i));
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<bvec> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "bvecArray") {
    int n;
    f.low_level_read(n);
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      f.low_level_read(v(i));
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<cvec> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "cvecArray") {
    int n;
    f.low_level_read(n);
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      f.low_level_read_hi(v(i));
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<std::string> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "stringArray") {
    int n;
    f.low_level_read(n);
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      f.low_level_read(v(i));
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<mat> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "matArray") {
    int n;
    f.low_level_read(n);
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      f.low_level_read_hi(v(i));
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<imat> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "imatArray") {
    int n;
    f.low_level_read(n);
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      f.low_level_read(v(i));
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<bmat> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "bmatArray") {
    int n;
    f.low_level_read(n);
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      f.low_level_read(v(i));
  }
  else
    it_error("Wrong type");

  return f;
}

it_ifile_old &operator>>(it_ifile_old &f, Array<cmat> &v)
{
  it_file_old::data_header h;

  f.read_data_header(h);
  if (h.type == "cmatArray") {
    int n;
    f.low_level_read(n);
    v.set_size(n, false);
    for (int i = 0; i < n; i++)
      f.low_level_read_hi(v(i));
  }
  else
    it_error("Wrong type");

  return f;
}

it_file_old &operator<<(it_file_old &f, char x)
{
  f.write_data_header("int8", sizeof(char));
  f.low_level_write(x);

  return f;
}

it_file_old &operator<<(it_file_old &f, bin x)
{
  f.write_data_header("bin", sizeof(bin));
  f.low_level_write(x);

  return f;
}

it_file_old &operator<<(it_file_old &f, short x)
{
  f.write_data_header("int16", sizeof(short));
  f.low_level_write(x);

  return f;
}

it_file_old &operator<<(it_file_old &f, int x)
{
  f.write_data_header("int32", sizeof(int));
  f.low_level_write(x);

  return f;
}

it_file_old &operator<<(it_file_old &f, float x)
{
  f.write_data_header("float32", sizeof(float));
  f.low_level_write(x);

  return f;
}

it_file_old &operator<<(it_file_old &f, double x)
{
  f.write_data_header("float64", sizeof(double));
  f.low_level_write(x);

  return f;
}

it_file_old &operator<<(it_file_old &f, std::complex<float> x)
{
  f.write_data_header("float32_complex", 2*sizeof(float));
  f.low_level_write(x);

  return f;
}

it_file_old &operator<<(it_file_old &f, std::complex<double> x)
{
  f.write_data_header("float64_complex", 2*sizeof(double));
  f.low_level_write(x);

  return f;
}

it_file_old &operator<<(it_file_old &f, const vec &v)
{
  if (f.get_low_precision())
    f.write_data_header("fvec", sizeof(int) + v.size() * sizeof(float));
  else
    f.write_data_header("dvec", sizeof(int) + v.size() * sizeof(double));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const ivec &v)
{
  f.write_data_header("ivec", (1 + v.size()) * sizeof(int));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const bvec &v)
{
  f.write_data_header("bvec", sizeof(int) + v.size() * sizeof(bin));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const cvec &v)
{
  if (f.get_low_precision())
    f.write_data_header("fcvec", sizeof(int) + v.size() * 2 * sizeof(float));
  else
    f.write_data_header("dcvec", sizeof(int) + v.size() * 2 * sizeof(double));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const std::string &str)
{
  f.write_data_header("string", sizeof(int) + int(str.size()) * sizeof(char));
  f.low_level_write(str);

  return f;
}

it_file_old &operator<<(it_file_old &f, const mat &m)
{
  if (f.get_low_precision())
    f.write_data_header("fmat", 2*sizeof(int) + m.rows()*m.cols()*sizeof(float));
  else
    f.write_data_header("dmat", 2*sizeof(int) + m.rows()*m.cols()*sizeof(double));
  f.low_level_write(m);

  return f;
}

it_file_old &operator<<(it_file_old &f, const imat &m)
{
  f.write_data_header("imat", (2 + m.rows()*m.cols()) * sizeof(int));
  f.low_level_write(m);

  return f;
}

it_file_old &operator<<(it_file_old &f, const bmat &m)
{
  f.write_data_header("bmat", 2*sizeof(int) + m.rows()*m.cols()*sizeof(bin));
  f.low_level_write(m);

  return f;
}

it_file_old &operator<<(it_file_old &f, const cmat &m)
{
  if (f.get_low_precision())
    f.write_data_header("fcmat", 2*sizeof(int) + 2*m.rows()*m.cols()*sizeof(float));
  else
    f.write_data_header("dcmat", 2*sizeof(int) + 2*m.rows()*m.cols()*sizeof(double));
  f.low_level_write(m);

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<float> &v)
{
  f.write_data_header("fArray", sizeof(int) + v.size() * sizeof(float));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<double> &v)
{
  if (f.get_low_precision())
    f.write_data_header("fArray", sizeof(int) + v.size() * sizeof(float));
  else
    f.write_data_header("dArray", sizeof(int) + v.size() * sizeof(double));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<int> &v)
{
  f.write_data_header("iArray", (1 + v.size()) * sizeof(int));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<bin> &v)
{
  f.write_data_header("bArray", sizeof(int) + v.size() * sizeof(bin));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<std::complex<float> > &v)
{
  f.write_data_header("fcArray", sizeof(int) + v.size() * 2 * sizeof(float));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<std::complex<double> > &v)
{
  if (f.get_low_precision())
    f.write_data_header("fcArray", sizeof(int) + v.size() * 2 * sizeof(float));
  else
    f.write_data_header("dcArray", sizeof(int) + v.size() * 2 * sizeof(double));
  f.low_level_write(v);

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<vec> &v)
{
  int i, sum_l = 0;

  // calculate total length of Array
  for (i = 0; i < v.size(); i++) {
    sum_l += v(i).size();
  }

  // write header
  f.write_data_header("vecArray", sizeof(int)*(1 + v.size()) + sum_l * sizeof(double));

  f.low_level_write(v.size());  // the length of the array

  // write one vector at a time (i.e. size and elements)
  for (i = 0; i < v.size(); i++)
    f.low_level_write(v(i));

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<ivec> &v)
{
  int i, sum_l = 0;

  // calculate total length of Array
  for (i = 0; i < v.size(); i++) {
    sum_l += v(i).size();
  }

  // write header
  f.write_data_header("ivecArray", sizeof(int)*(1 + v.size()) + sum_l * sizeof(int));

  f.low_level_write(v.size());  // the length of the array

  // write one vector at a time (i.e. size and elements)
  for (i = 0; i < v.size(); i++)
    f.low_level_write(v(i));

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<bvec> &v)
{
  int i, sum_l = 0;

  // calculate total length of Array
  for (i = 0; i < v.size(); i++) {
    sum_l += v(i).size();
  }

  // write header
  f.write_data_header("bvecArray", sizeof(int)*(1 + v.size()) + sum_l * sizeof(bin));

  f.low_level_write(v.size());  // the length of the array

  // write one vector at a time (i.e. size and elements)
  for (i = 0; i < v.size(); i++)
    f.low_level_write(v(i));

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<cvec> &v)
{
  int i, sum_l = 0;

  // calculate total length of Array
  for (i = 0; i < v.size(); i++) {
    sum_l += v(i).size();
  }

  // write header
  f.write_data_header("cvecArray", sizeof(int)*(1 + v.size()) + sum_l * sizeof(std::complex<double>));

  f.low_level_write(v.size());  // the length of the array

  // write one vector at a time (i.e. size and elements)
  for (i = 0; i < v.size(); i++)
    f.low_level_write(v(i));

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<std::string> &v)
{
  int i, sum_l = 0;

  // calculate total length of Array
  for (i = 0; i < v.size(); i++) {
    sum_l += int(v(i).size());
  }

  // write header
  f.write_data_header("stringArray", sizeof(int)*(1 + v.size()) + sum_l * sizeof(char));

  f.low_level_write(v.size());  // the length of the array

  // write one vector at a time (i.e. size and elements)
  for (i = 0; i < v.size(); i++)
    f.low_level_write(v(i));

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<mat> &v)
{
  int i, sum_l = 0;

  // calculate total length of Array
  for (i = 0; i < v.size(); i++) {
    sum_l += v(i)._datasize();
  }

  // write header
  f.write_data_header("matArray", sizeof(int)*(1 + 2*v.size()) + sum_l * sizeof(double));

  f.low_level_write(v.size());  // the length of the array

  // write one matrix at a time (i.e. size and elements)
  for (i = 0; i < v.size(); i++)
    f.low_level_write(v(i));

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<imat> &v)
{
  int i, sum_l = 0;

  // calculate total length of Array
  for (i = 0; i < v.size(); i++) {
    sum_l += v(i)._datasize();
  }

  // write header
  f.write_data_header("imatArray", sizeof(int)*(1 + 2*v.size()) + sum_l * sizeof(int));

  f.low_level_write(v.size());  // the length of the array

  // write one matrix at a time (i.e. size and elements)
  for (i = 0; i < v.size(); i++)
    f.low_level_write(v(i));

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<bmat> &v)
{
  int i, sum_l = 0;

  // calculate total length of Array
  for (i = 0; i < v.size(); i++) {
    sum_l += v(i)._datasize();
  }

  // write header
  f.write_data_header("bmatArray", sizeof(int)*(1 + 2*v.size()) + sum_l * sizeof(bin));

  f.low_level_write(v.size());  // the length of the array

  // write one matrix at a time (i.e. size and elements)
  for (i = 0; i < v.size(); i++)
    f.low_level_write(v(i));

  return f;
}

it_file_old &operator<<(it_file_old &f, const Array<cmat> &v)
{
  int i, sum_l = 0;

  // calculate total length of Array
  for (i = 0; i < v.size(); i++) {
    sum_l += v(i)._datasize();
  }

  // write header
  f.write_data_header("cmatArray", sizeof(int)*(1 + 2*v.size()) + sum_l * sizeof(std::complex<double>));

  f.low_level_write(v.size());  // the length of the array

  // write one matrix at a time (i.e. size and elements)
  for (i = 0; i < v.size(); i++)
    f.low_level_write(v(i));

  return f;
}

} // namespace itpp
