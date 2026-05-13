/*!
 * \file
 * \brief Binary file formats definitions
 * \author Tony Ottosson, Thomas Eriksson, Adam Piatyszek and Andy Panov
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2013  (see AUTHORS file for a list of contributors)
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

#ifndef BINFILE_H
#define BINFILE_H

#include <itpp/itexports.h>

#include <itpp/base/ittypes.h>
#include <itpp/base/binary.h>
#include <fstream>

#ifdef _MSC_VER
#pragma warning( disable : 4250 )
#endif

namespace itpp
{

/*!
  \brief Checks if a file named \c name already exists on the disk
  \ingroup itfile
*/
bool exist(const std::string& name);

/*!
  \brief Base class for binary file classes
  \ingroup itfile

  This class serves as a base class for the classes \c bofstream,
  \c bifstream, and \c bfstream. It controls the endianity (i.e. the
  byte order of multibyte numbers on the disk) of the inhereted classes.
*/
class ITPP_EXPORT bfstream_base
{
public:
  /*!
    \brief Definition of the endian data type

    The Endianness defines the order in which multibyte numbers are stored
    in the file. The two orders are called "Little Endian" (\c l_endian )
    and "Big Endian" (\c b_endian ).

    "Little Endian" means that the low-order byte of the number is stored
    at the lowest address (i.e. the little end comes first). "Big Endian"
    means that the high-order byte of the number is stored in memory at
    the lowest address (i.e. the big end comes first)
  */
  enum endian { l_endian, b_endian };

  /*!
    \brief Class Constructor

    \param e Defines the endianity of the class. Possible values are \c
    l_endian for little endian or \c b_endian for big endian. The default
    value is \c b_endian.
  */
  bfstream_base(endian e = b_endian);

  /*!
    \brief Returns the endianity of the class
  */
  endian get_endianity() const {
    if (switch_endianity) {
      if (native_endianity == l_endian)
        return b_endian;
      else
        return l_endian;
    }
    else
      return native_endianity;
  }

  /*!
    \brief Returns the native endianity for this computer architecture

    Intel processors use "Little Endian" byte ordering while e.g. Motorola
    processors use "Big Endian" byte ordering.
  */
  endian get_native_endianity() const { return native_endianity; }

  /*!
    \brief Set the endianity for this class
  */
  void set_endianity(endian e) {
    if (native_endianity == e)
      switch_endianity = false;
    else
      switch_endianity = true;
  }

  /*!
    \brief Set the endianity of this class to the native endianity for
    this computer architecture
  */
  void set_native_endianity() { switch_endianity = false; }

protected:
  //! Indicates if the endianity of the processed data needs to be changed
  bool switch_endianity;
  //! The native endianity for this computer architecture
  endian native_endianity;
};

namespace binfile_details
{
/*!
  \brief Ofstream Interface Facade for Binary Streams.
  \ingroup itfile

  This class implements std::ofstream facade to make ITPP binary file streams exportable from dll.
  This facade implements basic functionality only. It does not provide an access
  to the following stream facilities (all of them are useless for binary streams)
  1. locale(imbue) It does not make sence to change locale settings for binary streams.
  Changes in formatting or char conversion can result in compatibility problems
  with resulting binary files
  2. stream buffer (rdbuf). DLL and application can use different versions of runtime ,
  so it would be dangerous to use buffer created in DLL in application context
  3. stream insertion operators. It is assumed that stream insertion is defined in binary stream
  classes derived from this class
  4. formatting interface (copyfmt, fill, narrow, widen). This is not relevant to binary streams
  5. ios_base-related stuff. These things are excluded since they provide unnecessarily formatting
  facilities.
*/
  class ITPP_EXPORT Ofstream_Binfile_Facade
  {
    std::ofstream* _str;
    //following makes class non-copiable
    Ofstream_Binfile_Facade(const Ofstream_Binfile_Facade& o);
    void operator= (const Ofstream_Binfile_Facade& o);
  public:
    //! Default Constructor
    Ofstream_Binfile_Facade ( );
    //! Constructor from filename and stream mode
    explicit Ofstream_Binfile_Facade ( const char * filename,
      std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary);
    //! Open state
    bool is_open() {return _str->is_open();}
    //! Method to open corresponding file
    void open ( const char * filename,
      std::ios_base::openmode mode = std::ios_base::out | std::ios_base::binary )
    {_str->open(filename,mode);}
    //! Method to close corresponding file
    void close()
    {_str->close();}
    //! Output multiple characters
    Ofstream_Binfile_Facade& write (const char* c, std::streamsize n)
    {_str->write(c,n); return *this;}
    //! Output single char
    Ofstream_Binfile_Facade& put (const char c)
    {_str->put(c); return *this;};
    //! Get position
    std::streampos tellp()
    {return _str->tellp();}
    //! Set position
    Ofstream_Binfile_Facade& seekp (std::streampos pos)
    {_str->seekp(pos); return *this;}
    //! Set relative position
    Ofstream_Binfile_Facade& seekp (std::streamoff pos, std::ios_base::seekdir way)
    {_str->seekp(pos,way); return *this;}
    //! Flushes stream buffer
    Ofstream_Binfile_Facade& flush()
    {_str->flush(); return *this;}

    //! This method returns true is stream state is good
    bool good() const {return _str->good();}
    //! This method returns true if eof is reached
    bool eof() const {return _str->eof();}
    //! This method returns true if either failbit or badbit is set
    bool fail() const {return _str->fail();}
    //! This method returns true if badbit is set
    bool bad() const {return _str->bad();}

    //! Unary not operator to check the stream state
    bool operator!() const {return _str->fail();}
    //! Conversion to bool to validate stream state
    operator bool() const {return _str->good();}

    //! Method to read stream state flags
    std::ios_base::iostate rdstate() const {return _str->rdstate();}
    //! Method to set the stream state (combines already set flags with flags provide by user)
    void setstate (std::ios_base::iostate state) {_str->setstate(state);}
    //! Method to set stream state (overwrites stream state flags)
    void clear (std::ios_base::iostate state = std::ios_base::goodbit) {_str->clear(state);}
    //! Method to get the exceptions mask
    std::ios_base::iostate exceptions() const {return _str->exceptions();}
    //! Method to set the exceptions mask
    void exceptions (std::ios_base::iostate except) {_str->exceptions(except);}

    //! Destructor
    virtual ~Ofstream_Binfile_Facade();
  protected:
    //! Access to internal stream for derived classes
    std::ofstream* stream() {return _str;}
  };

/*!
  \brief Ifstream Interface Facade for Binary Streams.
  \ingroup itfile

  This class implements std::ofstream facade to make ITPP binary file streams exportable from dll.
  This facade implements basic functionality only. It does not provide an access
  to the following stream facilities (all of them are useless for binary streams)
  1. locale(imbue) It does not make sence to change locale settings for binary streams.
  Changes in formatting or char conversion can result in compatibility problems
  with resulting binary files
  2. stream buffer (rdbuf). DLL and application can use different versions of runtime ,
  so it would be dangerous to use buffer created in DLL in application context
  3. stream extraction operators. It is assumed that stream extraction is defined in binary stream
  classes derived from this class
  4. formatting interface (copyfmt, fill, narrow, widen). This is not relevant to binary streams
  5. ios_base-related stuff. These things are excluded since they provide unnecessarily formatting
  facilities.
*/
  class ITPP_EXPORT Ifstream_Binfile_Facade
  {
    std::ifstream* _str;
    //following makes class non-copiable
    Ifstream_Binfile_Facade(const Ifstream_Binfile_Facade& o);
    void operator= (const Ifstream_Binfile_Facade& o);
  public:
    //! Default Constructor
    Ifstream_Binfile_Facade ( );
    //! Constructor from filename and stream mode
    explicit Ifstream_Binfile_Facade ( const char * filename,
      std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary);
    //! Open state
    bool is_open()
    {return _str->is_open();}
    //! Method to open corresponding file
    void open ( const char * filename,
      std::ios_base::openmode mode = std::ios_base::in | std::ios_base::binary)
    {_str->open(filename,mode);}
    //! Method to close corresponding file
    void close() {_str->close();}
    //! Last extracted chars count
    std::streamsize gcount() const {return _str->gcount();}
    //! Get single char
    int get() {return _str->get();}
    //! Get single char
    Ifstream_Binfile_Facade& get(char& c) {_str->get(c); return *this;}
    //! Get multiple chars to c-string and add trailing 0
    Ifstream_Binfile_Facade& get (char* s, std::streamsize n)
    {_str->get(s,n); return *this;}
    //! Get multiple chars to c-string without trailing 0
    Ifstream_Binfile_Facade& get (char* s, std::streamsize n, char delim)
    {_str->get(s,n,delim); return *this;}
    //! Get multiple chars to c-string without trailing 0
    Ifstream_Binfile_Facade& getline (char* s, std::streamsize n )
    {_str->getline(s,n); return *this;}
    Ifstream_Binfile_Facade& getline (char* s, std::streamsize n, char delim )
    {_str->getline(s,n,delim); return *this;}
    //! Extract and ignore chars
    Ifstream_Binfile_Facade& ignore (std::streamsize n = 1, int delim = EOF)
    {_str->ignore(n,delim); return *this;}
    //! Peak single char from the top of the buffer
    int peek() {return _str->peek();}
    //! Read n chars from stream
    Ifstream_Binfile_Facade& read (char* s, std::streamsize n)
    {_str->read(s,n); return *this;}
    //! Read up to n available chars from stream
    std::streamsize readsome (char* s, std::streamsize n)
    {return _str->readsome(s,n);}
    //! This method attempts to put back single char
    Ifstream_Binfile_Facade& putback (char c)
    {_str->putback(c); return *this;}
    //! Unget last extracted char
    Ifstream_Binfile_Facade& unget()  {_str->unget(); return *this;}
    //! Get position
    std::streampos tellg() {return _str->tellg();}
    //! Set position
    Ifstream_Binfile_Facade& seekg (std::streampos pos)
    {_str->seekg(pos); return *this;}
    //! Set relative position
    Ifstream_Binfile_Facade& seekg (std::streamoff pos, std::ios_base::seekdir way)
    {_str->seekg(pos,way); return *this;}

    //! This method returns true is stream state is good
    bool good() const {return _str->good();}
    //! This method returns true if eof is reached
    bool eof() const {return _str->eof();}
    //! This method returns true if either failbit or badbit is set
    bool fail() const {return _str->fail();}
    //! This method returns true if badbit is set
    bool bad() const {return _str->bad();}

    //! Unary not operator to check the stream state
    bool operator!() const {return _str->fail();}
    //! Conversion to bool to validate stream state
    operator bool() const {return _str->good();}

    //! Method to read stream state flags
    std::ios_base::iostate rdstate() const {return _str->rdstate();}
    //! Method to set the stream state (combines already set flags with flags provide by user)
    void setstate (std::ios_base::iostate state) {_str->setstate(state);}
    //! Method to set stream state (overwrites stream state flags)
    void clear (std::ios_base::iostate state = std::ios_base::goodbit) {_str->clear(state);}
    //! Method to get the exceptions mask
    std::ios_base::iostate exceptions() const {return _str->exceptions();}
    //! Method to set the exceptions mask
    void exceptions (std::ios_base::iostate except) {_str->exceptions(except);}

    //! Destructor
    virtual ~Ifstream_Binfile_Facade();
  protected:
    //! Access to internal stream for derived classes
    std::ifstream* stream() {return _str;}
  };

/*!
  \brief Fstream Interface Facade for Binary Streams.
  \ingroup itfile

  This class implements std::fstream facade to make ITPP binary file streams exportable from dll.
  This facade implements basic functionality only. It does not provide an access
  to the following stream facilities (all of them are useless for binary streams)
  1. locale(imbue) It does not make sence to change locale settings for binary streams.
  Changes in formatting or char conversion can result in compatibility problems
  with resulting binary files
  2. stream buffer (rdbuf). DLL and application can use different versions of runtime ,
  so it would be dangerous to use buffer created in DLL in application context
  3. stream extraction operators. It is assumed that stream extraction is defined in binary stream
  classes derived from this class
  4. stream insertion operators. It is assumed that stream insertion is defined in binary stream
  classes derived from this class
  5. formatting interface (copyfmt, fill, narrow, widen). This is not relevant to binary streams
  6. ios_base-related stuff. These things are excluded since they provide unnecessarily formatting facilities.
*/
  class ITPP_EXPORT Fstream_Binfile_Facade
  {
    std::fstream* _str;
    //following makes class non-copiable
    Fstream_Binfile_Facade(const Fstream_Binfile_Facade& o);
    void operator= (const Fstream_Binfile_Facade& o);
  public:
    //! Default Constructor
    Fstream_Binfile_Facade ( );
    //! Constructor from filename and stream mode
    explicit Fstream_Binfile_Facade ( const char * filename,
      std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out | std::ios_base::binary);
    //! Open state
    bool is_open() {return _str->is_open();}
    //! Method to open corresponding file
    void open ( const char * filename,
      std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out | std::ios_base::binary)
    {_str->open(filename,mode);}
    //! Method to close corresponding file
    void close() {_str->close();}

    //! Output multiple characters
    Fstream_Binfile_Facade& write (const char* c, std::streamsize n)
    {_str->write(c,n); return *this;};
    //! Output single char
    Fstream_Binfile_Facade& put (const char c)
    {_str->put(c); return *this;};
    //! Get position
    std::streampos tellp() {return _str->tellp();}
    //! Set position
    Fstream_Binfile_Facade& seekp (std::streampos pos)
    {_str->seekp(pos); return *this;}
    //! Set relative position
    Fstream_Binfile_Facade& seekp (std::streamoff pos, std::ios_base::seekdir way)
    {_str->seekp(pos,way); return *this;}
    //! Flushes stream buffer
    Fstream_Binfile_Facade& flush() {_str->flush(); return *this;}
    //! Last extracted chars count
    std::streamsize gcount() const {return _str->gcount();}
    //! Get single char
    int get() {return _str->get();}
    //! Get single char
    Fstream_Binfile_Facade& get(char& c){_str->get(c); return *this;}
    //! Get multiple chars to c-string and add trailing 0
    Fstream_Binfile_Facade& get(char* s, std::streamsize n)
    {_str->get(s,n); return *this;}
    //! Get multiple chars to c-string without trailing 0
    Fstream_Binfile_Facade& get(char* s, std::streamsize n, char delim)
    {_str->get(s,n,delim); return *this;}
    //! Get multiple chars to c-string without trailing 0
    Fstream_Binfile_Facade& getline(char* s, std::streamsize n)
    {_str->getline(s,n); return *this;}
    Fstream_Binfile_Facade& getline(char* s, std::streamsize n, char delim)
    {_str->getline(s,n,delim); return *this;}
    //! Extract and ignore chars
    Fstream_Binfile_Facade& ignore (std::streamsize n = 1, int delim = EOF)
    {_str->ignore(n,delim); return *this;}
    //! Peak single char from the top of the buffer
    int peek() {return _str->peek();}
    //! Read n chars from stream
    Fstream_Binfile_Facade& read (char* s, std::streamsize n)
    {_str->read(s,n); return *this;}
    //! Read up to n available chars from stream
    std::streamsize readsome (char* s, std::streamsize n)
    {return _str->readsome(s,n);}
    //! This method attempts to put back single char
    Fstream_Binfile_Facade& putback (char c)
    {_str->putback(c); return *this;}
    //! Unget last extracted char
    Fstream_Binfile_Facade& unget()
    {_str->unget(); return *this;}
    //! Get position
    std::streampos tellg() {return _str->tellg();}
    //! Set position
    Fstream_Binfile_Facade& seekg (std::streampos pos)
    {_str->seekg(pos); return *this;}
    //! Set relative position
    Fstream_Binfile_Facade& seekg (std::streamoff pos, std::ios_base::seekdir way)
    {_str->seekg(pos,way); return *this;}

    //! This method returns true is stream state is good
    bool good() const {return _str->good();}
    //! This method returns true if eof is reached
    bool eof() const {return _str->eof();}
    //! This method returns true if either failbit or badbit is set
    bool fail() const {return _str->fail();}
    //! This method returns true if badbit is set
    bool bad() const {return _str->bad();}

    //! Unary not operator to check the stream state
    bool operator!() const {return _str->fail();}
    //! Conversion to bool to validate stream state
    operator bool() const {return _str->good();}

    //! Method to read stream state flags
    std::ios_base::iostate rdstate() const {return _str->rdstate();}
    //! Method to set the stream state (combines already set flags with flags provide by user)
    void setstate (std::ios_base::iostate state) {_str->setstate(state);}
    //! Method to set stream state (overwrites stream state flags)
    void clear (std::ios_base::iostate state = std::ios_base::goodbit)
    {_str->clear(state);}
    //! Method to get the exceptions mask
    std::ios_base::iostate exceptions() const {return _str->exceptions();}
    //! Method to set the exceptions mask
    void exceptions (std::ios_base::iostate except) {_str->exceptions(except);}

    //! Destructor
    virtual ~Fstream_Binfile_Facade();
  protected:
    //! Access to internal stream for derived classes
    std::fstream* stream() {return _str;}
  };


}

/*!
  \brief Binary Outfile Class
  \ingroup itfile
*/
class ITPP_EXPORT bofstream : public bfstream_base, public binfile_details::Ofstream_Binfile_Facade
{
public:
  /*!
    \brief Class constructor that opens a file and sets the endianity

    \param name The name of the file to open
    \param e Defines the endianity of the class. Possible values are
    \c l_endian for "Little Endian" or \c b_endian for "Big Endian". The
    default value is \c b_endian.
    Set \c truncate to true to discard file contents.
  */
  bofstream(const std::string& name, endian e = b_endian);

  //! Class Constructor
  bofstream();

  //! Class Destructor
  ~bofstream() { }

  /*!
    \brief Open a file for writing and set the endianity

    \param name The name of the file to open
    \param e Defines the endianity of the class (default value is
    \c b_endian )
    Set \c trunc to true to discard file contents.
  */
  void open(const std::string& name, bool trunc = false, endian e = b_endian);

  //! Writes a signed char variable to the binary output file
  bofstream& operator<<(char a);
  //! Writes a 8-bit signed integer variable to the binary file
  bofstream& operator<<(int8_t a);
  //! Writes a 8-bit unsigned integer variable to the binary file
  bofstream& operator<<(uint8_t a);
  //! Writes a 16-bit signed integer variable to the binary output file
  bofstream& operator<<(int16_t a);
  //! Writes a 16-bit unsigned integer variable to the binary output file
  bofstream& operator<<(uint16_t a);
  //! Writes a 32-bit signed integer variable to the binary output file
  bofstream& operator<<(int32_t a);
  //! Writes a 32-bit unsigned integer variable to the binary output file
  bofstream& operator<<(uint32_t a);
  //! Writes a 64-bit signed integer variable to the binary output file
  bofstream& operator<<(int64_t a);
  //! Writes a 64-bit unsigned ingeger variable to the binary output file
  bofstream& operator<<(uint64_t a);
  //! Writes a float variable to the binary output file
  bofstream& operator<<(float a);
  //! Writes a double variable to the binary output file
  bofstream& operator<<(double a);
  //! Writes a binary variable to the binary output file
  bofstream& operator<<(bin a);
  //! Writes a char* string to the binary output file
  bofstream& operator<<(const char* a);
  //! Writes a string variable to the binary output file
  bofstream& operator<<(const std::string& a);
};

/*!
  \brief Binary Infile Class
  \ingroup itfile
*/
class ITPP_EXPORT bifstream : public bfstream_base, public binfile_details::Ifstream_Binfile_Facade
{
public:
  /*!
    \brief Class constructor that opens a file and sets the endianity

    \param name The name of the file to open
    \param e Defines the endianity of the class. Possible values are
    \c l_endian for "Little Endian" or \c b_endian for "Big Endian". The
    default value is \c b_endian.
  */
  bifstream(const std::string& name, endian e = b_endian);

  //! Class Constructor
  bifstream();

  //! Class Destructor
  ~bifstream() { }

  /*!
    \brief Open a file for reading and set the endianity

    \param name The name of the file to open
    \param e Defines the endianity of the class (default value is
    \c b_endian )
  */
  void open(const std::string& name, endian e = b_endian);

  //! Returns the length in bytes of the file
  int length();

  //! Reads a signed char variable from the binary input file
  bifstream& operator>>(char& a);
  //! Reads a 8-bit signed integer variable from the binary input file
  bifstream& operator>>(int8_t& a);
  //! Reads a 8-bit unsigned integer variable from the binary input file
  bifstream& operator>>(uint8_t& a);
  //! Reads a 16-bit signed integer variable from the binary input file
  bifstream& operator>>(int16_t& a);
  //! Reads a 16-bit unsigned integer variable from the binary input file
  bifstream& operator>>(uint16_t& a);
  //! Reads a 32-bit signed integer variable from the binary input file
  bifstream& operator>>(int32_t& a);
  //! Reads a 32-bit unsigned integer variable from the binary input file
  bifstream& operator>>(uint32_t& a);
  //! Reads a 64-bit signed integer variable from the binary input file
  bifstream& operator>>(int64_t& a);
  //! Reads a 64-bit unsigned ingeger variable from the binary input file
  bifstream& operator>>(uint64_t& a);
  //! Reads a float variable from the binary input file
  bifstream& operator>>(float& a);
  //! Reads a double variable from the binary input file
  bifstream& operator>>(double& a);
  //! Reads a binary variable from the binary input file
  bifstream& operator>>(bin& a);
  //! Reads a char* string from the binary input file
  bifstream& operator>>(char* a);
  //! Reads a string variable from the binary input file
  bifstream& operator>>(std::string& a);
};

/*!
  \brief Binary in/out-file Class
  \ingroup itfile
*/
class ITPP_EXPORT bfstream : public bfstream_base, public binfile_details::Fstream_Binfile_Facade
{
public:
  /*!
    \brief Class constructor that opens a file and sets the endianity

    \param name The name of the file to open
    \param e Defines the endianity of the class. Possible values are
    \c l_endian for "Little Endian" or \c b_endian for "Big Endian".
    The default value is \c b_endian.
  */
  bfstream(const std::string& name, endian e = b_endian);

  //! Class Constructor
  bfstream();

  //! Class Destructor
  ~bfstream() { }

  /*!
    \brief Open a file for reading and writing and set the endianity

    \param name The name of the file to open
    \param trunc Rewrite the file if it exists (default value is \c false)
    \param e Defines the endianity of the class (default value is
    \c b_endian )
  */
  void open(const std::string& name, bool trunc = false, endian e = b_endian);

  /*!
    \brief Open a file for reading only and set the endianity

    \param name The name of the file to open
    \param e Defines the endianity of the class (default value is
    \c b_endian )
  */
  void open_readonly(const std::string& name, endian e = b_endian);

  //! Returns the length in bytes of the file
  int length();

  //! Writes an char variable to the binary file
  bfstream& operator<<(char a);
  //! Writes a 8-bit signed integer variable to the binary file
  bfstream& operator<<(int8_t a);
  //! Writes a 8-bit unsigned integer variable to the binary file
  bfstream& operator<<(uint8_t a);
  //! Writes a 16-bit signed integer variable to the binary file
  bfstream& operator<<(int16_t a);
  //! Writes a 16-bit unsigned integer variable to the binary file
  bfstream& operator<<(uint16_t a);
  //! Writes a 32-bit signed integer variable to the binary file
  bfstream& operator<<(int32_t a);
  //! Writes a 32-bit unsigned integer variable to the binary file
  bfstream& operator<<(uint32_t a);
  //! Writes a 64-bit signed integer variable to the binary file
  bfstream& operator<<(int64_t a);
  //! Writes a 64-bit unsigned ingeger variable to the binary file
  bfstream& operator<<(uint64_t a);
  //! Writes a float variable to the binary file
  bfstream& operator<<(float a);
  //! Writes a double variable to the binary file
  bfstream& operator<<(double a);
  //! Writes a binary variable to the binary file
  bfstream& operator<<(bin a);
  //! Writes a char* string to the binary file
  bfstream& operator<<(const char* a);
  //! Writes a string variable to the binary file
  bfstream& operator<<(const std::string& a);

  //! Reads a char variable from the binary file
  bfstream& operator>>(char& a);
  //! Reads a 8-bit signed integer variable from the binary file
  bfstream& operator>>(int8_t& a);
  //! Reads a 8-bit unsigned integer variable from the binary file
  bfstream& operator>>(uint8_t & a);
  //! Reads a 16-bit signed integer variable from the binary file
  bfstream& operator>>(int16_t& a);
  //! Reads a 16-bit unsigned integer variable from the binary file
  bfstream& operator>>(uint16_t& a);
  //! Reads a 32-bit signed integer variable from the binary file
  bfstream& operator>>(int32_t& a);
  //! Reads a 32-bit unsigned integer variable from the binary file
  bfstream& operator>>(uint32_t& a);
  //! Reads a 64-bit signed integer variable from the binary file
  bfstream& operator>>(int64_t& a);
  //! Reads a 64-bit unsigned ingeger variable from the binary file
  bfstream& operator>>(uint64_t& a);
  //! Reads a float variable from the binary file
  bfstream& operator>>(float& a);
  //! Reads a double variable from the binary file
  bfstream& operator>>(double& a);
  //! Reads a binary variable from the binary file
  bfstream& operator>>(bin& a);
  //! Reads a char* string from the binary file
  bfstream& operator>>(char* a);
  //! Reads a string variable from the binary file
  bfstream& operator>>(std::string& a);
};

} //namespace itpp


#endif // #ifndef BINFILE_H
