@rem File:   setenv.bat.py
@rem Brief:  Setting the environment for using Microsoft Visual Studio command line tools.
@rem Author: Bogdan Cristea
@rem         Initial version provided by B.Komazec
@rem
@rem Usage: setenv.bat x64
@rem
@rem This batch file calls MS Visual Studio provided batch files (vcvars32.bat or vcvars64.bat)
@rem in order to setup the environment variables needed for Microsoft Visual Studio
@rem command line tools. Note that the paths to vcvars32.bat or vcvars64.bat are hard
@rem coded and needed to be provided if needed. By default these paths correspond to
@rem Microsoft Visual Studio 2010 installation.
@rem
@rem -------------------------------------------------------------------------
@rem
@rem Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
@rem
@rem This file is part of IT++ - a C++ library of mathematical, signal
@rem processing, speech processing, and communications classes and functions.
@rem
@rem IT++ is free software: you can redistribute it and/or modify it under the
@rem terms of the GNU General Public License as published by the Free Software
@rem Foundation, either version 3 of the License, or (at your option) any
@rem later version.
@rem
@rem IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
@rem WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
@rem FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
@rem details.
@rem
@rem You should have received a copy of the GNU General Public License along
@rem with IT++.  If not, see <http://www.gnu.org/licenses/>.
@rem
@rem -------------------------------------------------------------------------

@echo off

@if "%1"=="x86" goto set_x86
@if "%1"=="x64" goto set_x64
@if "%1"=="" goto error

:set_x86

call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\vcvars32.bat"

goto test_bin_locations

:set_x64

call "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\amd64\vcvars64.bat"

goto test_bin_locations

:test_bin_locations
@echo on
where nmake
where cl.exe
where link.exe
@echo off
goto:eof

:error
@echo Usage: setenv.bat [x86^|x64]

goto:eof
