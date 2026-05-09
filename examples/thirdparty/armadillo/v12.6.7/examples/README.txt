How to compile example1.cpp


** Linux and macOS **

If you have installed Armadillo via the CMake installer:
  g++ example1.cpp -o example1 -std=c++11 -O2 -larmadillo
  
Otherwise, if you want to use Armadillo without installation:
  g++ example1.cpp -o example1 -std=c++11 -O2 -I /home/user/armadillo-12.2.0/include -DARMA_DONT_USE_WRAPPER -lopenblas
  
The above command assumes that the armadillo archive was unpacked into /home/user/
The command needs to be adjusted if the archive was unpacked into a different directory,
and/or for each specific version of Armadillo (ie. "12.2.0" needs to be changed).
  
If you don't have OpenBLAS, on Linux change -lopenblas to -lblas -llapack
and on macOS change -lopenblas to -framework Accelerate


** Windows **
  
Open "example1_win64.sln" or "example1_win64.vcxproj" with Visual Studio.
The example1_win64 project needs to be compiled as a 64 bit program.
Make sure the active solution platform is set to x64, instead of win32.
