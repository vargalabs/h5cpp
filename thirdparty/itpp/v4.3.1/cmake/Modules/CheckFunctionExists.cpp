#ifdef CHECK_FUNCTION_EXISTS

//only for std:isfinite, std::isinf and std::isnan

#include <cmath>

#ifdef __CLASSIC_C__
int main(){
  int ac;
  char*av[];
#else
int main(int ac, char*av[]){
#endif
  double arg = 0.0;
  CHECK_FUNCTION_EXISTS(arg);
  if(ac > 1000)
    {
    return *av[0];
    }
  return 0;
}

#else  /* CHECK_FUNCTION_EXISTS */

#  error "CHECK_FUNCTION_EXISTS has to specify the function"

#endif /* CHECK_FUNCTION_EXISTS */
