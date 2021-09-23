#include"runge_kutta.h"
#include"ode_shrod.h"

// As this program comes from rewriting Fortran, where global variables were used,
// the same method has been used here, so that the code does not need to be 
// changed extensively
long double energy, height, para_cons;
std::array<long double,4> param;
enum shape potential_shape;

int main(){

   read_param();
   runge_kutta();
   find_coeff();

   return 0;
}