#pragma once

#include<array>
#include<complex>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<numbers>
#include<sstream>
#include<unordered_map>


#include"complex_2D.h"

enum shape {triangle, parabolic_traingle, rectangle, double_rectangle, trapezoid};

// globals (fallback from Fortran)
extern long double energy, height, para_cons;
extern std::array<long double,4> param;
extern shape potential_shape;

// map to relate shape strings to enum, for switch case usage
static const std::unordered_map<std::string, shape> table = {
   {"triangle", shape::triangle},
   {"parabolic_triangle", shape::parabolic_traingle},
   {"rectangle", shape::rectangle},
   {"double_rectangle", shape::double_rectangle},
   {"trapezoid", shape::trapezoid}
};

//======================================================================
// Reads the form and dimentions of the given potential from a file,    
// input2.dat                                               
//----------------------------------------------------------------------
// This subroutine need to be called before using the Runge Kutta       
// method, if the ode is giong to be the shrodinger equation            
//======================================================================
void read_param(){

   const std::string input_file = "../input2.dat";

   std::ifstream input(input_file);
   if (!input.is_open()){
      std::cerr << "Error in opening " << input_file << std::endl;
      std::exit(EXIT_FAILURE);
   }

   // read potential shape
   std::string shape_string;
   input >> shape_string;
   auto it = table.find(shape_string);
   if (it != table.end()){ potential_shape = it->second; }
   else {
      std::cerr << "Invalid potential name\n";
      std::exit(EXIT_FAILURE); 
   }

   // move to next line
   input.ignore(1000,'\n');
   
   input >> height >> param[0] >> param[1];
   if (potential_shape == shape::double_rectangle || potential_shape == shape::trapezoid){
      input >> param[2] >> param[3];
   }

   // calculated for the parabolic triangle case (ignore otherwise)
   para_cons = height/(param[0]*param[0] - 2.0l*param[0]*param[1] + param[1]*param[1]);

   input.ignore(1000,'\n');
   input >> energy;
   input.close();
}

//======================================================================
// This function take the constants read in by read_pot and then        
// calculates a value for the potential at a particular value of x      
//----------------------------------------------------------------------
// x is in units of Bohr (atomic units)                                 
//======================================================================
long double potential(const double& x){

   long double output;

   switch (potential_shape){
      //   ___/\___
      case triangle:  
         if (x >= param[0] && x <= param[1]){
            output = height*x/(param[0]-param[1]) + height*param[1]/(param[1]-param[0]);
         }
         else { output = 0.0l; }
         break;

      //   ___()___
      case parabolic_traingle:
         if (x >= param[0] && x <= param[1]){
            output = para_cons*(x*x - param[1]*(2.0*x + param[1]));
         }
         else { output = 0.0l; }
         break;
      
      //   ___|``|___
      case rectangle:
         if (x >= param[0] && x <= param[1]){ output = height; }
         else { output = 0.0l; }
         break;

      //   ___|``|___|``|___
      case double_rectangle:
         if ((x >= param[0] && x <= param[1]) || (x >= param[2] && x <= param[3])){ output = height; }
         else { output = 0.0l; }
         break;

      //   ___/```\___
      case trapezoid:
         if      (x >= param[0] && x <= param[1]){ output = height*(x-param[0])/(param[1]-param[0]); }
         else if (x >= param[1] && x <= param[2]){ output = height; }
         else if (x >= param[2] && x <= param[3]){ output = height + height*(x-param[2])/(param[2]-param[3]); }
         else    { output = 0.0l; }         
   }

   return output;
}

//======================================================================
// Defines the ode under study. This is the Schrödinger Eqn             
//----------------------------------------------------------------------
// Units:                                                               
// The TISE has been written using atomic units, therefore              
// x is in Hartrees                                                     
// E is in Bohrs                                                        
//======================================================================
complex_2D<long double> schrod(
   const long double& x,
   const complex_2D<long double>& psi){

   // output
   complex_2D<long double> output;

   output.first = psi.second;
   output.second = 2.0l*(potential(x)-energy)*psi.first;

   return output;
}

//==========================================================================
// This subroutine goes through the newly created output.dat file and gets  
// the values of psi at x=0 and x=pi/2*k. Then, the values are used to      
// calculate the Reflection coefficient (and hence Transmission)               
//==========================================================================
void find_coeff(){

   const std::string input_file = "output.dat";

   std::ifstream input(input_file);
   if (!input.is_open()){
      std::cerr << "Error in opening " << input_file << std::endl;
      std::exit(EXIT_FAILURE);
   }

   long double x, Re1, Im1, Re2, Im2;
   complex_2D<long double> psi;

   long double k = std::sqrt(2.0l*energy);

   while (input >> x >> Re1 >> Im1 >> Re2 >> Im2){

      if (std::abs(x) > 1.0e-6l){
         psi.first = std::complex<long double>(Re1,Im1);
      }

      if (std::abs(x) > std::numbers::pi*0.5l/k){
         psi.second = std::complex<long double>(Re1,Im1);
      }
   }

   input.close();

   // calculate Reflection coeff
   std::complex<long double> At = (psi.first - std::complex<long double>(0.0l,1.0l)*psi.second)*0.5l;
   std::complex<long double> Ar = (std::complex<long double>(0.0l,1.0l)*psi.second + psi.first)*0.5l;
   long double refl = (std::abs(Ar)*std::abs(Ar))/(std::abs(At)*std::abs(At));

   std::cout << "Coeff: " << 1.0l-refl << std::endl;
}