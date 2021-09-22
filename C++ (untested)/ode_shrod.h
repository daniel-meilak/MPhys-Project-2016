#pragma once

#include<array>
#include<complex>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<unordered_map>


#include"complex_2D.h"

enum shape {triangle, parabolic_traingle, rectangle, double_rectangle, trapezoid};

// globals (fallback from Fortran)
extern double energy, height, para_cons;
extern std::array<double,4> param;
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

   const std::string input_file = "input2,dat";

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
   
   input >> height >> param[0] >> param[1];
   if (potential_shape == shape::double_rectangle || potential_shape == shape::trapezoid){
      input >> param[2] >> param[3];
   }

   // calculated for the parabolic triangle case (ignore otherwise)
   para_cons = height/(param[0]*param[0] - 2.0*param[0]*param[1] + param[1]*param[1]);

   input >> energy;
   input.close();
}

//======================================================================
// Defines the ode under study. This is the Schrödinger Eqn             
//----------------------------------------------------------------------
// Units:                                                               
// The TISE has been written using atomic units, therefore              
// x is in Hartrees                                                     
// E is in Bohrs                                                        
//======================================================================
complex_2D<double> schrod(
   const double& x,
   const complex_2D<double>& psi){

   // output
   complex_2D<double> output;

   output.first = psi.second;
   output.second = 2.0*(potential(x)-energy)*psi.first;

   return output;
}

//======================================================================
// This function take the constants read in by read_pot and then        
// calculates a value for the potential at a particular value of x      
//----------------------------------------------------------------------
// x is in units of Bohr (atomic units)                                 
//======================================================================
double potential(const double& x){

   double output;

   switch (potential_shape){
      //   ___/\___
      case triangle:  
         if (x >= param[0] && x <= param[1]){
            output = height*x/(param[0]-param[1]) + height*param[1]/(param[1]-param[0]);
         }
         else { output = 0.0; }
         break;

      //   ___()___
      case parabolic_traingle:
         if (x >= param[0] && x <= param[1]){
            output = para_cons*(x*x - param[1]*(2.0*x + param[1]));
         }
         else { output = 0.0; }
         break;
      
      //   ___|``|___
      case rectangle:
         if (x >= param[0] && x <= param[1]){ output = height; }
         else { output = 0.0; }
         break;

      //   ___|``|___|``|___
      case double_rectangle:
         if ((x >= param[0] && x <= param[1]) || (x >= param[2] && x <= param[3])){ output = height; }
         else { output = 0.0; }
         break;

      //   ___/```\___
      case trapezoid:
         if      (x >= param[0] && x <= param[1]){ output = height*(x-param[0])/(param[1]-param[0]); }
         else if (x >= param[1] && x <= param[2]){ output = height; }
         else if (x >= param[2] && x <= param[3]){ output = height + height*(x-param[2])/(param[2]-param[3]); }
         else    { output = 0.0; }         
   }

   return output;
}
