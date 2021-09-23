#pragma once

#include<algorithm>
#include<array>
#include<cmath>
#include<complex>
#include<cstdlib>
#include<iomanip>
#include<iostream>
#include<fstream>
#include<ranges>
#include<string>
#include<string_view>

#include"complex_2D.h"
#include"ode_shrod.h"

// globals (fallback from Fortran)
extern long double energy, height;

//===========================================================================
// This is an algorithm for the fourth order Runge Kutta method. It is a 4th 
// order ODE solver. However, this method is written so that it starts on    
// b and moves backwards, i.e., opposite to the normal direction             
//---------------------------------------------------------------------------
// Notes: 
//
//  The input file must exist when the program is run.                                                  
//===========================================================================
void runge_kutta(){

   // future note: constexpr when g++ implements C++20 fully
   const std::string input_file  = "../input1.dat";
   const std::string output_file = "output.dat";
   
   // parameters
   int n_iter;
   long double x, step_size, start, finish;
   std::complex<long double> wave_n;
   complex_2D<long double> psi, k1, k2, k3, k4; 

   // read input parameters
   std::ifstream input(input_file);
   if (!input.is_open()){
      std::cerr << "Error in opening " << input_file << std::endl;
      std::exit(EXIT_FAILURE);
   }
   input >> start >> finish;
   input.ignore(1000,'\n');
   input >> n_iter;
   input.close();

   // calculate starting psi and psi'
   wave_n     = std::complex<long double>(0.0l,std::sqrt(2.0l*energy));
   psi.first  = std::exp(wave_n*finish);
   psi.second = psi.first*wave_n;
   step_size  = (finish-start)/n_iter;
   x          = finish;

   // start writing output data (overwrite old output)
   std::ofstream output(output_file, std::ios::trunc);
   if (!output.is_open()){
      std::cerr << "Error in opening " << output_file << std::endl;
      std::exit(EXIT_FAILURE);
   }
   output << std::setprecision(15) << x << "\t" << psi.first.real() << "\t" << psi.first.imag()
          << "\t" << psi.second.real() << "\t" << psi.second.imag() << "\n";

   for (int i=1; i<=n_iter; i++){

      // 4 function calls, each depends on the preceding one
      k1 = step_size*schrod(x           , psi       );
      k2 = step_size*schrod(x-height/2.0l, psi-k1/2.0l);
      k3 = step_size*schrod(x-height/2.0l, psi-k2/2.0l);
      k4 = step_size*schrod(x-height    , psi-k3    );

      psi -= (k1 + 2.0l*k2 + 2.0l*k3 + k4)/6.0l;
      x   -= step_size;

      output << std::setprecision(15) << x << "\t" << psi.first.real() << "\t" << psi.first.imag()
             << "\t" << psi.second.real() << "\t" << psi.second.imag() << "\n";
   }

   output.close();
}