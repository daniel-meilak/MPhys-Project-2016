#pragma once

#include<algorithm>
#include<array>
#include<cmath>
#include<complex>
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<ranges>
#include<string>
#include<string_view>

#include"complex_2D.h"
#include"ode_shrod.h"

// globals (fallback from Fortran)
extern double energy, height;

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
   const std::string input_file  = "input1.dat";
   const std::string output_file = "output.dat";
   
   // parameters
   int n_iter;
   double x, height, start, finish;
   std::complex<double> wave_n;
   complex_2D<double> psi, k1, k2, k3, k4; 

   // read input parameters
   std::ifstream input(input_file);
   if (!input.is_open()){
      std::cerr << "Error in opening " << input_file << std::endl;
      std::exit(EXIT_FAILURE);
   }
   input >> start >> finish;
   input >> n_iter;
   input.close();

   // calculate starting psi and psi'
   wave_n     = std::complex<double>(0.0,std::sqrt(2.0*energy));
   psi.first  = std::exp(wave_n*finish);
   psi.second = psi.first*wave_n;
   height     = finish-start/n_iter;
   x          = finish;

   // start writing output data (overwrite old output)
   std::ofstream output(output_file, std::ios::trunc);
   if (!output.is_open()){
      std::cerr << "Error in opening " << output_file << std::endl;
      std::exit(EXIT_FAILURE);
   }
   output << x << "\t" << psi.first.real() << "\t" << psi.first.imag() << "\t"
          << psi.second.real() << "\t" << psi.second.imag() << "\n";

   for (int i=1; i<=n_iter; i++){

      // 4 function calls, each depends on the preceding one
      k1 = height*schrod(x           , psi       );
      k2 = height*schrod(x-height/2.0, psi-k1/2.0);
      k3 = height*schrod(x-height/2.0, psi-k2/2.0);
      k4 = height*schrod(x-height    , psi-k3    );

      psi -= (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
      x    = x - height;

      output << x << "\t" << psi.first.real() << "\t" << psi.first.imag() << "\t"
          << psi.second.real() << "\t" << psi.second.imag() << "\n";
   }

   output.close();
}