#pragma once

#include<array>
#include<complex>

using namespace std::complex_literals;


//===========================================================================
// While fortran has built in vector,scalar maths, it must be done manually
// in C++. Simplified to work on size 2 arrays.
//===========================================================================
template<typename T>
struct complex_2D{

   std::complex<T> first;
   std::complex<T> second;

   // spaceship operator defines "==","!=",">","<",">=","<="
   auto operator<=>(const complex_2D&) const = default;

   // constructors
   complex_2D(): first(0.+0i), second(0.+0i){};
   template<typename S1, typename S2>
   complex_2D(const S1& first, const S2& second): first(first), second(second){};

   // += and -= operators
   complex_2D& operator+=(const complex_2D& p){
      this->fisrt  += p.first;
      this->second += p.second;
   }
   complex_2D& operator-=(const complex_2D& p){
      this->first  -= p.first;
      this->second -= p.second;
   }

};

// add 2D complex vectors
template<typename T>
complex_2D<T> operator+(const complex_2D<T>& lhs, const complex_2D<T>& rhs){
   return {lhs.first + rhs.first, lhs.second + rhs.second};
}

// subtract 2D complex vectors
template<typename T>
complex_2D<T> operator-(const complex_2D<T>& lhs, const complex_2D<T>& rhs){
   return {lhs.first - rhs.first, lhs.second - rhs.second};
}

// multiply vector by scalar
template<typename T, typename S>
complex_2D<T> operator*(const complex_2D<T>& complex, const S& scalar){
   return {complex.first*scalar, complex.second*scalar};
}
template<typename T, typename S>
complex_2D<T> operator*(const S& scalar, const complex_2D<T>& complex){
   return complex * scalar;
}

// divide vector by scalar
template<typename T, typename S>
complex_2D<T> operator/(const complex_2D<T>& complex, const S& scalar){
   return {complex.first/scalar, complex.second/scalar};
}