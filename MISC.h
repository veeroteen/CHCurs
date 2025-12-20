//MISC.h
#pragma once
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>

void genGil(size_t n,std::string &fileName)
{
   std::ofstream file(fileName);
   file << n << std::endl;
   for(size_t i = 0; i < n; i++)
   {
      for(size_t j = 0; j < n;j++)
      {
         file << std::setprecision(16) << 1.0 / double((i + j + 1)) << " ";
      }
      file << i+1;
      if (i != n - 1) 
      {
         file << std::endl;
      }
   }

}

template <typename T>
concept Field = requires(T a, T b, std::istream & is, std::ostream & os) {
   { a + b } -> std::same_as<T>;
   { a - b } -> std::same_as<T>;
   { a *b } -> std::same_as<T>;
   { a / b } -> std::same_as<T>;

   { a == b } -> std::convertible_to<bool>;
   { a != b } -> std::convertible_to<bool>;
   { a > b }-> std::convertible_to<bool>;
   { a < b }-> std::convertible_to<bool>;
   { a = b } -> std::same_as<T &>;
   { std::abs(a) };
   { os << a } -> std::same_as<std::ostream &>;
   { is >> a } -> std::same_as<std::istream &>;
};

template <Field T>
T calcPoly(const std::vector<double> &poly, T x, T y, T z)
{
   //T res = poly[0];
   
   //for (size_t i = 0; i < (poly.size() - 1) / 3; i++)
   {
      if (poly.size() != 3)
      {
         T res = poly[0];//delete
         res += poly[1 + i * 3] * pow(x, i + 1);
         res += poly[2 + i * 3] * pow(y, i + 1);
         res += poly[3 + i * 3] * pow(z, i + 1);
      }
      else
      {
         res += sin(x + y + z);
      
      
      }
   }
   return res;
}



template<Field T>
T pow(T a, unsigned b)
{
   T res = (b == 0 ? 1 : a);

   for (unsigned i = 1; i < b; i++)
   {
      res *= a;
   }
   return res;
}

