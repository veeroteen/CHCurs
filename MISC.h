//MISC.h
#pragma once
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

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
   T res = poly[0];
   
   for (size_t i = 0; i < (poly.size() - 1) / 3; i++)
   {
      res += poly[1 + i * 3] * pow(x, i + 1);
      res += poly[2 + i * 3] * pow(y, i + 1);
      res += poly[3 + i * 3] * pow(z, i + 1);

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

template<Field T>
void polyToStr(std::string &res,std::vector<T> &poly)
{
   res.clear();

   res += std::to_string(poly[0]);
   res += std::string("+" + std::to_string(poly[1]) + "x");
   res += std::string("+" + std::to_string(poly[2]) + "y");
   res += std::string("+" + std::to_string(poly[3]) + "z");
}