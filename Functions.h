//Function.h
#pragma once
#include "MISC.h"
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <array>
#include "Calculator/exprtk/exprtk.hpp"

template <Field T>
class GlobalDict
{
private:
   inline static std::map<std::string, T> variables;
   inline static exprtk::symbol_table<T> table;
   inline static bool compiled = false;
public:
   static void setVal(const std::map<std::string, T> &args)
   {
      if (!compiled)
      {
         variables["x"] = 0;
         variables["y"] = 0;
         variables["z"] = 0;
         variables["u"] = 0;
         for (auto &a : variables)
         {
            table.add_variable(a.first, a.second);
         }
         table.add_constants();
         compiled = true;
      }
      for (auto &a : args)
      {
         variables[a.first] = a.second;
      }
   }
   static void initT(T &t)
   {
      table.add_variable("t", t);
   }
   static void regTable(exprtk::expression<T> &expr)
   {
      expr.register_symbol_table(table);
   }

};

template <Field T>
class Fun
{

public:
   virtual T evaluate(const std::map<std::string, T> &args) = 0;
};

template <Field T>
class StaticFun : public Fun<T>
{
   std::function<T(const std::map<std::string, T> &)> _fun;
public:
   StaticFun(std::function<T(const std::map<std::string, T> &)> fun)
   {
      _fun = fun;
   }
   T evaluate(const std::map<std::string, T> &args) override
   {
      return _fun(args);
   }

};

template <Field T>
class StringFun : public Fun<T>
{
   std::string _fun;
   bool compiled = false;
   exprtk::expression<T> expression;
   std::map<std::string, T> variables;
public:
   StringFun(const std::string &fun)
   {
      _fun = fun;
   }
   StringFun() = default;
   T evaluate(const std::map<std::string, T> &args) override
   {
      if (compiled)
      {
         GlobalDict<T>::setVal(args);
      }
      else
      {
         GlobalDict<T>::setVal(args);
         GlobalDict<T>::regTable(expression);
         exprtk::parser<T> parser;
         parser.compile(_fun, expression);
         compiled = true;
      }
      return expression.value();
   }
   T evaluate(const std::array<T, 3> &node)
   {
      std::map<std::string, T> map;
      map["x"] = node[0];
      map["y"] = node[1];
      map["z"] = node[2];
      return evaluate(map);
   }

   std::string &getStringFun()
   {
      return _fun;
   }
};


