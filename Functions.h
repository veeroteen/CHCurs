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
      if(compiled)
      {
         for (auto &a : variables)
         {
            a.second = args.at(a.first);
         }
      }
      else
      {
         for(auto &a : args)
         {
            variables[a.first] = a.second;
         }
         exprtk::symbol_table<T> table;
         for (auto &a : variables)
         {
            table.add_variable(a.first, a.second);
         }
         table.add_constants();
         expression.register_symbol_table(table);
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
      evaluate(map);
   }

   std::string& getStringFun()
   {
      return _fun;
   }
};