//Gauss
#pragma once
#include <fstream>
#include <math.h>
#include <iostream>
#include <random>
#include <ctime>
#include "MISC.h"

template<Field T>
class GaussAlgoth
{
private:


   std::vector<std::vector<T>*> *matrix;
   unsigned findMax(unsigned j)
   {
      auto &matr = *matrix;
      std::pair<unsigned, T> max = { 0,0 };
      for (unsigned i = j; i < matr.size(); i++)
      {
         if (std::abs((*(matr[i]))[j]) > max.second)
         {
            max = { i,std::abs((*(matr[i]))[j]) };

         }
      }
      return max.first;
   }

public:
   GaussAlgoth(std::vector<std::vector<T> *> *matrix): matrix(matrix) {}

   void transform(std::vector<T>& vec)
   {
      auto &matr = *matrix;
      for (unsigned i = 0; i < matr.size() - 1; i++)
      {
         unsigned j = findMax(i);
         auto buff = matr[j];
         matr[j] = matr[i];
         matr[i] = buff;
         T bv = vec[j];
         vec[j] = vec[i];
         vec[i] = bv;

         for (j = i+1; j < matr.size(); j++) {
            T mult = (*(matr[j]))[i] / (*(matr[i]))[i];
            for (unsigned k = i; k < matr.size(); k++)
            {
               (*(matr[j]))[k] = (*(matr[j]))[k] - (*(matr[i]))[k] * mult;
            }
            vec[j] = vec[j] - vec[i] * mult;
         }

      }
   }

   void transform(std::vector<std::vector<T>> &B)
   {
      auto &matr = *matrix;
      for (unsigned i = 0; i < matr.size() - 1; i++)
      {
         unsigned j = findMax(i);

         std::swap(matr[i], matr[j]);


         for(auto & col : B)
         {
            std::swap(col[i], col[j]);
         }


         for (j = i + 1; j < matr.size(); j++) {
            T mult = (*(matr[j]))[i] / (*(matr[i]))[i];
            for (unsigned k = i; k < matr.size(); k++)
            {
               (*(matr[j]))[k] = (*(matr[j]))[k] - (*(matr[i]))[k] * mult;
            }
            for(auto &col : B)
            {
               col[j] = col[j] - col[i] * mult;
            }
         }
      }
      
   }
   void calcX(std::vector<T> &vec)
   {
      auto &matr = *matrix;
      for(unsigned i = vec.size()-1; i < vec.size();i--)
      {
         T buff = 0;
         for(unsigned j = vec.size() - 1; j > i;j--)
         {
            buff += (*matr[i])[j] * vec[j];
         }
         vec[i] = (vec[i] - buff) / (*matr[i])[i];
      }
   }

   std::vector<T> calcB(std::vector<T>& X)
   {
      auto &matr = *matrix;
      std::vector<T> B(X.size());
      for (unsigned i = 0; i < X.size(); i++)
      {
         T buff = 0;
         for (unsigned j = 0; j < X.size(); j++)
         {
            T b = (*matr[i])[j];
            b = X[j];
            buff += (*matr[i])[j] * X[j];
         }
         B[i] = buff;
      }
      return B;
   }

   void OutGaussRAW(std::ostream& out)
   {
      auto &matr = *matrix;
      for (unsigned i = 0; i < matr.size(); i++)
      {
         for (unsigned j = 0; j < matr.size(); j++)
         {
            out << (*(matr[i]))[j] << " ";

         }
         out << std::endl;
      }
   }

};