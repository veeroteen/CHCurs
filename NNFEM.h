#pragma once
#include "FEM.h"


template<Field T, ElemType El>
class NNFEM : public FEM <T, El >
{
   using FEM<T, El>::elements;
   using FEM<T, El>::nodes;
   using FEM<T, El>::u;
   using FEM<T, El>::f;
   using FEM<T, El>::size;
   using FEM<T, El>::gamma;

   using FEM<T, El>::SLAU;
   using FEM<T, El>::Neumann;
   using FEM<T, El>::Robin;
   using FEM<T, El>::Dirih;
   using FEM<T, El>::setNodes;
   using FEM<T, El>::setElements;
   using FEM<T, El>::setF;
private:
   void NewtonAdd(CMatrix<T> &matrix)
   {
      auto &dir = *matrix.di;
      auto &fr = *f;

      auto &ilr = *matrix.il;
      auto &jlr = *matrix.jl;
      auto &llr = *matrix.ll;

      auto &iur = *matrix.iu;
      auto &jur = *matrix.ju;
      auto &lur = *matrix.lu;
      std::vector<StringFun<T>> gammaTerm;
      bool gam = false;
      std::string buff = gamma.getStringFun();
      auto a = buff.find("u");
      if (a != buff.npos)
      {
         std::string base2 = buff;
         buff.replace(a, 1, "(u+(1E-8))");
         base2.replace(a, 1, "(u-(1E-8))");
         buff = "((" + buff + ")-(" + base2 + "))/" + "(2 * (1E-8))";
         gammaTerm.emplace_back(buff);
         buff = "u";
         gammaTerm.emplace_back(buff);
         gam = true;
      }
      else
      {
         gam = false;
      }

      for (size_t e = 0; e < elements.size(); e++)
      {
         bool lam = false;
         std::string base = "";
         std::vector<StringFun<T>> lambdaTerm;
         std::vector<T> local(El::GetNodesCount() * El::GetNodesCount(), 0);

         base = elements[e].h.getStringFun();
         auto a = base.find("u");

         if (a != base.npos)
         {
            std::string base2 = base;
            base.replace(a, 1, "(u+(1E-8))");
            base2.replace(a, 1, "(u-(1E-8))");
            base = "((" + base + ")-(" + base2 + "))/" + "(2 * (1E-8))";
            lambdaTerm.emplace_back(base);
            lam = true;
         }
         else
         {
            lam = false;
         }

         if (gam || lam)
         {
            if (lam)
            {
               std::array<T, El::GetDim()> gradu{};
               for (size_t k = 0; k < El::GetNodesCount(); k++)
               {
                  for (size_t g = 0; g < El::GetDim(); g++)
                  {
                     gradu[g] += u[elements[e].vertices[k]] * elements[e].basis[k][g + 1];
                  }
               }
               for (size_t i = 0; i < El::GetNodesCount(); i++)
               {
                  T res = 0;
                  for (size_t k = 0; k < El::GetDim(); k++)
                  {
                     res += gradu[k] * elements[e].basis[i][k + 1];
                  }

                  for (size_t j = 0; j < El::GetNodesCount(); j++)
                  {
                     polyToStr(base, elements[e].basis[j].basis);
                     lambdaTerm.emplace_back(base);
                     local[j + i * El::GetNodesCount()] = res * integrate(elements[e], elements[e].vertices, nodes, elements[e].V, lambdaTerm);
                     lambdaTerm.pop_back();
                  }

               }
            }
            if (gam)
            {
               for (size_t i = 0; i < El::GetNodesCount(); i++)
               {
                  polyToStr(base, elements[e].basis[i].basis);
                  gammaTerm.emplace_back(base);
                  for (size_t j = 0; j < El::GetNodesCount(); j++)
                  {
                     polyToStr(base, elements[e].basis[j].basis);
                     gammaTerm.emplace_back(base);
                     local[j + i * El::GetNodesCount()] = integrate(elements[e], elements[e].vertices, nodes, elements[e].V, gammaTerm);
                     gammaTerm.pop_back();
                  }
                  gammaTerm.pop_back();
               }
            }
            for (size_t i = 0; i < elements[e].vertices.size(); i++)
            {
               dir[elements[e].vertices[i]] += local[i + i * El::GetNodesCount()];
               for (size_t j = 0; j < i; j++)
               {
                  for (size_t it = ilr[elements[e].vertices[i]]; it < ilr[elements[e].vertices[i] + 1]; it++)
                  {
                     if (jlr[it] == elements[e].vertices[j])
                     {
                        llr[it] += local[j + i * El::GetNodesCount()];
                     }
                     if (jur[it] == elements[e].vertices[j])
                     {
                        lur[it] += local[i + j * El::GetNodesCount()];
                     }
                  }
               }

            }
         }
      }

   }

public:
   using FEM<T, El>::printU;
   NNFEM(std::string &confPath) : FEM<T, El>(confPath) {};
   void Dirih(size_t count, std::vector<T> &dir, std::vector<T> &llr, std::vector<T> &lur, std::vector<size_t> &ilr, std::vector<size_t> &jlr, std::vector<size_t> &iur, std::vector<size_t> &jur)
   {
      std::ifstream file(FEM<T, El>::dirihPath);
      if (file.is_open())
      {
         auto &fr = *f;
         size_t node;
         T value;
         for (size_t i = 0; i < count; i++)
         {
            file >> node >> value;
            fr[node] = value;
            dir[node] = 1;
            u[node] = value;
            for (size_t it = ilr[node]; it < ilr[node + 1]; it++)
            {
               fr[jlr[it]] -= value * llr[it];
               llr[it] = 0;
            }

            for (size_t p = node + 1; p < dir.size(); p++)
            {
               for (size_t it = ilr[p]; it < ilr[p + 1] && jlr[it] <= node; it++)
               {
                  if (jlr[it] == node)
                  {
                     fr[p] = fr[p] - value * llr[it];
                     llr[it] = 0;
                  }
               }

            }

            for (size_t it = iur[node]; it < iur[node + 1]; it++)
            {
               lur[it] = 0;
            }

            for (size_t p = node + 1; p < dir.size(); p++)
            {
               for (size_t it = iur[p]; it < iur[p + 1] && jur[it] <= node; it++)
               {
                  if (jur[it] == node)
                  {
                     lur[it] = 0;
                  }
               }

            }

         }
         file.close();
      }
   }
   void Solve(unsigned solven)
   {
      StringFun<T> fl;
      setF(fl);

      f = new std::vector<T>(size);
      nodes.resize(size);
      u.resize(size);
      setNodes();

      std::vector<std::set<size_t>> dict(size);
      size_t llcount = 0;
      setElements(this->elemCount, dict, llcount);


      std::vector<size_t> *il = new std::vector<size_t>(size + 1, 0);
      std::vector<size_t> *jl = new std::vector<size_t>(llcount, 0);
      auto &ilr = *il;
      auto &jlr = *jl;
      auto &fr = *f;
      size_t jli = 0;
      for (size_t a = 0; a < dict.size(); a++)
      {
         ilr[a + 1] = dict[a].size() + ilr[a];
         for (auto &n : dict[a])
         {
            jlr[jli] = n;
            jli++;
         }
      }
      std::vector<T> *di = new std::vector<T>(size, 0);
      std::vector<T> *ll = new std::vector<T>(llcount, 0);
      auto &llr = *ll;
      auto &dir = *di;

      T delta = T(0);
      size_t count = 0;

      while (true)
      {
         std::vector<StringFun<T>> polys;
         polys.push_back(gamma);
         for (size_t e = 0; e < elements.size(); e++)
         {
            std::vector<T> local(El::GetNodesCount() * (El::GetNodesCount() + 1) / 2, 0);
            for (size_t i = 0; i < El::GetNodesCount(); i++)
            {

               std::string base = "";
               polyToStr(base, elements[e].basis[i].basis);
               polys.emplace_back(base);

               for (size_t j = 0; j <= i; j++)
               {
                  polyToStr(base, elements[e].basis[j].basis);
                  polys.emplace_back(base);
                  std::vector<StringFun<T>> grads;
                  std::string base = std::to_string(elements[e].scalGrad(i, j));
                  grads.emplace_back(base);
                  grads.emplace_back(elements[e].h);
                  local[((i * (i + 1)) / 2) + j] += integrate(elements[e], elements[e].vertices, nodes, elements[e].V, grads);
                  local[((i * (i + 1)) / 2) + j] += integrate(elements[e], elements[e].vertices, nodes, elements[e].V, polys);
                  polys.pop_back();
               }
               polys.push_back(fl);
               fr[elements[e].vertices[i]] += integrate(elements[e], elements[e].vertices, nodes, elements[e].V, polys);
               polys.pop_back();
               polys.pop_back();
            }

            for (size_t i = 0; i < elements[e].vertices.size(); i++)
            {
               dir[elements[e].vertices[i]] += local[((i * (i + 1)) / 2) + i];
               for (size_t j = 0; j <= i; j++)
               {
                  if (elements[e].vertices[j] < elements[e].vertices[i])
                  {
                     for (size_t it = ilr[elements[e].vertices[i]]; it < ilr[elements[e].vertices[i] + 1]; it++)
                     {
                        if (jlr[it] == elements[e].vertices[j])
                        {
                           llr[it] += local[((i * (i + 1)) / 2) + j];
                        }
                     }
                  }
               }

            }
         }

         Neumann(this->neumann);
         Robin(this->robin, dir, llr, ilr, jlr);

         std::vector<T> *lu = new std::vector<T>();
         std::vector<size_t> *iu = new std::vector<size_t>(), *ju = new std::vector<size_t>();
         *lu = llr;
         *iu = ilr;
         *ju = jlr;
         CMatrix<T> matrix(il, jl, iu, ju, ll, lu, di, false);
         auto &lur = *lu;
         auto &iur = *iu;
         auto &jur = *ju;
         Dirih(this->dirih, dir, llr, lur, ilr, jlr, iur, jur);
         std::vector<T> Re(u.size(), 0);
         matrix.multiplyA(u, Re);
         diff(Re, fr, Re);
         NewtonAdd(matrix);
         Dirih(this->dirih, dir, llr, lur, ilr, jlr, iur, jur);

         for (auto &a : Re)
         {
            a = -a;
         }
         delta = sqrt(scalar(Re, Re));
         std::cout << delta << std::endl;
         if (delta < 1E-14)
         {
            break;
         }
         BCGStab<T> *LAU = new BCGStab<T>(il, jl, iu, ju, ll, lu, di, &Re);

         LAU->Solve(1);
         std::vector<T> du = LAU->getX();

         for (size_t i = 0; i < nodes.size(); i++)
         {
            u[i] += du[i];
         }
         for (auto &a : llr)
         {
            a = T(0);
         }
         for (auto &a : dir)
         {
            a = T(0);
         }
         for (auto &a : fr)
         {
            a = T(0);
         }
         count++;


      }

      std::cout << count << std::endl;
   }
};


