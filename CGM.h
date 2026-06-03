#pragma once
#include "CS.h"

template <Field T>
class CGM : public ThreeStageBase<T>
{
	using ThreeStageBase<T>::matrix;
	using ThreeStageBase<T>::f;
	using ThreeStageBase<T>::x;
	using ThreeStageBase<T>::symmetry;
	using ThreeStageBase<T>::A;
	using ThreeStageBase<T>::incompLU;
	using ThreeStageBase<T>::incompChol;
	using ThreeStageBase<T>::eps;
	using ThreeStageBase<T>::iterC;

private:


	void iterationSym(std::vector<T> &r, std::vector<T> &p, std::vector<T> &z, CMatrix<T> &M)
	{
		std::vector<T> temp(z);
		A->multiplyA(temp); // Az_k-1

		T b = scalar(p, r);
		T a = b / scalar(temp, z);

		for (size_t i = 0; i < f->size(); i++)
		{
			(*x)[i] = (*x)[i] + a * z[i];
			r[i] = r[i] - a * temp[i];
		}

		matrix.forwSolutionCSR(*M.il, *M.jl, r, p, *M.ll);
		matrix.revrsSolutionCSC(*M.il, *M.jl, p, p, *M.ll);

		b = scalar(p, r) / b;
		for (size_t i = 0; i < f->size(); i++)
		{
			z[i] = p[i] + b * z[i];
		}

	}
	void iterationDiag(std::vector<T> &r, std::vector<T> &p, std::vector<T> &z, std::vector<T> &diag)
	{
		std::vector<T> temp(z);
		A->multiplyA(temp);
		T b = scalar(p, r);
		T a = b / scalar(temp, z);

		for (size_t i = 0; i < f->size(); i++)
		{
			(*x)[i] = (*x)[i] + a * z[i];
			r[i] = r[i] - a * temp[i];
		}
		A->diagSolve(diag, r, p);

		b = scalar(p, r) / b;
		for (size_t i = 0; i < f->size(); i++)
		{
			z[i] = p[i] + b * z[i];
		}

	}
	void iteration(std::vector<T> &r, std::vector<T> &z)
	{
		std::vector<T> temp(z);
		A->multiplyA(temp);

		T rr = scalar(r, r);
		T a = rr / scalar(temp, z);

		for (size_t i = 0; i < z.size(); i++)
		{
			(*x)[i] = (*x)[i] + a * z[i];
			r[i] = r[i] - a * temp[i];
		}

		T b = scalar(r, r) / rr;
		for (size_t i = 0; i < z.size(); i++)
		{
			z[i] = r[i] + b * z[i];
		}
	}
	void iterationLU(std::vector<T> &r, std::vector<T> &z)
	{
		std::vector<T> tmp(z);
		matrix.revrsSolutionCSC(*matrix.iu, *matrix.ju, tmp, tmp, *matrix.lu, false);
		A->multiplyA(tmp);
		matrix.forwSolutionCSR(*matrix.il, *matrix.jl, tmp, tmp, *matrix.ll);
		matrix.revrsSolutionCSC(*matrix.il, *matrix.jl, tmp, tmp, *matrix.ll);
		A->multiplyAT(tmp);
		matrix.forwSolutionCSR(*matrix.iu, *matrix.ju, tmp, tmp, *matrix.lu, false);


		T rr = scalar(r, r);
		T a = rr / scalar(tmp, z);

		for (size_t i = 0; i < z.size(); i++)
		{
			(*x)[i] = (*x)[i] + a * z[i];
			r[i] = r[i] - a * tmp[i];
		}
		T b = scalar(r, r) / rr;
		for (size_t i = 0; i < z.size(); i++)
		{
			z[i] = r[i] + b * z[i];
		}


	}

public:
	CGM(std::vector<size_t> *il, std::vector<size_t> *jl, std::vector<size_t> *iu,std::vector<size_t>*ju, std::vector<T> *ll, std::vector<T> *lu, std::vector<T> *di, std::vector<T> *f) : ThreeStageBase<T>(il, jl,iu,ju, ll, lu, di, f)
	{

	}
	CGM(std::istream &input) : ThreeStageBase<T>(input)
	{
	}
	CGM(std::string &path) : ThreeStageBase<T>(path)
	{
	}

	void factorizate()
	{
		if (symmetry)
		{
			incompChol();
		}
		else
		{
			incompLU();
		}
	}

	void multiplyTriangle(std::vector<T> &x, std::vector<T> &al, std::vector<size_t> &ia, std::vector<size_t> &ja, std::vector<T> &di)
	{
		std::vector tx(x);

		for (size_t i = 0; i < x.size(); i++)
		{
			T buff = (di)[i] * tx[i];
			for (size_t k = (ia)[i]; k < (ia)[i + 1]; k++)
			{
				x[ja[k]] += (al)[k] * tx[i];
			}
		}
	}

	void Solve(unsigned par)
	{
		this->nullifyX();
		size_t iter = 0;
		switch (par)
		{
			case 1:
			{
				if (symmetry)
				{

					std::vector<T> r(*f);
					std::vector<T> z(r);

					for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
					{
						iteration(r, z);
						iter++;
					}
				}
				break;
			}
			case 2:
			{

				std::vector<T> r(*f);
				std::vector<T> z(r.size(), 0);
				matrix.diagSolve(*matrix.di, r, z);

				std::vector<T> p(z);

				for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
				{
					iterationDiag(r, p, z, *matrix.di);
					iter++;
				}
				break;
			}
			case 3:
			{
				factorizate();

				if (symmetry)
				{
					std::vector<T> r(*f);
					std::vector<T> z(r.size(), 0);
					matrix.forwSolutionCSR(*matrix.il, *matrix.jl, r, z, *matrix.ll);
					matrix.revrsSolutionCSC(*matrix.il, *matrix.jl, z, z, *matrix.ll);
					std::vector<T> p(z);

					for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
					{
						iterationSym(r, p, z, matrix);
						iter++;
					}
				}
				else
				{
					std::vector<T> r(*f);
					matrix.forwSolutionCSR(*matrix.il, *matrix.jl, r, r, *matrix.ll);
					matrix.revrsSolutionCSC(*matrix.il, *matrix.jl, r, r, *matrix.ll);
					A->multiplyAT(r);
					matrix.forwSolutionCSR(*matrix.iu, *matrix.ju, r, r, *matrix.lu, false);
					std::vector<T> z(r);

					multiplyTriangle(*x, *matrix.lu, *matrix.iu, *matrix.ju, *matrix.di);


					for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
					{

						iterationLU(r, z);
						iter++;

						std::vector<T> tmp(x->size(), 0);
						matrix.revrsSolutionCSC(*matrix.iu, *matrix.ju, *x, tmp, *matrix.lu, false);

					}

					matrix.revrsSolutionCSC(*matrix.iu, *matrix.ju, *x, *x, *matrix.lu, false);
				}

				break;
			}

		}
	}
	void printTL()
	{
		for (size_t i = 0; i < A->di->size(); i++)
		{
			size_t count = i;
			size_t j = 0;
			for (size_t il = (*matrix.il)[i]; il < (*matrix.il)[i + 1]; il++)
			{
				for (; j < (*matrix.jl)[il]; j++)
				{
					std::cout << 0 << " ";
					count--;
				}
				std::cout << std::setprecision(8) << (*matrix.ll)[il] << " ";
				j++;
				count--;

			}
			for (size_t t = 0; t < count; t++)
			{
				std::cout << 0 << " ";
			}
			std::cout << std::setprecision(8) << (*matrix.di)[i] << " ";
			std::cout << std::endl;
		}



	}
	void printV(std::vector<T> &p)
	{
		for (auto &a : p)
		{
			std::cout << std::setprecision(16) << a << std::endl;
		}

	}
};

