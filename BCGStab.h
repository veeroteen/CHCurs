#pragma once
#include "CS.h"

template <Field T>
class BCGStab :public ThreeStageBase<T>
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
	void iteration(std::vector<T> &r, std::vector<T> &p, std::vector<T> &v, std::vector<T> &r0, T &rop, T &a, T &w)
	{
		T ro = scalar(r0, r);
		T b = ro / rop * a / w;

		for (size_t i = 0; i < p.size(); i++)
		{
			p[i] = r[i] + b * (p[i] - w * v[i]);
		}

		A->multiplyA(p, v);
		a = ro / scalar(r0, v);
		std::vector<T> s = v;
		for (auto &var : s)
		{
			var *= a;
		}
		diff(r, s, s);

		if (sqrt(scalar(s, s) / scalar(*f, *f)) < eps)
		{
			for (size_t i = 0; i < x->size(); i++)
			{
				(*x)[i] += a * r[i];
			}

			r = s;
			rop = ro;

			return;
		}

		std::vector<T> t(s.size(), 0);
		A->multiplyA(s, t);
		w = scalar(t, s) / scalar(t, t);

		for (size_t i = 0; i < x->size(); i++)
		{
			(*x)[i] = (*x)[i] + a * p[i] + w * s[i];
			r[i] = s[i] - w * t[i];
			rop = ro;
		}
	}
	void iterationDiag(std::vector<T> &r, std::vector<T> &p, std::vector<T> &v, std::vector<T> &r0, T &rop, T &a, T &w)
	{
		T ro = scalar(r0, r);
		T b = ro / rop * a / w;

		for (size_t i = 0; i < p.size(); i++)
		{
			p[i] = r[i] + b * (p[i] - w * v[i]);
		}
		std::vector<T>z(p.size(), 0);
		A->diagSolve(*A->di, p, z);

		A->multiplyA(z, v);

		a = ro / scalar(r0, v);
		std::vector<T> s = v;
		for (auto &var : s)
		{
			var *= a;
		}
		diff(r, s, s);

		std::vector<T> t(s.size(), 0);
		std::vector<T> y(s.size(), 0);
		A->diagSolve(*A->di, s, y);
		A->multiplyA(y, t);
		w = scalar(t, s) / scalar(t, t);

		for (size_t i = 0; i < x->size(); i++)
		{
			(*x)[i] = (*x)[i] + a * z[i] + w * y[i];
			r[i] = s[i] - w * t[i];
			rop = ro;
		}
	}

	void iterationLU(std::vector<T> &r, std::vector<T> &p, std::vector<T> &v, std::vector<T> &r0, T &rop, T &a, T &w)
	{
		T ro = scalar(r0, r);
		T b = ro / rop * a / w;

		for (size_t i = 0; i < p.size(); i++)
		{
			p[i] = r[i] + b * (p[i] - w * v[i]);
		}
		std::vector<T>z(p.size(), 0);
		matrix.forwSolutionCSR(*matrix.il, *matrix.jl, p, z, *matrix.ll);
		matrix.revrsSolutionCSC(*matrix.iu, *matrix.ju, z, z, *matrix.lu, false);
		A->multiplyA(z, v);

		a = ro / scalar(r0, v);
		std::vector<T> s = v;
		for (auto &var : s)
		{
			var *= a;
		}
		diff(r, s, s);

		std::vector<T> t(s.size(), 0);
		std::vector<T> y(s.size(), 0);
		matrix.forwSolutionCSR(*matrix.il, *matrix.jl, s, y, *matrix.ll);
		matrix.revrsSolutionCSC(*matrix.iu, *matrix.ju, y, y, *matrix.lu, false);
		A->multiplyA(y, t);
		w = scalar(t, s) / scalar(t, t);

		for (size_t i = 0; i < x->size(); i++)
		{
			(*x)[i] = (*x)[i] + a * z[i] + w * y[i];
			r[i] = s[i] - w * t[i];
		}
		rop = ro;
	}

	
public:
	BCGStab(std::vector<size_t> *il, std::vector<size_t> *jl, std::vector<size_t> *iu, std::vector<size_t> *ju, std::vector<T> *ll, std::vector<T> *lu, std::vector<T> *di, std::vector<T> *f) : ThreeStageBase<T>(il, jl, iu, ju, ll, lu, di, f)
	{

	}
	BCGStab(std::istream &input) : ThreeStageBase<T>(input)
	{
	}
	BCGStab(std::string &path) : ThreeStageBase<T>(path,false)
	{
	}

	void Solve(unsigned par)
	{
		size_t iter = 0;
		switch (par)
		{
			case 1:
			{

				std::vector<T> r(*f);
				T a = 1, w = 1;
				T ro = 1;
				std::vector<T> p(r.size(), 0);
				std::vector<T> v(r.size(), 0);
				std::vector<T> r0 = r;

				for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
				{
					iteration(r,p,v, r0,ro,a,w);

					iter++;
				}
				break;
			}
			case 2:
			{
				std::vector<T> r(*f);
				T a = 1, w = 1;
				T ro = 1;
				std::vector<T> p(r.size(), 0);
				std::vector<T> v(r.size(), 0);
				std::vector<T> r0 = r;

				for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
				{
					iterationDiag(r, p, v, r0, ro, a, w);

					iter++;
				}
				break;
			}
			case 3:
			{
				incompLU();
				
				std::vector<T> r(*f);
				T a = 1, w = 1;
				T ro = 1;
				std::vector<T> p(r.size(), 0);
				std::vector<T> v(r.size(), 0);
				std::vector<T> r0 = r;

				for (size_t i = 0; i < iterC && sqrt(scalar(r, r) / scalar(*f, *f)) > eps; i++)
				{
					iterationLU(r, p, v, r0, ro, a, w);

					iter++;
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
			std::cout << std::setprecision(8) << 1 << " ";
			std::cout << std::endl;
		}

	}
	void printV(std::vector<T> &p)
	{
		for (auto &a : p)
		{
			std::cout << std::setprecision(16) << a << std::endl;
		}
		std::cout << std::endl;
	}

};