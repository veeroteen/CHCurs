#pragma once
#include <array>
#include <string>
#include <vector>
double u(std::array<double,4> &cords)
{
	return cords[0] +cords[1] + cords[2];

}
void ustr(std::string &fun)
{
	fun = "x+y+z+sin(t)";
}
void fu(std::string &fun)
{
	fun = "cos(t)";
}

void gu(std::array<double,3> &normal,std::string &poly)
{
	double a = 0;
	a += normal[0] * 1;
	a += normal[1] * 1;
	a += normal[2] * 1;
	poly = std::to_string(a);
	return;
}

double u2(std::array<double, 4> &cords)
{
	return cords[0]* cords[0] + cords[1]* cords[1] + cords[2]* cords[2];

}
void u2str(std::string &fun)
{
	fun = "x^2 + y^2 + z^2";
}

void gu2(std::array<double, 3> &normal, std::string &poly)
{
	if (normal[0] != 0)
	{
		poly += std::to_string(normal[0]) + "*2*x";
	}
	if (normal[1] != 0)
	{
		poly += std::to_string(normal[1]) + "*2*y";
	}
	if (normal[2] != 0)
	{
		poly += std::to_string(normal[2]) + "*2*z";
	}
}

double dgu2(std::array<double, 4> &cords)
{
	return -6;
}

double usin(std::array<double, 4> &cords)
{
	return sin(cords[0] + cords[1] + cords[2]);

}
void ustrsin(std::string &fun)
{
	fun = "sin(x + y + z)";
}

void gusin(std::array<double, 3> &normal, std::string &poly)
{
	if (normal[0] != 0) 
	{
		poly += std::to_string(normal[0]) + "*cos(x+y+z)";
	}
	if (normal[1] != 0)
	{
		poly += std::to_string(normal[1]) + "*cos(x+y+z)";
	}
	if (normal[2] != 0)
	{
		poly += std::to_string(normal[2]) + "*cos(x+y+z)";
	}

	return;
}
double dgusin(std::array<double, 4> &cords)
{
	return 3*sin(cords[0] + cords[1] + cords[2]);
}

double u3(std::array<double, 4> &cords)
{
	return cords[0] * cords[0]* cords[0] + cords[1]* cords[1] * cords[1] + cords[2] * cords[2]* cords[2];

}
void u3str(std::string &fun)
{
	fun = "x^3 + y^3 + z^3";
}

void gu3(std::array<double, 3> &normal, std::string &poly)
{
	if (normal[0] != 0)
	{
		poly += std::to_string(normal[0]) + "*3x^2";
	}
	if (normal[1] != 0)
	{
		poly += std::to_string(normal[1]) + "*3y^2";
	}
	if (normal[2] != 0)
	{
		poly += std::to_string(normal[2]) + "*3z^2";
	}
}

double dgu3(std::array<double, 4> &cords)
{
	return -6*(cords[0] + cords[1] + cords[2]);
}