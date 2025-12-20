#pragma once
#include <array>
#include <vector>
double u(std::array<double,3> &cords)
{
	return cords[0] + cords[1] + cords[2];

}

size_t gu(std::array<double,3> &normal,std::vector<double> &poly)
{
	poly.resize(1);
	poly[0] += normal[0] * 1;
	poly[0] += normal[1] * 1;
	poly[0] += normal[2] * 1;
	return 1;
}
double dgu(std::array<double,3> &cords)
{
	return 0;
}

double u2(std::array<double, 3> &cords)
{
	return cords[0]* cords[0] + cords[1]* cords[1] + cords[2]* cords[2];

}

size_t gu2(std::array<double, 3> &normal, std::vector<double> &poly)
{
	poly.resize(4);
	poly[0] = 0;
	poly[1] = 2 * normal[0];
	poly[2] = 2 * normal[1];
	poly[3] = 2 * normal[2];
	return 4;
}

double dgu2(std::array<double, 3> &cords)
{
	return -6;
}