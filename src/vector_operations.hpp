#pragma once

#include<iostream>
#include<cmath>
#include<vector>

std::vector<double>operator*(const std::vector<double> & x, const double & num); //vector * number

std::vector<double>operator*(const double & num, const std::vector<double> & x); //number * vector

double operator*(const std::vector<double> & x, const std::vector<double> & y); //vector * vector

std::vector<double>operator+(const std::vector<double> & x, const std::vector<double> & y); // vector + vector

std::vector<double>operator-(const std::vector<double> & x, const std::vector<double> & y); // vector - vector

double vlength(const std::vector<double> & x); // lenght of vector

void show(const std::vector<double> & result); // just shows elements of vector