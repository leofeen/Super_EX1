#pragma once

#include <matrix.hpp>
#include <numbers>
#include <string>
#include <memory>
#include <iostream>
#include <cmath>

class Solver
{
public: 
    Solver() = default;
    
    Solver(size_t M, size_t N, double (*boundary_condition)(double, double), int (*check_point_inside_region)(double, double), double x_min, double x_max, double y_min, double y_max);

    Matrix<double> solve(double precision);

private:
    double x_min;
    double h_x;
    size_t M;

    double y_min;
    double h_y;
    size_t N;

    double eps;

    Matrix<double> Q;
    Matrix<double> F;
    Matrix<double> S;
    Matrix<double> B;
    Matrix<double> D;

    Matrix<double> operator_A(Matrix<double> iteration);

    double get_x(size_t i);
    double get_y(size_t j);

    // double (*density)(double, double) = nullptr;  // We have f = 1. Why bother?

    double (*boundary_condition)(double, double) = nullptr;
    int (*check_point_inside_region)(double, double) = nullptr;

    double scalar_product(Matrix<double> A, Matrix<double> B);
    double norm(Matrix<double> A);

    double dichotomy_solve_x(double x_start, double x_end, double y);
    double dichotomy_solve_y(double y_start, double y_end, double x);
    void init_S();
};
