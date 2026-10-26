#include <solver.hpp>
#include <iostream>


int check_point_inside_region(double x, double y)
{
    if ((abs(x) >= 2) || (abs(y) >= 2) || (x*x + (y-2)*(y-2) < 1))
    {
        return -1;
    }
    if (abs(x*x + (y-2)*(y-2) - 1) < 1e-6)
    {
        return 0;
    }
    return 1;
}

double boundary_condition(double x, double y)
{
    return sin(x + y);
}

int main()
{
    size_t M = 80;
    size_t N = M;
    double x_min = -2.1;
    double x_max = 2.1;
    double y_min = -2.1;
    double y_max = 2.1;

    Solver solver(M, N, boundary_condition, check_point_inside_region, x_min, x_max, y_min, y_max);
    std::cout << "start solving" << std::endl;
    auto result = solver.solve(1e-8);

    auto file = std::ofstream("./result.txt");
    file << M << ' ' << N << ' ' << x_min << ' ' << x_max << ' ' << y_min << ' ' << y_max << '\n';
    file.close();

    result.write_to_file("./result.txt");
}