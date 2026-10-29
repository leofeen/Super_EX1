#include <solver.hpp>
#include <iostream>
#include <chrono>
#include <omp.h>


#define NUM_THREADS 16
#define N_POINTS 80


int check_point_inside_region(double x, double y)
{
    if ((std::abs(x) >= 2) || (std::abs(y) >= 2) || (x*x + (y-2)*(y-2) < 1))
    {
        return -1;
    }
    if (std::abs(x*x + (y-2)*(y-2) - 1) < 1e-6)
    {
        return 0;
    }
    return 1;
}

double boundary_condition(double x, double y)
{
    return sin(x + y);
}

int main(int argc, char* argv[])
{
    // omp_set_num_threads(std::stoi(argv[2]));
    size_t M, N;
    int index = 0;
    if (argc == 4)
    {
        M = std::stoi(argv[1]);
        N = std::stoi(argv[2]);
        index = std::stoi(argv[3]);
    }
    else
    {
        M = N_POINTS;
        N = M;
        omp_set_num_threads(NUM_THREADS);
    }

    double x_min = -2.1;
    double x_max = 2.1;
    double y_min = -2.1;
    double y_max = 2.1;

    size_t iteration_count = 1;

    #ifdef _OPENMP
    auto start_time = omp_get_wtime();
    #else
    auto start_time = std::chrono::high_resolution_clock::now();
    #endif

    Solver solver(M, N, boundary_condition, check_point_inside_region, x_min, x_max, y_min, y_max);
    std::cout << "start solving" << std::endl;
    auto result = solver.solve(1e-8, &iteration_count);

    #ifdef _OPENMP
    auto end_time = omp_get_wtime();
    auto duration = end_time - start_time;
    #else
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time) / 1e6;
    #endif

    // auto file = std::ofstream("./result_main.txt");
    // file << duration << ' ' << iteration_count << '\n';
    // file << M << ' ' << N << ' ' << x_min << ' ' << x_max << ' ' << y_min << ' ' << y_max << '\n';
    // file.close();

    char filename[40];
    std::sprintf(filename, "./result_%zu_%zu_%d.txt", M, N, index);
    FILE* file = fopen(filename, "w");
    std::fprintf(file, "%f %zu\n", duration, iteration_count);
    std::fprintf(file, "%zu %zu %f %f %f %f\n", M, N, x_min, x_max, y_min, y_max);
    fclose(file);

    result.write_to_file(filename);
}