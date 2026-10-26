#include "solver.hpp"


Solver::Solver(size_t M, size_t N, double (*boundary_condition)(double, double), int (*check_point_inside_region)(double, double), double x_min, double x_max, double y_min, double y_max) :
M(M), N(N), x_min(x_min), y_min(y_min), boundary_condition(boundary_condition), check_point_inside_region(check_point_inside_region)
{
    if ((boundary_condition == nullptr) || (check_point_inside_region == nullptr))
    {
        std::cout << "nullptr passed" << '\n';
        std::exit(1);
    }

    this->h_x = (x_max - x_min) / M;
    this->h_y = (y_max - y_min) / N;
    this->eps = pow(std::max(this->h_x, this->h_y), 2);

    this->init_S();

    this->Q = (Matrix<double>(M+1, N+1, this->h_x*this->h_y) - this->S) / (this->h_x*this->h_y * pow(this->eps, 2));
    // this->F = this->S / (this->h_x*this->h_y);
    this->F = Matrix<double>(M+1, N+1, 1); 

    this->B = Matrix<double>(N+1, M+1);
    #pragma omp parallel for simd
    for (size_t i = 1; i <= M - 1; i++)
    {
        for (size_t j = 1; j <= N - 1; j++)
        {
            this->B[i, j] = F[i, j] + Q[i, j]*this->boundary_condition(this->get_x(i), this->get_y(j));
        }
    }

    this->D = this->Q + Matrix<double>(this->M+1, this->N+1, 2.0/pow(this->h_x, 2) + 2.0/pow(this->h_y, 2));
}

Matrix<double> Solver::solve(double precision, size_t* iteration_count)
{
    double beta;
    Matrix<double> current_iteration(this->M+1, this->N+1);

    Matrix<double> next_r(this->M+1, this->N+1);
    auto current_r = this->B - this->operator_A(current_iteration);

    Matrix<double> next_z(this->M+1, this->N+1);
    auto current_z = current_r.hadamard_division(this->D);

    auto descent = current_z;

    auto step = this->scalar_product(current_z, current_r) / this->scalar_product(this->operator_A(descent), descent);

    auto next_iteration = current_iteration + step*descent;
    *iteration_count = 1;

    while (this->norm(next_iteration - current_iteration) >= precision)
    {
        current_iteration = next_iteration;

        next_r = current_r - step*this->operator_A(descent);
        next_z = next_r.hadamard_division(this->D);
        beta = this->scalar_product(next_z, next_r) / this->scalar_product(current_z, current_r);

        current_r = next_r;
        current_z = next_z;

        descent = current_z + beta*descent;
        step = this->scalar_product(current_z, current_r) / this->scalar_product(this->operator_A(descent), descent);

        next_iteration = current_iteration + step*descent;
        *iteration_count += 1;
    }

    return next_iteration;
}

Matrix<double> Solver::operator_A(Matrix<double> iteration)
{
    auto result = Matrix<double>(iteration.get_rows(), iteration.get_columns());
    #pragma omp parallel for simd
    for (size_t i = 1; i < result.get_rows()-1; i++)
    {
        for (size_t j = 1; j < result.get_columns()-1; j++)
        {
            result[i, j] = -(iteration[i+1, j] - 2*iteration[i, j] + iteration[i-1, j]) / pow(this->h_x, 2) - (iteration[i, j+1] - 2*iteration[i, j] + iteration[i, j-1]) / pow(this->h_y, 2) + Q[i, j]*iteration[i, j];
        }
    }

    return result;
}

double Solver::get_x(size_t i)
{
    return this->x_min + i*this->h_x;
}

double Solver::get_y(size_t j)
{
    return this->y_min + j*this->h_y;
}

double Solver::scalar_product(Matrix<double> A, Matrix<double> B)
{
    double result = 0;
    // #pragma omp parallel for reduction(+ : result)
    for (size_t i = 1; i <= this->M - 1; i++)
    {
        for (size_t j = 1; j <= this->N - 1; j++)
        {
            result += this->h_x * this->h_y * A[i, j] * B[i, j];
        }
    }

    return result;
}

double Solver::norm(Matrix<double> A)
{
    return sqrt(this->scalar_product(A, A));
}

double Solver::dichotomy_solve_x(double x_start, double x_end, double y)
{
    double x_left = x_start;
    double x_right = x_end;
    
    int f_left = this->boundary_condition(x_start, y);
    if (f_left == 0) {return x_left;}

    int f_right = this->boundary_condition(x_end, y);
    if (f_right == 0) {return x_right;}

    double x_centre = (x_start + x_end) / 2.0;
    int f_centre;
    while (x_right - x_left >= 1e-6)
    {
        x_centre = (x_left + x_right) / 2.0;
        f_centre = this->boundary_condition(x_centre, y);
        if (f_centre == 0) return x_centre;

        if (f_left * f_centre < 0)
        {
            f_right = f_centre;
            x_right = x_centre;
        }

        if (f_centre * f_right < 0)
        {
            f_left = f_centre;
            x_left = x_centre;
        }
    }

    return x_centre;
}

double Solver::dichotomy_solve_y(double y_start, double y_end, double x)
{
    double y_left = y_start;
    double y_right = y_end;
    
    int f_left = this->boundary_condition(x, y_start);
    if (f_left == 0) {return y_left;}

    int f_right = this->boundary_condition(x, y_end);
    if (f_right == 0) {return y_right;}

    double y_centre = (y_start + y_end) / 2.0;
    int f_centre;
    while (y_right - y_left >= 1e-6)
    {
        y_centre = (y_left + y_right) / 2.0;
        f_centre = this->boundary_condition(x, y_centre);
        if (f_centre == 0) return y_centre;

        if (f_left * f_centre < 0)
        {
            f_right = f_centre;
            y_right = y_centre;
        }

        if (f_centre * f_right < 0)
        {
            f_left = f_centre;
            y_left = y_centre;
        }
    }

    return y_centre;
}

void Solver::init_S()
{
// This function checks for all 14 cases of (P_ij intersects D) 
    // double x_left, x_right, y_down, y_up, zero1, zero2, s_aux;
    // auto inside_region = Matrix<int>(2, 2);
    this->S = Matrix<double>(this->M+1, this->N+1);
    #pragma omp parallel for
    for (size_t i = 1; i <= this->M - 1; i++)
    {
        for (size_t j = 1; j <= this->N - 1; j++)
        {
            double x_left, x_right, y_down, y_up, zero1, zero2, s_aux;
            auto inside_region = Matrix<int>(2, 2);

            x_left = this->get_x(i) - 0.5*this->h_x;
            x_right = this->get_x(i) + 0.5*this->h_x;

            y_down = this->get_y(j) - 0.5*this->h_y;
            y_up = this->get_y(j) + 0.5*this->h_y;

            inside_region[0, 0] = this->check_point_inside_region(x_left, y_down);
            inside_region[0, 1] = this->check_point_inside_region(x_left, y_up);
            inside_region[1, 1] = this->check_point_inside_region(x_right, y_up);
            inside_region[1, 0] = this->check_point_inside_region(x_right, y_down);

            if ((inside_region[0, 0] >= 0) && (inside_region[0, 1] >= 0) && (inside_region[1, 1] >= 0) && (inside_region[1, 0] >= 0))
            {
                this->S[i, j] = this->h_x*this->h_y;
                continue;
            }

            if ((inside_region[0, 0] <= 0) && (inside_region[0, 1] <= 0) && (inside_region[1, 1] <= 0) && (inside_region[1, 0] <= 0))
            {
                this->S[i, j] = 0;
                continue;
            }

            for (int i1 = 0; i1 <= 1; i1++)
            {
                for (int j1 = 0; j1 <= 1; j1++)
                {
                    inside_region[i1, j1] += 1 - abs(inside_region[i1, j1]); // Turn all 0 into 1, dichotomy will deal with it automatically
                }
            }

            if ((inside_region[0, 0]*inside_region[1, 0] == -1) && (inside_region[0, 1]*inside_region[1, 1] == -1))
            {
                zero1 = this->dichotomy_solve_x(x_left, x_right, y_down);
                zero2 = this->dichotomy_solve_x(x_left, x_right, y_up);

                s_aux = this->h_y * (zero1 - x_left + zero2 - x_left)/2.0;

                if (inside_region[0, 0] == 1)
                {
                    this->S[i, j] = s_aux;
                    continue;
                }

                this->S[i, j] = this->h_x*this->h_y - s_aux;
                continue;
            }

            if ((inside_region[0, 0]*inside_region[0, 1] == -1) && (inside_region[1, 0]*inside_region[1, 1] == -1))
            {
                zero1 = this->dichotomy_solve_y(y_down, y_up, x_left);
                zero2 = this->dichotomy_solve_y(y_down, y_up, x_right);

                s_aux = this->h_x * (zero1 - y_down + zero2 - y_down)/2.0;

                if (inside_region[0, 0] == 1)
                {
                    this->S[i, j] = s_aux;
                    continue;
                }

                this->S[i, j] = this->h_x*this->h_y - s_aux;
                continue;
            }

            if ((inside_region[0, 0]*inside_region[0, 1] == -1) && (inside_region[0, 0]*inside_region[1, 0] == -1))
            {
                zero1 = this->dichotomy_solve_x(x_left, x_right, y_down);
                zero2 = this->dichotomy_solve_y(y_down, y_up, x_left);

                s_aux = 0.5 * (zero1 - x_left) * (zero2 - y_down);

                if (inside_region[0, 0] == 1)
                {
                    this->S[i, j] = s_aux;
                    continue;
                }

                this->S[i, j] = this->h_x*this->h_y - s_aux;
                continue;
            }

            if ((inside_region[0, 0]*inside_region[0, 1] == -1) && (inside_region[0, 1]*inside_region[1, 1] == -1))
            {
                zero1 = this->dichotomy_solve_x(x_left, x_right, y_up);
                zero2 = this->dichotomy_solve_y(y_down, y_up, x_left);

                s_aux = 0.5 * (zero1 - x_left) * (y_up - zero2);

                if (inside_region[0, 1] == 1)
                {
                    this->S[i, j] = s_aux;
                    continue;
                }

                this->S[i, j] = this->h_x*this->h_y - s_aux;
                continue;
            }

            if ((inside_region[1, 0]*inside_region[1, 1] == -1) && (inside_region[0, 1]*inside_region[1, 1] == -1))
            {
                zero1 = this->dichotomy_solve_x(x_left, x_right, y_up);
                zero2 = this->dichotomy_solve_y(y_down, y_up, x_right);

                s_aux = 0.5 * (x_right - zero1) * (y_up - zero2);

                if (inside_region[1, 1] == 1)
                {
                    this->S[i, j] = s_aux;
                    continue;
                }

                this->S[i, j] = this->h_x*this->h_y - s_aux;
                continue;
            }

            if ((inside_region[1, 0]*inside_region[1, 1] == -1) && (inside_region[1, 0]*inside_region[0, 0] == -1))
            {
                zero1 = this->dichotomy_solve_x(x_left, x_right, y_down);
                zero2 = this->dichotomy_solve_y(y_down, y_up, x_right);

                s_aux = 0.5 * (x_right - zero1) * (zero2 - y_down);

                if (inside_region[1, 0] == 1)
                {
                    this->S[i, j] = s_aux;
                    continue;
                }

                this->S[i, j] = this->h_x*this->h_y - s_aux;
                continue;
            }
        }
    }
}
