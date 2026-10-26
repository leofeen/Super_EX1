#pragma once

#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <string>

template <class T>
class Matrix
{
private:
    size_t rows;
    size_t columns;

public:
    std::vector<T> data;

    size_t get_rows() const noexcept {return rows;}
    size_t get_columns() const noexcept {return columns;}
    
    T& operator[] (size_t row, size_t column)
    {
        return data[row*this->columns + column];
    };

    const T& operator[] (size_t row, size_t column) const
    {
        return data[row*this->columns + column];
    };

    T get(size_t row, size_t column)
    {
        return (*this)[row, column];
    }

    void set(size_t row, size_t column, T value)
    {
        (*this)[row, column] = value;
    }

    Matrix() = default;

    Matrix(size_t rows, size_t columns, T default_val = 0)
    {
        data.resize(rows*columns);
        this->rows = rows;
        this->columns = columns;

        // #pragma omp parallel for simd
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < columns; j++)
            {
                (*this)[i, j] = default_val;
            }
        }
    };


    Matrix<T> operator+ (const Matrix<T>& other)
    {
        if (columns != other.columns && rows != other.rows)
        {
            std::cout << "Two matricies don\'t have right shapes for summation" << std::endl;
            exit(1);
        }

        auto result = Matrix<T>(rows, columns);
        // #pragma omp parallel for simd
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < columns; j++)
            {
                result[i, j] = (*this)[i, j] + other[i, j];
            }
        }

        return result;
    }

    Matrix<T> operator- (const Matrix<T>& other)
    {
        if (columns != other.columns && rows != other.rows)
        {
            std::cout << "Two matricies don\'t have right shapes for finding difference" << std::endl;
            exit(1);
        }

        auto result = Matrix<T>(rows, columns);

        // #pragma omp parallel for simd
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < columns; j++)
            {
                result[i, j] = (*this)[i, j] - other[i, j];
            }
        }

        return result;
    }

    Matrix<T> operator* (const T& scalar)
    {
        auto result = Matrix<T>(rows, columns);

        // #pragma omp parallel for simd
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < columns; j++)
            {
                result[i, j] = (*this)[i, j] * scalar;
            }
        }

        return result;
    }

    Matrix<T> operator/ (const T& scalar)
    {
        return (*this) * (1.0 / scalar);
    }

    Matrix<T> operator* (const Matrix<T>& other)
    {
        if (columns != other.rows)
        {
            std::cout << "Two matricies don\'t have right shapes for multiplying" << std::endl;
            exit(1);
        }

        auto result = Matrix<T>(rows, other.columns);

        // #pragma omp parallel for reduction(+ : result)
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < other.columns; j++)
            {
                for (size_t k = 0; k < columns; k++)
                {
                    result[i, j] += (*this)[i, k] * other[k ,j];
                }
            }
        }

        return result;
    }

    Matrix<T> hadamard_product(const Matrix<T>& other)
    {
        if ((columns != other.columns) || (rows != other.rows))
        {
            std::cout << "Two matricies don\'t have right shapes for Hadamard product" << '\n';
            exit(1);
        }

        auto result = Matrix<T>(rows, columns);
        // #pragma omp parallel for simd
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < columns; j++)
            {
                result[i, j] = (*this)[i, j] * other[i, j];
            }
        }

        return result;
    }

    Matrix<T> hadamard_division(const Matrix<T>& other)
    {
        if ((columns != other.columns) || (rows != other.rows))
        {
            std::cout << "Two matricies don\'t have right shapes for Hadamard division" << '\n';
            exit(1);
        }

        auto result = Matrix<T>(rows, columns);
        // #pragma omp parallel for simd
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < columns; j++)
            {
                result[i, j] = (*this)[i, j] / other[i, j];
            }
        }

        return result;
    }

    static Matrix<T> hadamard_product(const Matrix<T>& A, const Matrix<T>& B)
    {
        return A.hadamard_product(B);
    }

    static Matrix<T> hadamard_division(const Matrix<T>& A, const Matrix<T>& B)
    {
        return A.hadamard_division(B);
    }

    void write_to_file(std::string filename)
    {
        auto file = std::ofstream(filename, std::ios::out | std::ios::app);
        
        if (!file.is_open())
        {
            std::cout << "Couldn\'t open the file" << '\n';
            return;
        }

        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < columns-1; j++)
            {
                file << (*this)[i, j] << ' ';
            }
            file << (*this)[i, columns-1] << '\n';
        }

        file.close();
    }
};

template <class T>
Matrix<T> operator* (const T& scalar, const Matrix<T>& mat)
{
    int rows = mat.get_rows();
    int columns = mat.get_columns();

    auto result = Matrix<T>(rows, columns);
    // #pragma omp parallel for simd
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            result[i, j] = mat[i, j] * scalar;
        }
    }

    return result;
}

template <class T>
std::ostream& operator<< (std::ostream& os, const Matrix<T>& mat)
{
    int rows = mat.get_rows();
    int columns = mat.get_columns();

    for (size_t i = 0; i < rows; i++)
    {
        if (rows == 1) {os << "( ";}
        else if (i == 0) {os << "/ ";}
        else if (i == rows - 1) {os << "\\ ";}
        else {os << "| ";}

        for (size_t j = 0; j < columns; j++)
        {
            os << mat[i, j] << ' ';
        }

        if (rows == 1) {os << ")" << std::endl;}
        else if (i == 0) {os << "\\" << std::endl;}
        else if (i == rows - 1) {os << "/" << std::endl;}
        else {os << "|" << std::endl;}
    }

    return os;
}