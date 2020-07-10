#ifndef __LWE_CONTAINER__
#define __LWE_CONTAINER__

#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

namespace LWE {
    template <typename T, uint64_t LENGTH> class Vector {
    public:
        std::array<T, LENGTH> vec;

        Vector() : vec() {}

        Vector(const Vector &other) = default;

        inline T &operator[](const uint64_t &i) { return this->vec[i]; }

        inline const T &operator[](const uint64_t &i) const {
            return this->vec[i];
        }

        inline Vector &operator=(const Vector &other) {
            for (uint64_t i = 0; i < LENGTH; i++)
                this->vec[i] = other[i];
            return *this;
        }

        inline Vector &operator+=(const Vector &other) {
            for (uint64_t i = 0; i < LENGTH; i++)
                this->vec[i] += other[i];
            return *this;
        }

        inline Vector operator+(const Vector &other) const {
            Vector what(*this);
            what += other;
            return what;
        }

        inline Vector &operator-=(const Vector &other) {
            for (uint64_t i = 0; i < LENGTH; i++)
                this->vec[i] -= other[i];
            return *this;
        }

        inline Vector operator-(const Vector &other) const {
            Vector what(*this);
            what -= other;
            return what;
        }

        inline Vector &operator*=(const T &var) {
            for (uint64_t i = 0; i < LENGTH; i++)
                this->vec[i] = var * this->vec[i];
            return *this;
        }

        inline Vector operator*(const T &var) const {
            Vector what(*this);
            what *= var;
            return what;
        }

        inline Vector &operator*=(const uint64_t &var) {
            this->operator*=(T(var));
            return *this;
        }

        inline Vector operator*(const uint64_t &var) const {
            Vector what(*this);
            what *= var;
            return what;
        }

        template <typename Tp> inline Vector &operator*=(const Tp &var) {
            for (auto &i : this->vec)
                i = var.lift_ring_multiply(i);
            return *this;
        }

        template <typename Tp> inline Vector operator*(const Tp &var) const {
            Vector what(*this);
            what *= var;
            return what;
        }

        inline Vector &rescale(__uint128_t modulus_scale, uint64_t p_prime) {
            for (auto &i : this->vec)
                i.rescale(modulus_scale, p_prime);
            return *this;
        }

        inline Vector &discrete_gaussian() {
            T::discrete_gaussian_sequence(this->vec);
            return *this;
        }

        inline Vector &bounded(const __uint128_t &bound) {
            T::bounded_sequence(bound, this->vec);
            return *this;
        }

        inline Vector &pm_bounded(const __uint128_t &bound) {
            T::pm_bounded_sequence(bound, this->vec);
            return *this;
        }

        inline Vector &random_element() {
            T::random_element_sequence(this->vec);
            return *this;
        }
    };

    template <typename T, uint64_t ROW, uint64_t COLUMN> class Matrix {
    public:
        std::vector<std::array<T, COLUMN>> mat;

        Matrix() : mat(ROW) {}

        Matrix(const Matrix &other) = default;

        inline std::array<T, COLUMN> &operator[](const uint64_t &i) {
            return this->mat[i];
        }

        inline const std::array<T, COLUMN> &
        operator[](const uint64_t &i) const {
            return this->mat[i];
        }

        inline Matrix<T, ROW, COLUMN> &operator*=(const T &var) {
            for (uint64_t i = 0; i < ROW; i++)
                for (uint64_t j = 0; j < COLUMN; j++)
                    this->mat[i][j] *= var;
            return *this;
        }

        inline Vector<T, ROW> operator*(const Vector<T, COLUMN> &v) const {
            Vector<T, ROW> res;
            for (uint64_t i = 0; i < ROW; i++)
                for (uint64_t j = 0; j < COLUMN; j++)
                    res[i] += (this->mat[i][j] * v[j]);
            return res;
        }

        template <uint64_t C2>
        inline Matrix<T, ROW, C2>
        operator*(const Matrix<T, COLUMN, C2> &o) const {
            Matrix<T, ROW, C2> res;
            for (uint64_t i = 0; i < ROW; i++)
                for (uint64_t j = 0; j < C2; j++)
                    for (uint64_t k = 0; k < COLUMN; k++)
                        res[i][j] += (this->mat[i][k] * o[k][j]);
            return res;
        }

        inline Matrix &operator+=(const Matrix &o) {
            for (uint64_t i = 0; i < ROW; i++)
                for (uint64_t j = 0; j < COLUMN; j++)
                    this->mat[i][j] += o[i][j];
            return *this;
        }

        inline Matrix operator+(const Matrix &o) const {
            Matrix what(*this);
            what += o;
            return what;
        }

        inline Matrix &discrete_gaussian() {
            for (uint64_t i = 0; i < ROW; i++)
                T::discrete_gaussian_sequence(this->mat[i]);
            return *this;
        }

        inline Matrix &bounded(const T &bound) {
            for (uint64_t i = 0; i < ROW; i++)
                T::bounded_sequence(bound, this->mat[i]);
            return *this;
        }

        inline Matrix &random_element() {
            for (uint64_t i = 0; i < ROW; i++)
                T::random_element_sequence(this->mat[i]);
            return *this;
        }
    };
}

#endif
