#pragma once

#include <cstddef>
#include <memory>

namespace mexutil {

template <typename T>
class Buffer1D {
   public:
    explicit Buffer1D(std::size_t size) : values_(new T[size]) {}

    T& operator[](std::size_t index) { return values_[index]; }
    operator T*() noexcept { return values_.get(); }

   private:
    std::unique_ptr<T[]> values_;
};

template <typename T, std::size_t Columns>
class Buffer2D {
   public:
    explicit Buffer2D(std::size_t rows) : values_(new T[rows][Columns]) {}

    T* operator[](std::size_t row) { return values_[row]; }

    using CArrayPointer = T (*)[Columns];
    operator CArrayPointer() noexcept { return values_.get(); }

   private:
    std::unique_ptr<T[][Columns]> values_;
};

template <typename T, std::size_t Rows, std::size_t Columns>
class Buffer3D {
   public:
    explicit Buffer3D(std::size_t depth)
        : values_(new T[depth][Rows][Columns]) {}

    using RowPointer = T (*)[Columns];
    RowPointer operator[](std::size_t depth) { return values_[depth]; }

    using CArrayPointer = T (*)[Rows][Columns];
    operator CArrayPointer() noexcept { return values_.get(); }

   private:
    std::unique_ptr<T[][Rows][Columns]> values_;
};

}  // namespace mexutil
