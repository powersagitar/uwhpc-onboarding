#pragma once

#include <cstddef>
#include <tuple>
#include <vector>

// Starter Grid for the 2D heat-diffusion problem.
//
// The evaluation harness uses operator() to set initial conditions and to read
// results; it never touches your internal storage. Keep this interface,
// everything else is yours.
class Grid {
private:
  std::size_t rows_;
  std::size_t cols_;
  std::vector<double> grid_;

public:
  Grid(std::size_t rows, std::size_t cols)
      : rows_(rows), cols_(cols), grid_(rows * cols, 0) {}

  double &operator()(std::size_t i, std::size_t j) {
    return grid_[i * cols_ + j];
  }

  double operator()(std::size_t i, std::size_t j) const {
    return grid_[i * cols_ + j];
  }

  std::tuple<std::size_t, std::size_t> dims() const { return {rows_, cols_}; }
};

void apply_stencil_on_boundaries(const Grid &old_grid, Grid &new_grid) {
  const auto [rows, cols] = old_grid.dims();

  // top & bottom
  for (size_t col = 0; col < cols; ++col) {
    new_grid(0, col) = old_grid(0, col);
    new_grid(rows - 1, col) = old_grid(rows - 1, col);
  }

  // left & right
  for (size_t row = 0; row < rows; ++row) {
    new_grid(row, 0) = old_grid(row, 0);
    new_grid(row, cols - 1) = old_grid(row, cols - 1);
  }
}

void apply_stencil_on_interior(const Grid &old_grid, Grid &new_grid) {
  const auto [rows, cols] = old_grid.dims();

  for (std::size_t row = 1; row < rows - 1; ++row) {
    for (std::size_t col = 1; col < cols - 1; ++col) {
      new_grid(row, col) =
          0.5 * old_grid(row, col) +
          0.125 * (old_grid(row - 1, col) + old_grid(row + 1, col) +
                   old_grid(row, col - 1) + old_grid(row, col + 1));
    }
  }
}

// Apply the five-point stencil over all interior points, copying the boundary
// values unchanged from old_grid to new_grid. Implement your solution here.
void apply_stencil(const Grid &old_grid, Grid &new_grid) {
  apply_stencil_on_boundaries(old_grid, new_grid);
  apply_stencil_on_interior(old_grid, new_grid);
}
