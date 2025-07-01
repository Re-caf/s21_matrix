# Project: Matrix Operations Library (s21_matrix)

## Part 1: Matrix Operations Implementation
- **Task**: Implement matrix operations library with prefix `s21_`
- **Requirements**:
  - C11 standard, GCC compiler
  - Static library: `s21_matrix.a` (header: `s21_matrix.h`)
  - POSIX.1-2017 compliance, Google C++ Style
  - 80%+ test coverage via Check library
  - Makefile targets: `all`, `clean`, `test`, `gcov_report`

## Core Functionality
- **Matrix Operations**:
  - `create_matrix` / `remove_matrix` (memory management)
  - `eq_matrix` (comparison with 1e-6 precision)
  - Arithmetic: `sum_matrix`, `sub_matrix`, `mult_matrix`, `mult_number`
  - Advanced: `transpose`, `determinant`, `calc_complements`, `inverse_matrix`

## Implementation Details
- **Structure**:
  ```c
  typedef struct matrix_struct {
      double** matrix;
      int rows;
      int columns;
  } matrix_t;
