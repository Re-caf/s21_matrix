#include "s21_matrix.h"

matrix_t s21_create_minor(matrix_t *mat, int row, int col) {
  matrix_t minor;
  s21_create_matrix(mat->rows - 1, mat->columns - 1, &minor);
  int m = 0;  // индекс для строк новой матрицы
  for (int i = 0; i < mat->rows; i++) {
    if (i == row) continue;  // пропускаем строку
    int n = 0;  // индекс для столбцов новой матрицы
    for (int j = 0; j < mat->columns; j++) {
      if (j == col) continue;  // пропускаем столбец
      minor.matrix[m][n] = mat->matrix[i][j];
      n++;
    }
    m++;
  }

  return minor;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error = 0;

  if (A == NULL || A->matrix == NULL || A->columns <= 0 || A->rows <= 0)
    error = 1;

  if (!error && A->rows != A->columns) error = 2;

  if (!error) s21_create_matrix(A->rows, A->columns, result);

  for (int i = 0; !error && i < A->rows; i++) {
    for (int j = 0; j < A->columns && !error; j++) {
      matrix_t minor = s21_create_minor(A, i, j);
      double minor_det;
      if (s21_determinant(&minor, &minor_det) != 0) {
        s21_remove_matrix(&minor);
        error = 1;
      }

      result->matrix[i][j] = (i + j) % 2 == 0 ? minor_det : -minor_det;
      s21_remove_matrix(&minor);
    }
  }

  return error;
}

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = 0;
  if (rows <= 0 || columns <= 0) error = 1;
  if (!error) {
    result->rows = rows;
    result->columns = columns;

    result->matrix = malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++)
      result->matrix[i] = malloc(columns * sizeof(double));

    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = 0.0;
      }
    }
  }
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i] != NULL) {
        free(A->matrix[i]);
      }
    }
    free(A->matrix);
    A->matrix = NULL;  //пооочемуу
  }
}

int s21_determinant(matrix_t *A, double *result) {
  int error = 0;
  int finish = 0;

  if (A->rows != A->columns || A->columns <= 0 || A->rows <= 0 || A == NULL ||
      A->matrix == NULL)
    error = 1;

  if (A->rows != A->columns) error = 2;

  if (!error && A->rows == 1) {
    *result = A->matrix[0][0];
    finish = 1;
  }

  if (!error && finish != 1 && A->rows == 2) {
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    finish = 1;
  }
  if (!error && finish != 1) *result = 0;

  for (int i = 0; i < A->rows && !error && finish != 1; i++) {
    matrix_t minor = s21_create_minor(A, 0, i);
    if (minor.matrix == NULL) {
      s21_remove_matrix(&minor);
      error = 1;
    }
    double minor_det;
    if (!error && s21_determinant(&minor, &minor_det) != 0) {
      s21_remove_matrix(&minor);
      error = 1;
    }
    if (!error)
      *result += ((i % 2 == 0) ? 1 : -1) * A->matrix[0][i] * minor_det;
    s21_remove_matrix(&minor);
  }

  return error;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int result = SUCCESS;
  if (A->rows > 0 && B->rows > 0 && A->columns > 0 && B->columns > 0 &&
      A->rows == B->rows && A->columns == B->columns) {
    for (int row = 0; row < A->rows; row++) {
      for (int col = 0; col < A->columns; col++) {
        if (fabs(A->matrix[row][col] - B->matrix[row][col]) >= 1E-7) {
          result = FAILURE;
        }
      }
    }
  } else {
    result = FAILURE;
  }
  return result;
}

double s21_sum(matrix_t *A, matrix_t *B, int i, int j) {
  double sum = 0.0;
  for (int k = 0; k < A->columns; k++) sum += A->matrix[i][k] * B->matrix[k][j];
  return sum;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = 0;

  if (A->columns != B->rows || A->rows != B->columns)
    error = 2;
  else if (A == NULL || A->matrix == NULL || A->columns <= 0 || A->rows <= 0)
    error = 1;
  else if (B->columns <= 0 || B->rows <= 0)
    error = 1;
  else if (!error)
    s21_create_matrix(A->rows, B->columns, result);
  for (int i = 0; i < A->rows && !error; i++) {
    for (int j = 0; j < B->columns; j++) {
      result->matrix[i][j] = s21_sum(A, B, i, j);
    }
  }

  return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = 0;

  if (A == NULL || A->matrix == NULL || A->columns <= 0 || A->rows <= 0)
    error = 1;
  else if (A->rows != A->columns)
    error = 2;
  if (!error) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return error;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = 0;

  if (A == NULL || A->matrix == NULL || B == NULL || B->matrix == NULL)
    error = 1;

  else if (B->columns <= 0 || B->rows <= 0 || A->columns <= 0 || A->rows <= 0)
    error = 1;
  else if (A->rows != B->rows || A->columns != B->columns)
    error = 2;

  else if (!error) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return error;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = 0;
  if (A == NULL || A->matrix == NULL || B == NULL || B->matrix == NULL)
    error = 1;

  else if (B->columns <= 0 || B->rows <= 0 || A->columns <= 0 || A->rows <= 0)
    error = 1;
  else if (A->rows != B->rows || A->columns != B->columns)
    error = 2;

  else if (!error) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }
  return error;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error = 0;

  if (A == NULL || A->matrix == NULL)
    error = 1;

  else if (A->columns <= 0 || A->rows <= 0)
    error = 1;

  else if (A->rows <= 0 || A->columns <= 0)
    error = 1;
  else if (!error) {
    s21_create_matrix(A->columns, A->rows, result);

    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
    return 0;
  }
  return error;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = 0;
  double determin = 0.0;

  if (A == NULL || A->matrix == NULL || A->columns <= 0 || A->rows <= 0)
    error = 1;

  else if (A->rows != A->columns)
    error = 2;

  else if (!s21_determinant(A, &determin) && determin) {
    matrix_t matrix1 = {NULL, 0, 0};
    matrix_t matrix2 = {NULL, 0, 0};

    s21_calc_complements(A, &matrix1);
    s21_transpose(&matrix1, &matrix2);
    s21_mult_number(&matrix2, 1 / determin, result);
    s21_remove_matrix(&matrix1);
    s21_remove_matrix(&matrix2);
  } else
    error = 2;

  return error;
}