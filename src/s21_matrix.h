#ifndef S21_MATRIX_H
#define S21_MATRIX_H

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SUCCESS 1
#define FAILURE 0
#define K 10000000.0
#define SUCCESS 1
#define FAILURE 0
#define EPS 1e-6

typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_create_matrix(int rows, int columns, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
matrix_t s21_create_minor(matrix_t *mat, int row, int col);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
double s21_sum(matrix_t *A, matrix_t *B, int i, int j);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
#endif