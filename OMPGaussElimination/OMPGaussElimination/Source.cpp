#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

#define SIZE 5000
//input matrix
//входна€ матрица
long double matrix[SIZE][SIZE], b[SIZE], x[SIZE], y[SIZE];
//block sizes
//размеры блоков дл€ тестов
int block_sizes[] = { 0,2,5,10,20,50,100 };
//start_time, end_time - start and end time of the experement
//врем€ выполнени€ эксперимента
//diff_time - time of an experement
double start_time, end_time, diff_time;
//output file
//выходной файл
ofstream fout;

//Return min value
//возвращает минимальное значение
int min(int x, int y) {
	return (x < y) ? x : y;
}

//Print given matrix
//вывод матрицы
//size - size of current experiment matrix
void print(int size) {
	cout << "\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << "\n";
	}
}

//Printing head for the experement
//печатает заголовок эксперимента
void print_head(int block_size = 0) {
	fout << (block_size == 0 ? "–азмер блока = размеру матрицы" : "–азмер блока равен ") << block_size << "\t\n";
	cout << (block_size == 0 ? "–азмер блока = размеру матрицы" : "–азмер блока равен ") << block_size << "\t\n";
	fout << "–азмер матрицы" << "\t" << "1 поток" << "\t" << omp_get_max_threads() << "потоков" << "\t\n";
	cout << "–азмер матрицы" << "\t" << "1 поток" << "\t" << omp_get_max_threads() << "потоков" << "\t\n";
}

//Generate experiment matrix
//генерирует случайную входную матрицу(одинакова дл€ каждого эксперимента)
//size - size of current experiment matrix
void generate_matrix(int size) {
	srand(size);
	long double diag_dominant = 0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i != j) {
				matrix[i][j] = rand() % 100 + 1;
				diag_dominant += abs(matrix[i][j]);
			}
		}
		matrix[i][i] = diag_dominant + 1;
		diag_dominant = 0;
		b[i] = rand() % 2;
	}
}

//Recursive parallel block lu decomposition
//рекурсивное блочное lu разложение
//matr - given matrix
//block_size - size of a block
//size - size of the matrix
//s - start index for decomposition
//nthreads - num of threads
void lu_decompose(long double matr[SIZE][SIZE], int block_size, int size, int s = 0, int nthreads = 1) {
	int i = 0, j = 0, k = 0, ii = 0, jj = 0, kk = 0;
	long double factor;
#pragma omp parallel num_threads(nthreads) private(i,j,k,ii,jj,kk,factor)
	{
		//Solving for L11 and U11
		//нахождение матриц L11 и U11
		for (j = s; j < block_size - 1 + s; j++) {
#pragma omp for
			for (i = j + 1; i < block_size + s; i++) {
				factor = matr[i][j] / matr[j][j];
				for (k = s; k < block_size + s; k++) {
					if (k < j)continue;
					matr[i][k] = matr[i][k] - (matr[j][k] * factor);
				}
				matr[i][j] = factor;
			}
		}
		//Solving for U12
		//нахождение матрицы U12
		for (j = block_size + s; j < size + s; j++) {
#pragma omp for
			for (i = 1 + s; i < block_size + s; i++) {
				for (k = s; k < i; k++) {
					matr[i][j] = matr[i][j] - matr[i][k] * matr[k][j];
				}
			}
		}
		//Solving for L12
		//нахождение матрицы L12
#pragma omp for 
		for (i = s + block_size; i < size + s; i++) {
			for (j = s; j < block_size + s; j++) {
				factor = matr[j][j];
				for (k = s; k < j; k++) {
					matr[i][j] = matr[i][j] - matr[i][k] * matr[k][j];
				}
				matr[i][j] /= factor;
			}
		}
		//Solving for reductive A22
		//нахождение редуцированной матрицы A22
#pragma omp for
		for (ii = block_size + s; ii < size + s; ii += block_size) {
			for (jj = block_size + s; jj < size + s; jj += block_size) {
				for (kk = s; kk < block_size + s; kk += block_size) {
					for (i = ii; i < min(ii + block_size, size + s); i++) {
						for (j = jj; j < min(jj + block_size, size + s); j++) {
							for (k = kk; k < min(kk + block_size, block_size + s); k++) {
								matr[i][j] -= matr[i][k] * matr[k][j];
							}
						}
					}
				}
			}
		}
	}
	//Recursive lu factorization for A22
	//рекурсивный вызов дл€ блока A22
	if (size - block_size > 0) {
		lu_decompose(matr, (block_size < size - block_size) ? block_size : size - block_size, size - block_size, block_size + s, nthreads);
	}
}

//Solve lu decomosed matrix 
//решение разложенной матрицы
//size - size of current experiment matrix
void solve(long double matr[SIZE][SIZE], int size) {
	int i = 0, k = 0;
	//Solving for Ly=b
	for (i = 0; i < size; i++) {
		for (k = 0; k < i; k++) {
			b[i] = b[i] - matr[i][k] * y[k];
		}
		y[i] = b[i];
	}
	//Solving for Ux=y
	for (i = size - 1; i >= 0; i--) {
		for (k = size - 1; k > i; k--) {
			y[i] = y[i] - matr[i][k] * x[k];
		}
		x[i] = y[i] / matr[i][i];
	}
}

//Perform experement 
//выполнение эксперемента
void perform_experement(int blck_size = 0) {
	//block_size - size of block for lu decomposition
	//k - step for matrix size
	int block_size = 0, k = 0;
	for (k = 1000; k <= 5000; k += 1000) {
		fout << k << "\t";
		cout << k << "\t";
		//in one thread
		//исполнение в один поток
		generate_matrix(k);
		block_size = (blck_size == 0) ? k : blck_size;
		start_time = omp_get_wtime();
		lu_decompose(matrix, block_size, k);
		solve(matrix, k);
		end_time = omp_get_wtime();
		diff_time = difftime(end_time, start_time);
		fout << diff_time << "\t";
		cout << diff_time << "\t";
		//in n threads
		//исполнение использую все потоки
		generate_matrix(k);
		start_time = omp_get_wtime();
		lu_decompose(matrix, block_size, k, 0, omp_get_max_threads());
		solve(matrix, k);
		end_time = omp_get_wtime();
		diff_time = difftime(end_time, start_time);
		fout << diff_time << "\t";
		cout << diff_time << "\t";
		//
		fout << endl;
		cout << endl;
	}
	fout << endl;
	cout << endl;
}

int main() {
	fout.open("results.xls");
	for (int i = 0; i < 6; i++) {
		print_head(block_sizes[i]);
		perform_experement(block_sizes[i]);
	}
	fout.close();
	return 0;
}