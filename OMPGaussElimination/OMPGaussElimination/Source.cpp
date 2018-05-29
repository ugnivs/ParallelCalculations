#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

//input experiment configuration file
//конфигурация эксперимента
ifstream inconfig;
//output file
//выходной файл
ofstream fout;

//Print given matrix
//вывод матрицы
void print(double** matrix, double* b, int size) {
	cout << "\n";
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << "\t" << b[i] << " ";
		cout << "\n";
	}
}

//Printing head for the experement
//печатает заголовок эксперимента
void print_head(int matr_size) {
	fout << "Matrix size = " << matr_size << "\n";
	cout << "Matrix size = " << matr_size << "\n";
	fout << "Block size" << "\t" << "Time" << "\t\t" << "ln(Block size)" << "\t" << "ln(time)" << "\t\t\n";
	cout << "Block size" << std::setw(10) << "Time" << std::setw(20) << "ln(Block size)" << std::setw(15) << "ln(time)" << "\n";
}

//Alloc n,m matrix
//Выделение памяти под матрицу n,m
double** alloc_matrix(int n, int m) {
	double** result = new double*[n];
	for (int i = 0; i < n; i++)
		result[i] = new double[m];
	return result;
}

void free_matrix(double** matrix, int size) {
	for (int i = 0; i < size; i++)
		delete[] matrix[i];
	delete[] matrix;
}

//Generate random n,m matrix and vector b
//Генерация случайно матрицы n,m и вектора b
double** generate_matrix(int size, double* b) {
	double** matrix = alloc_matrix(size, size);
	srand(size);
	double diag_dominant = 0;
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
	return matrix;
}

int min(int x, int y) {
	return (x < y) ? x : y;
}

//LU decomposition of the matrix
//LU разложение матрицы
double** lu_decompose(double** matrix, int size, int block_size, int s = 0) {
	int i = 0, j = 0, k = 0, ii = 0, jj = 0, kk = 0;
	double factor;
	//Solving for L11 and U11
	//нахождение матриц L11 и U11
	for (j = s; j < block_size - 1 + s; j++) {
		for (i = j + 1; i < block_size + s; i++) {
			factor = matrix[i][j] / matrix[j][j];
			for (k = j; k < block_size + s; k++) {
				matrix[i][k] = matrix[i][k] - (matrix[j][k] * factor);
			}
			matrix[i][j] = factor;
		}
	}
	//Solving for U12
	//нахождение матрицы U12
	for (j = block_size + s; j < size + s; j++) {
		for (i = 1 + s; i < block_size + s; i++) {
			for (k = s; k < i; k++) {
				matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
			}
		}
	}
	//Solving for L12
	//нахождение матрицы L12
	for (i = s + block_size; i < size + s; i++) {
		for (j = s; j < block_size + s; j++) {
			factor = matrix[j][j];
			for (k = s; k < j; k++) {
				matrix[i][j] -= matrix[i][k] * matrix[k][j];
			}
			matrix[i][j] /= factor;
		}
	}
	//Solving for reductive A22
	//нахождение редуцированной матрицы A22
	for (ii = block_size + s; ii < size + s; ii += block_size)
		for (jj = block_size + s; jj < size + s; jj += block_size)
			for (kk = s; kk < block_size + s; kk += block_size)
				for (i = ii; i < min(ii + block_size, size + s); i++)
					for (j = jj; j < min(jj + block_size, size + s); j++)
						for (k = kk; k < min(kk + block_size, block_size + s); k++)
							matrix[i][j] -= matrix[i][k] * matrix[k][j];
	//Recursive lu factorization for A22
	//рекурсивный вызов для блока A22
	int next_size = size - block_size;
	if (next_size > 0) {
		return lu_decompose(matrix, next_size, (block_size < next_size) ? block_size : next_size, block_size + s);
	}
	return matrix;
}

//Solve lu decomosed matrix 
//решение разложенной матрицы
void solve(double** matrix, double* b, double*x, int size, int block_size) {
	int i = 0, k = 0;
	lu_decompose(matrix, size, block_size);
	double* y = new double[size];
	//Solving for Ly=b
	for (i = 0; i < size; i++) {
		for (k = 0; k < i; k++) {
			b[i] = b[i] - matrix[i][k] * y[k];
		}
		y[i] = b[i];
	}
	//Solving for Ux=y
	for (i = size - 1; i >= 0; i--) {
		for (k = size - 1; k > i; k--) {
			y[i] = y[i] - matrix[i][k] * x[k];
		}
		x[i] = y[i] / matrix[i][i];
	}
	delete[] y;
}

int main() {
	setlocale(LC_ALL, "Rus");
	inconfig.open("config.txt");
	fout.open("results.xls");
	int matr_size, start_block_size, small_step, middle_block_size, big_step, last_block_size;
	double start_watch, end_watch;
	double* b;
	double* x;
	double** matrix;
	inconfig >> matr_size >> start_block_size >> small_step >> middle_block_size >> big_step >> last_block_size;
	print_head(matr_size);
	//Block Gauss
	while (start_block_size <= last_block_size&&start_block_size != matr_size) {
		b = new double[matr_size];
		x = new double[matr_size];
		matrix = generate_matrix(matr_size, b);
		start_block_size = start_block_size == 0 ? matr_size : start_block_size;
		start_block_size = matr_size;
		start_watch = omp_get_wtime();
		solve(matrix, b, x, matr_size, start_block_size);
		end_watch = omp_get_wtime() - start_watch;
		fout << start_block_size << "\t\t" << end_watch << "\t\t" << std::log(start_block_size) << "\t\t" << std::log(end_watch) << "\t\t\n";
		cout << start_block_size << std::setw(20) << end_watch << std::setw(15) << std::log(start_block_size) << std::setw(20) << std::log(end_watch) << "\n";
		if (start_block_size == last_block_size)
			start_block_size = 0;
		else
			start_block_size += start_block_size < middle_block_size ? small_step : big_step;
		free_matrix(matrix, matr_size);
		delete[] b;
		delete[] x;
	}
	inconfig.close();
	fout.close();
	system("pause");
	return 0;
}