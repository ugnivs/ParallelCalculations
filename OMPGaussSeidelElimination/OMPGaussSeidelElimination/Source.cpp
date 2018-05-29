#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
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
//LU ���������� �������
double** lu_decompose(double** matrix, int size, int block_size, int s = 0) {
	int i = 0, j = 0, k = 0, ii = 0, jj = 0, kk = 0;
	double factor;
	//Solving for L11 and U11
	//���������� ������ L11 � U11
	for (j = s; j < block_size - 1 + s; j++) {
		for (i = j + 1; i < block_size + s; i++) {
			factor = matrix[i][j] / matrix[j][j];
			for (k = j; k < block_size + s; k++) {
				matrix[i][k] -= (matrix[j][k] * factor);
			}
			matrix[i][j] = factor;
		}
	}
	//Solving for U12
	//���������� ������� U12
	for (j = block_size + s; j < size + s; j++) {
		for (i = 1 + s; i < block_size + s; i++) {
			for (k = s; k < i; k++) {
				matrix[i][j] -= matrix[i][k] * matrix[k][j];
			}
		}
	}
	//Solving for L12
	//���������� ������� L12
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
	//���������� �������������� ������� A22
	for (ii = block_size + s; ii < size + s; ii += block_size)
		for (jj = block_size + s; jj < size + s; jj += block_size)
			for (kk = s; kk < block_size + s; kk += block_size)
				for (i = ii; i < min(ii + block_size, size + s); i++)
					for (j = jj; j < min(jj + block_size, size + s); j++)
						for (k = kk; k < min(kk + block_size, block_size + s); k++)
							matrix[i][j] -= matrix[i][k] * matrix[k][j];
	//Recursive lu factorization for A22
	//����������� ����� ��� ����� A22
	int next_size = size - block_size;
	if (next_size > 0) {
		return lu_decompose(matrix, next_size, (block_size < next_size) ? block_size : next_size, block_size + s);
	}
	return matrix;
}

//Inversing matrix
//Вычисление обратной матрицы
double** inverse(double** matrix, int size, int block_size) {
	int i, j, k;
	double s;

	lu_decompose(matrix, size, block_size);

	//Inversing U
	//Нахождение обратной матрицы U
	for (i = size - 1; i >= 0; i--) {
		for (j = size - 1; j >= i; j--) {
			if (i == j) {
				for (k = size - 1; k > j; k--)
					matrix[i][k] /= matrix[i][j];
				matrix[i][j] = 1 / matrix[i][j];
			}
			else {
				for (k = size - 1; k > j; k--)
					matrix[i][k] -= matrix[i][j] * matrix[j][k];
				matrix[i][j] = -matrix[i][j] * matrix[j][j];
			}
		}
	}
	//Inversing L
	//Нахождение обратной матрицы L
	for (i = 1; i < size; i++) {
		for (j = 0; j < i; j++) {
			for (k = 0; k < j; k++)
				matrix[i][k] -= matrix[i][j] * matrix[j][k];
			matrix[i][j] = -matrix[i][j];
		}
	}
	//Multiplying U inversion with L inversion
	//Нахождение обратной матрицы исходной, путем перемножения обратных матриц U и L
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			s = 0;
			for (k = std::fmax(i, j); k < size; k++) {
				s += (j != k) ? matrix[i][k] * matrix[k][j] : matrix[i][j];
			}
			matrix[i][j] = s;
		}
	}
	return matrix;
}

//Array of diagonal inversions
//Массив обратных диагональных матриц
double*** d_inversions(double** matrix, int size, int block_size) {
	int i, j, k, num_blocks = size / block_size;
	double*** dInv = new double**[num_blocks];
	for (k = 0; k < num_blocks; k++) {
		dInv[k] = alloc_matrix(block_size, block_size);
		int block_start = block_size*k;
		int block_end = block_size*(k + 1);
		for (i = block_start; i < block_end; i++)
			for (j = block_start; j < block_end; j++)
				dInv[k][i - block_start][j - block_start] = matrix[i][j];
		inverse(dInv[k], block_size, block_size);
	}
	return dInv;
}

//Block Seidel elimination
//Блочный метод Зейделя
void solve_block(double** matrix, double* b, int size, int block_size, int K) {
	int i, j, i_b, j_b, k = 0, num_blocks = size / block_size;
	double*** dInv = d_inversions(matrix, size, block_size);
	double* mult = new double[size];
	double* x = new double[size];
	double* x0 = new double[size];

	for (i = 0; i < size; i++)
		x0[i] = x[i] = 0;

	while (k < K) {
		for (i = 0; i < size; i++) {
			mult[i] = b[i];
			x0[i] = x[i];
			x[i] = 0;
		}
		for (i_b = 0; i_b < size; i_b += block_size) {
			int i_block_end = i_b + block_size;
			int j_block_end = 0;

			for (j_b = i_b + block_size; j_b < size; j_b += block_size) {
				j_block_end = j_b + block_size;
				for (i = i_b; i < i_block_end; i++)
					for (j = j_b; j < j_block_end; j++)
						mult[i] -= matrix[i][j] * x0[j];
			}

			for (j_b = 0; j_b < i_b; j_b += block_size) {
				j_block_end = j_b + block_size;
				for (i = i_b; i < i_block_end; i++)
					for (j = j_b; j < j_block_end; j++)
						mult[i] -= matrix[i][j] * x[j];
			}

			int block_num = i_b / block_size;
			for (i = i_b; i < i_block_end; i++)
				for (j = i_b; j < i_block_end; j++)
					x[i] += dInv[block_num][i - i_b][j - i_b] * mult[j];
		}
		k++;
	}
	delete[]mult;
	delete[]x;
	delete[]x0;
	for (i = 0; i < num_blocks; i++)
		free_matrix(dInv[i], block_size);
	delete[]dInv;
}

//Classic Seidel elimination
//Точечный метод Зейделя
void solve(double** matrix, double* b, int size, int K) {
	int i, j, k = 0;
	double mult, dInv;
	double* x = new double[size];
	double* x0 = new double[size];

	for (i = 0; i < size; i++)
		x0[i] = x[i] = 0;

	while (k < K)
	{
		for (i = 0; i < size; i++) {
			mult = 0;
			dInv = 1 / matrix[i][i];
			for (j = 0; j < i; j++) {
				mult += matrix[i][j] * x[j];
			}
			x[i] -= dInv*mult;
			mult = 0;
			for (j = i + 1; j < size; j++) {
				mult += matrix[i][j] * x0[j];
			}
			x[i] -= dInv*(mult - b[i]);
		}
		for (i = 0; i < size; i++) {
			x0[i] = x[i];
			x[i] = 0;
			mult = 0;
		}
		k++;
	}
	delete[]x;
	delete[]x0;
}

int main() {
	setlocale(LC_ALL, "Rus");
	inconfig.open("config.txt");
	fout.open("results.xls");
	int matr_size, iter_count, start_block_size, small_step, middle_block_size, big_step, last_block_size;
	double start_watch, end_watch;
	inconfig >> matr_size >> iter_count >> start_block_size >> small_step >> middle_block_size >> big_step >> last_block_size;
	print_head(matr_size);
	//Block Seidel
	double* b = new double[matr_size];
	double** matrix = generate_matrix(matr_size, b);
	while (start_block_size <= last_block_size) {
		if (matr_size%start_block_size == 0) {
			start_watch = clock();
			solve_block(matrix, b, matr_size, start_block_size, iter_count);
			end_watch = (clock() - start_watch) / 1000;
			fout << start_block_size << "\t\t" << end_watch << "\t\t" << std::log(start_block_size) << "\t\t" << std::log(end_watch) << "\t\t\n";
			cout << start_block_size << std::setw(20) << end_watch << std::setw(15) << std::log(start_block_size) << std::setw(20) << std::log(end_watch) << "\n";
		}
		start_block_size += start_block_size < middle_block_size ? small_step : big_step;
	}
	//Classic
	start_watch = clock();
	solve(matrix, b, matr_size, iter_count);
	end_watch = (clock() - start_watch) / 1000;
	fout << start_block_size << "\t\t" << end_watch << "\t\t" << std::log(start_block_size) << "\t\t" << std::log(end_watch) << "\t\t\n";
	cout << 1 << std::setw(20) << end_watch << std::setw(15) << std::log(1) << std::setw(20) << std::log(end_watch) << "\n";
	//
	free_matrix(matrix, matr_size);
	delete[] b;
	inconfig.close();
	fout.close();
	system("pause");
	return 0;
}