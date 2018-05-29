#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

#define SIZE 5000
//input matricurrentValues
//входная матрица
double matrix[SIZE][SIZE], b[SIZE];
//solution matrix
//матрица решения
double solution[SIZE][SIZE], b_s[SIZE], currentValues_s[SIZE], previousValues_s[SIZE];
//start_time, end_time - start and end time of the experement
//время выполнения эксперимент
//diff_time - time of an experement
double start_time, end_time, diff_time;
//output file
ofstream fout;

//Print given matrix
//вывод матрицы
//size - size of current experiment matrix
void print(double matr[SIZE][SIZE]) {
	cout << "\n";
	for (int i = 0; i < SIZE; i++)
	{
		for (int j = 0; j < SIZE; j++)
		{
			cout << matr[i][j] << " ";
		}
		cout << "\n";
	}
	for (int j = 0; j < SIZE; j++)
		cout << b[j] << " ";
}

//Printing head for the experement
//печатает заголовок эксперимента
void print_head(int block_size = 0) {
	fout << "Block size" << "\t" << "time" << "\t" << "ln(Block size)" << "\t" << "ln(time)" << "\t";
	cout << "Block size" << "\t" << "time" << "\t" << "ln(Block size)" << "\t" << "ln(time)" << "\t";
}

//Return min value
//возвращает минимальное значение
int min(int x, int y) {
	return (x < y) ? x : y;
}

//Generate experiment matrix
//генерирует случайную входную матрицу(одинакова для каждого эксперимента)
//size - size of current experiment matrix
void generate_matrix(int size) {
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
}

//Reset experiment
void reset(int size) {
	for (int i = 0; i < size; i++) {
		b_s[i] = b[i];
		previousValues_s[i] = 0;
		currentValues_s[i] = 0;
		for (int j = 0; j < size; j++) {
			solution[i][j] = matrix[i][j];
		}
	}
}

//solving
//решение методом просты итераций
void solve(double matr[SIZE][SIZE], int size, int block_size) {
	int iter = 0;
	while (iter <= 100)
	{
		int i = 0, j = 0, ii = 0, jj = 0;
		for (ii = 0; ii < size; ii += block_size) {
			for (jj = 0; jj < size; jj += block_size) {
				for (i = ii; i < min(ii + block_size, size); i++) {
					if (jj == 0)
						currentValues_s[i] = b[i];
					for (j = jj; j < min(jj + block_size, size); j++) {
						if (i != j)
							currentValues_s[i] -= matr[i][j] * previousValues_s[j];
					}
					if (jj + block_size == size)
						currentValues_s[i] /= matr[i][i];
				}
			}
		}

		for (int i = 0; i < size; i++) {
			previousValues_s[i] = currentValues_s[i];
		}
		iter++;
	}
}

//Perform experement 
void perform_experement(int matr_size = 0) {
	generate_matrix(matr_size);
	fout << matr_size << "\t\n";
	cout << matr_size << "\t\n";
	//block_size - size of block for lu decomposition
	//k - step for matrix size
	int block_size = 10;
	//in one threads
	//исполнение одним потоком
	while (block_size <= 2000) {
		reset(matr_size);
		start_time = omp_get_wtime();
		solve(solution, matr_size, block_size);
		end_time = omp_get_wtime();
		diff_time = difftime(end_time, start_time);
		fout << block_size << "\t" << diff_time << "\t\t" << log(block_size) << "\t" << log(diff_time) << "\t";
		cout << block_size << "\t" << diff_time << "\t\t" << log(block_size) << "\t" << log(diff_time) << "\t";
		fout << endl;
		cout << endl;
		if (block_size == 2000) {
			block_size = matr_size;
			continue;
		}
		block_size += (block_size < 50) ? 10 : 50;
	}
	fout << endl;
	cout << endl;
}

int main()
{
	fout.open("results.xls");
	print_head(SIZE);
	perform_experement(SIZE);
	fout.close();
	return 0;
}