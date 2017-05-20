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
long double matrix[SIZE][SIZE], b[SIZE], currentValues[SIZE], previousValues[SIZE];
//block sizes
int block_sizes[] = { 0,2,5,10,20,50,100 };
//definite
long double eps;
//start_time, end_time - start and end time of the experement
//diff_time - time of an experement
double start_time, end_time, diff_time;
//output file
ofstream fout;

//Print given matrix
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
void print_head(int block_size = 0) {
	fout << (block_size == 0 ? "Размер блока = размеру матрицы" : "Размер блока равен ") << block_size << "\t\n";
	cout << (block_size == 0 ? "Размер блока = размеру матрицы" : "Размер блока равен ") << block_size << "\t\n";
	fout << "Размер матрицы" << "\t" << "1 поток" << "\t" << omp_get_max_threads() << "потоков" << "\t\n";
	cout << "Размер матрицы" << "\t" << "1 поток" << "\t" << omp_get_max_threads() << "потоков" << "\t\n";
}

//Return min value
int min(int x, int y) {
	return (x < y) ? x : y;
}

//Generate experiment matrix
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
		previousValues[i] = 0;
	}
}

void solve(long double matr[SIZE][SIZE], int size, int block_size, int nthreads = 1) {
	while (true)
	{
		int i = 0, j = 0, ii = 0, jj = 0;
#pragma omp parallel for num_threads(nthreads)
		for (ii = 0; ii < size; ii += block_size) {
			for (jj = 0; jj < size; jj += block_size) {
				for (i = ii; i < min(ii + block_size, size); i++) {
					if (jj == 0)
						currentValues[i] = b[i];
					for (j = jj; j < min(jj + block_size, size); j++) {
						if (i != j)
							currentValues[i] -= matr[i][j] * previousValues[j];
					}
					if (jj + block_size == size)
						currentValues[i] /= matr[i][i];
				}
			}
		}

		long double error = 0.0;

		for (int i = 0; i < size; i++)
		{
			error += abs(currentValues[i] - previousValues[i]);
		}

		if (error < eps)
		{
			break;
		}

		for (int i = 0; i < size; i++) {
			previousValues[i] = currentValues[i];
		}
	}
}

//Perform experement 
void perform_experement(int blck_size = 0) {
	//block_size - size of block for lu decomposition
	//k - step for matrix size
	int block_size = 0, k = 0;
	eps = pow(10, -6);
	for (k = 100; k <= 500; k += 100) {
		fout << k << "\t";
		cout << k << "\t";
		//in one threads
		generate_matrix(k);
		block_size = (blck_size == 0) ? k : blck_size;
		start_time = omp_get_wtime();
		solve(matrix, k, block_size);
		end_time = omp_get_wtime();
		diff_time = difftime(end_time, start_time);
		fout << diff_time << "\t";
		cout << diff_time << "\t";
		//in n threads
		generate_matrix(k);
		start_time = omp_get_wtime();
		solve(matrix, k, block_size, omp_get_max_threads());
		end_time = omp_get_wtime();
		diff_time = difftime(end_time, start_time);
		fout << diff_time << "\t";
		cout << diff_time << "\t";
		//
		fout << endl;
		cout << endl;
	}
}

int main()
{
	fout.open("results.xls");
	for (int i = 0; i < 6; i++) {
		print_head(block_sizes[i]);
		perform_experement(block_sizes[i]);
	}
	fout.close();
	return 0;
}