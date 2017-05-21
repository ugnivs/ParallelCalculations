#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <time.h>
#include <fstream>

#include <Windows.h>
#undef max

using std::vector;
std::fstream fout;

// Класс "матрица", представляющий двумерный массив
// вещественных чисел. Размеры матрицы задаются при
// создании и впоследствии не меняются. Элементы
// инициализируются нулями.

class Matrix {
public:
	// Конструкторы
	Matrix()
		: n(0)
		, m(0)
	{}
	Matrix(int n, int m)
		: n(n)
		, m(m)
		, data(vector<double>(n * m))
	{}

	// Доступ к элементу по двум индексам
	const double& operator ()(int i, int j) const {
		return data[i * m + j];
	}
	double& operator ()(int i, int j) {
		return data[i * m + j];
	}

	// Вывод в поток
	friend inline std::ostream& operator <<(std::ostream& out, const Matrix& matrix) {
		for (int i = 0; i < matrix.n; ++i) {
			for (int j = 0; j < matrix.m; ++j) {
				out << matrix(i, j) << ' ';
			}
			out << '\n';
		}
		return out;
	}

	// Операция нахождения обратной матрицы.
	// Используется самая простая реализация метода Гаусса-Жордана
	// без выбора главного элемента. Вырожденность матрицы никак не
	// не обрабатывается. Функция применима только к квадратным матрицам.
	Matrix invert() const {
		assert(m == n);
		Matrix tmp(n, 2 * n);

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				tmp(i, j) = (*this)(i, j);
			}
			tmp(i, n + i) = 1.0;
		}

		for (int k = 0; k < n; ++k) {
			const double inverse = 1.0 / tmp(k, k);
			for (int j = k; j < n + n; ++j) {
				tmp(k, j) *= inverse;
			}

			for (int i = 0; i < n; ++i) {
				if (i != k) {
					const double c = tmp(i, k);
					for (int j = k; j < n + n; ++j) {
						tmp(i, j) -= c * tmp(k, j);
					}
				}
			}
		}
		return tmp.submatrix(0, n, n, n);
	}


	// Функиция возвращает часть матрицы с верхним левым элементом
	// (i, j) размера subN x subM
	Matrix submatrix(int i, int j, int subN, int subM) const {
		Matrix sub(subN, subM);
		for (int ii = 0; ii < subN; ++ii) {
			for (int jj = 0; jj < subM; ++jj) {
				sub(ii, jj) = (*this)(i + ii, j + jj);
			}
		}
		return sub;
	}

	// Функция вычисляет произведение матриц a и b и записывает результат в
	// матрицу c, которая должна иметь необходимый размер.
	friend inline void Product(Matrix& c, const Matrix& a, const Matrix& b) {
		assert(a.m == b.n);
		assert(c.GetN() == a.n && c.GetM() == b.m);
		for (int i = 0; i < a.n; ++i) {
			for (int j = 0; j < b.m; ++j) {
				c(i, j) = 0;
				for (int k = 0; k < a.m; ++k) {
					c(i, j) += a(i, k) * b(k, j);
				}
			}
		}
	}

	// По возможности следует применять явно функцию Product, чтобы
	// избегать лишних выделений памяти на хранение промежуточного результата.
	// Однако оператор умножения в общем случае может быть удобнее.
	friend inline Matrix operator *(const Matrix& a, const Matrix& b) {
		Matrix c(a.n, b.m);
		Product(c, a, b);
		return c;
	}

	// Декремент
	Matrix& operator -=(const Matrix& obj) {
		assert(n == obj.n && m == obj.m);
		for (size_t i = 0; i < data.size(); ++i) {
			data[i] -= obj.data[i];
		}
		return *this;
	}

	// Функции возвращают размеры матрицы
	int GetN() const {
		return n;
	}
	int GetM() const {
		return m;
	}

private:
	int n, m;
	vector<double> data;
};

// Для измерения времени
class Timer {
private:
	double start_time;
	double end_time;
public:
	Timer(double current = omp_get_wtime())
		: start_time(current)
		, end_time(current)
	{}
	void LogTotal() {
		end_time = omp_get_wtime();
		const double diff_time = difftime(end_time, start_time);
		start_time = end_time;
		fout << diff_time << "\t";
		std::cout << diff_time << "\t";
	}
};

class Solver {
private:
	const int N; // размер исходной матрицы до блочного разбиения
	const int r; // размер блока (обязан быть делителем N)
	const int n; // размер блочной матрицы (в блоках), то есть N / r

	vector< vector<Matrix> > aBlocks; // блоки, организованные в виде двумерного массива
	vector<Matrix> aBlocksInverted;   // обратные к диагональным блокам
	vector<Matrix> bBlocks;           // блочное разбиение вектора (правой части)


public:
	Solver(int N, int r, const Matrix& a, const Matrix& b)
		: N(N)
		, r(r)
		, n(N / r)
		, aBlocks(vector< vector<Matrix> >(n, vector<Matrix>(n))) // матрица из Matrix(r, r)
		, aBlocksInverted(vector<Matrix>(n))                      // массив из Matrix(r, r)
		, bBlocks(vector<Matrix>(n))                              // массив из Matrix(r, 1)
	{
		assert((a.GetN() == N) && (a.GetM() == N) &&
			(b.GetN() == N) && (b.GetM() == 1) && (N % r == 0));

		// Осуществляем разбиение на блоки.
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				aBlocks[i][j] = a.submatrix(i * r, j * r, r, r);
			}
		}

		for (int i = 0; i < n; ++i) {
			bBlocks[i] = b.submatrix(i * r, 0, r, 1);
		}
	}

	// Обычный последовательный метод Зейделя.
	Matrix SerialSeidel() {
		Timer timer;
		for (int i = 0; i < n; ++i) {
			aBlocksInverted[i] = aBlocks[i][i].invert();
		}

		vector<Matrix> xBlocks = bBlocks;
		vector<Matrix> xBlocksPrev(n);

		Matrix tmp(r, 1);

		while (!converge(xBlocks, xBlocksPrev, pow(10, -6))) {
			xBlocksPrev.swap(xBlocks);
			for (int i = 0; i < n; ++i) {
				xBlocks[i] = bBlocks[i];
				for (int j = i + 1; j < n; ++j) {
					Product(tmp, aBlocks[i][j], xBlocksPrev[j]);
					xBlocks[i] -= tmp;
					// Эквивалентно записи "xBlocks[i] -= aBlocks[i][j] * xBlocksPrev[j];",
					// но быстрее работает.
				}
			}

			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < i; ++j) {
					Product(tmp, aBlocks[i][j], xBlocks[j]);
					xBlocks[i] -= tmp;
					// "xBlocks[i] -= aBlocks[i][j] * xBlocks[j];"
				}
				xBlocks[i] = aBlocksInverted[i] * xBlocks[i];
			}
		}

		const Matrix& result = merge(xBlocks);
		timer.LogTotal();
		return result;
	}

	Matrix ParallelSeidel(int threadCount) {

		Timer timer;
		// Параллельно ищем обратные матрицы.
#pragma omp parallel for num_threads(threadCount)
		for (int i = 0; i < n; ++i) {
			aBlocksInverted[i] = aBlocks[i][i].invert();
		}

		// К этому моменту созданные OpenMP вспомогательные потоки завершатся,
		// далее вновь программа работает в один поток.

		vector<Matrix> xBlocks = bBlocks;
		vector<Matrix> xBlocksPrev(n, Matrix(r, 1));

		Matrix tmp(r, 1);

		while (!converge(xBlocks, xBlocksPrev, pow(10, -6))) {
			xBlocksPrev.swap(xBlocks);

#pragma omp parallel for firstprivate(tmp) num_threads(threadCount)
			// * firstprivate(tmp) означает, что каждый поток получит свою копию
			// исходной переменной tmp. Если указать private(tmp), то у каждого потока
			// будет своя tmp, созданная конструктором по умолчанию.
			// * schedule(dynamic, 1) означает, что работа между потоками будет распределяться
			// динамически: как только поток закончит вычисление для своего i,
			// он запросит новое i из тех, что ещё не отданы другим потокам.
			for (int i = 0; i < n; ++i) {
				xBlocks[i] = bBlocks[i];
				for (int j = i + 1; j < n; ++j) {
					Product(tmp, aBlocks[i][j], xBlocksPrev[j]);
					xBlocks[i] -= tmp;
				}
			}

			for (int t = 1; t < 2 * n; ++t) {
#pragma omp parallel for firstprivate(tmp) num_threads(threadCount)
				for (int j = fmax(0, t - n); j < t / 2; ++j) {
					const int i = t - j - 1;
					Product(tmp, aBlocks[i][j], xBlocks[j]);
					xBlocks[i] -= tmp;
				}
				if (t % 2 == 1) {
					const int i = t / 2;
					xBlocks[i] = aBlocksInverted[i] * xBlocks[i];
				}
			}
		}
		const Matrix& result = merge(xBlocks);
		timer.LogTotal();
		return result;
	}

private:
	// Функция, которая из блоков склеивает столбец исходного размера.
	Matrix merge(const vector<Matrix>& xBlocks) const {
		Matrix x(N, 1);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < r; ++j) {
				x(i * r + j, 0) = xBlocks[i](j, 0);
			}
		}
		return x;
	}
	bool converge(const vector<Matrix>& xBlocks, const vector<Matrix>& xpBlocks, double eps) {
		Matrix x, xp;
		x = merge(xBlocks);
		xp = merge(xpBlocks);
		int s = x.GetN();
		double norm = 0;
		for (int i = 0; i < N; i++)
		{
			norm += (x(i, 0) - xp(i, 0))*(x(i, 0) - xp(i, 0));
		}
		if (sqrt(norm) >= eps)
			return false;
		return true;
	}
};


// Код, генерирующий исходную матрицу
// и вычисляющий норму невязки решения.
class TestInstanceGenerator {
public:
	TestInstanceGenerator(int n)
		: n(n)
	{
	}

	Matrix GetA() const {
		srand(n);
		double diag_dominant = 0;
		Matrix a(n, n);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				a(i, j) = rand() % 100 + 1;
				diag_dominant += abs(a(i, j));
			}
			a(i, i) = diag_dominant + 1;
			diag_dominant = 0;
		}
		return a;
	}

	Matrix GetX() const {
		Matrix x(n, 1);
		for (int i = 0; i < n; ++i) {
			x(i, 0) = 1.0;
		}
		return x;
	}

	Matrix GetB() const {
		srand(n);
		Matrix x(n, 1);
		for (int i = 0; i < n; ++i) {
			x(i, 0) = rand() % 2;
		}
		return x;
	}

	double GetError(const Matrix& result) const {
		double maxDelta = 0;
		const Matrix correctResult = GetX();
		for (int i = 0; i < n; ++i) {
			maxDelta = fmax(maxDelta, std::abs(result(i, 0) - correctResult(i, 0)));
		}
		return maxDelta;
	}

private:
	int n;
};


class ExperementPerformer
{
public:
	ExperementPerformer();
	~ExperementPerformer();

	void print_head(int block_size = 0) {
		fout << (block_size == 0 ? "Размер блока = размеру матрицы" : "Размер блока равен ") << block_size << "\t\n";
		std::cout << (block_size == 0 ? "Размер блока = размеру матрицы" : "Размер блока равен ") << block_size << "\t\n";
		fout << "Размер матрицы" << "\t" << "1 поток" << "\t" << omp_get_max_threads() << "потоков" << "\t\n";
		std::cout << "Размер матрицы" << "\t" << "1 поток" << "\t" << omp_get_max_threads() << "потоков" << "\t\n";
	}

	void perform_experement(int blck_size = 0) {
		//block_size - size of block for lu decomposition
		//k - step for matrix size
		int block_size = 0, k = 0;
		Solver *solver;
		Matrix A, B;
		TestInstanceGenerator *generator;
		for (k = 1000; k <= 5000; k += 1000) {
			generator = new TestInstanceGenerator(k);
			A = generator->GetA();
			B = generator->GetB();
			fout << k << "\t";
			std::cout << k << "\t";
			//in one thread			
			block_size = (blck_size == 0) ? k : blck_size;
			solver = new Solver(k, block_size, A, B);
			if (block_size == k) {
				solver->SerialSeidel();
				continue;
			}
			solver->ParallelSeidel(1);
			//in n threads		
			solver = new Solver(k, block_size, A, B);
			solver->ParallelSeidel(omp_get_max_threads());
			//
			fout << std::endl;
			std::cout << std::endl;
		}
		fout << std::endl;
		std::cout << std::endl;
	}
};

ExperementPerformer::ExperementPerformer()
{
	fout.open("results.xls");
}

ExperementPerformer::~ExperementPerformer()
{
	fout.close();
}
int main() {
	int block_sizes[] = { 0,2,5,10,20,50,100 };
	ExperementPerformer *performer = new ExperementPerformer();
	for (int i = 0; i < 6; i++) {
		performer->print_head(block_sizes[i]);
		performer->perform_experement(block_sizes[i]);
	}
	system("pause");
	return 0;
};
