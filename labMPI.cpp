//Lab MPI

//#include "mpi.h"
//#include <iostream>
//#include <random>
//#include <vector>
//
//int main(int argc, char** argv)
//{
//	srand(time(NULL));
//	int rank, size;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	int n = 5;
//	//std::cin >> n;
//
//	int n1 = n / size;
//	int res = 0;
//	int res1 = 0;
//
//	std::vector<int> a(n);
//	std::vector<int> b(n);
//	std::vector<int> a1(n);
//	std::vector<int> b1(n);
//
//	if (rank == 0) {
//		for (int i = 0; i < n; i++)
//		{
//			a[i] = std::rand() % 2;
//			b[i] = std::rand() % 2;
//		}
//	}
//
//	//  Распределение данных с нулевого процесса
//	MPI_Scatter(&a[0], n1, MPI_INT, &a1[0], n1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Scatter(&b[0], n1, MPI_INT, &b1[0], n1, MPI_INT, 0, MPI_COMM_WORLD);
//
//	for (int i = 0; i < n1; ++i) {
//		res1 += a1[i] * b1[i];
//	}
//
//	//Сбор данных
//	MPI_Gather(&res1, n1, MPI_INT, &res, n1, MPI_INT, 0, MPI_COMM_WORLD);
//	if (rank == 0) {
//		std::cout << "Vector a: ";
//		for (int i = 0; i < n; ++i) {
//			std::cout << a[i] << ' ';
//		}
//
//		std::cout << std::endl;
//		std::cout << "Vector b: ";
//		for (int i = 0; i < n; ++i) {
//			std::cout << b[i] << ' ';
//		}
//		std::cout << std::endl;
//		std::cout << "Result: " << res << std::endl;
//	}
//	MPI_Finalize();
//	return 0;
//}

#include <iostream>
#include "mpi.h"
#include <cstdlib>

int main(int argc, char** argv)
{
	srand(time(NULL));
	MPI_Init(&argc, &argv);
	MPI_Status status;

	// Получаем размер группы и индекс текущего процесса
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// размеры матрицы (L- линии, M-столбцы)
	int L=3; 
	int M=2;
	
	// пересылаем размеры матрицы 
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	int n1 = L / size;
	int r = L % size;

	double* a = new double[L * M];
	double* b = new double[L * M];
	double* c = new double[L * M];

	//количество строк, которые каждый процесс должен обработать
	int local_n = (rank < r) ? (n1 + 1) : n1;

	double* a1 = new double[local_n * M];
	double* b1 = new double[local_n * M];
	double* c1 = new double[local_n * M];

	//количестве элементов, отправляемых каждому процессу, и смещениe в исходном массиве
	int* sendcounts = new int[size];
	int* displs = new int[size];

	// заполнение матриц
	if (rank == 0) {
		for (int i = 0; i < L * M; i++)
		{
			a[i] = rand() % 10;
			b[i] = rand() % 10;
		}
	}

	int sum = 0;
	for (int i = 0; i < size; i++)
	{
		sendcounts[i] = ((i < r) ? (n1 + 1) : n1) * M;
		displs[i] = sum;
		sum += sendcounts[i];
	}
	
	MPI_Scatterv(&a[0], sendcounts, displs, MPI_DOUBLE, a1, local_n * M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(&b[0], sendcounts, displs, MPI_DOUBLE, b1, local_n * M, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = 0; i < local_n * M; i++)
	{
		c1[i] = a1[i] + b1[i];
	}

	MPI_Gatherv(c1, local_n * M, MPI_DOUBLE, c, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		std::cout << "Matrix a:" << std::endl;
		for (int i = 0; i < L * M; i++)
		{
			std::cout << a[i] << " ";
			if ((i + 1) % M == 0) {
				std::cout << std::endl;
			}
		}
		std::cout << "Matrix b:" << std::endl;
		for (int i = 0; i < L * M; i++)
		{
			std::cout << b[i] << " ";
			if ((i + 1) % M == 0) {
				std::cout << std::endl;
			}
		}
		std::cout << "Matrix a + b:" << std::endl;
		for (int i = 0; i < L * M; i++)
		{
			std::cout << c[i] << " ";
			if ((i + 1) % M == 0) {
				std::cout << std::endl;
			}
		}
	}

	MPI_Finalize();

	return 0;
}
