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


//3


//#include <iostream>
//#include <cstdlib>
//#include "mpi.h"
//
//int main(int argc, char** argv)
//{
//	srand(time(NULL));
//	MPI_Init(&argc, &argv);
//	MPI_Status status;
//
//	int rank, size;
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	double startTime = MPI_Wtime();
//
//	int K = atoi(argv[1]);
//	int L = atoi(argv[2]);
//	int M = atoi(argv[3]);
//
//	int whole = K / size;
//	int remainder = K % size;
//	int local_k = (rank < remainder) ? (whole + 1) : whole;
//	int local_l = (rank < remainder) ? (whole + 1) : whole;
//	double* a = new double[K * L];
//	double* b = new double[L * M];
//	double* c = new double[K * M];
//	double* local_a = new double[local_k * L];
//	double* local_c = new double[local_k * M];
//
//	if (rank == 0) {
//		for (int i = 0; i < K * L; i++)
//		{
//			a[i] = rand() % 2 - 1;
//		}
//		for (int i = 0; i < L * M; i++)
//		{
//			b[i] = rand() % 2 - 1;
//		}
//	}
//
//	MPI_Bcast(b, L * M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//	int sum = 0;
//	int* sendcounts = new int[size];
//	int* displs = new int[size];
//
//	for (int i = 0; i < size; i++)
//	{
//		sendcounts[i] = ((i < remainder) ? (whole + 1) : whole) * L;
//		displs[i] = sum;
//		sum += sendcounts[i];
//	}
//
//	MPI_Scatterv(a, sendcounts, displs, MPI_DOUBLE, local_a, local_k * L, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//	for (int i = 0; i < local_k; ++i) {
//		for (int j = 0; j < M; ++j) {
//			local_c[i * M + j] = 0;
//			for (int k = 0; k < L; ++k) {
//				local_c[i * M + j] += local_a[i * L + k] * b[k * M + j];
//			}
//		}
//	}
//
//	sum = 0;
//	for (int i = 0; i < size; i++)
//	{
//		sendcounts[i] = ((i < remainder) ? (whole + 1) : whole) * M;
//		displs[i] = sum;
//		sum += sendcounts[i];
//	}
//
//	MPI_Gatherv(local_c, local_k * M, MPI_DOUBLE, c, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//	if (rank == 0) {
//		std::cout << "Matrix a:" << std::endl;
//		for (int i = 0; i < K * L; i++)
//		{
//			std::cout << a[i] << " ";
//			if ((i + 1) % L == 0) {
//				std::cout << std::endl;
//			}
//		}
//		std::cout << "Matrix b:" << std::endl;
//		for (int i = 0; i < L * M; i++)
//		{
//			std::cout << b[i] << " ";
//			if ((i + 1) % M == 0) {
//				std::cout << std::endl;
//			}
//		}
//		std::cout << "Matrix a * b:" << std::endl;
//		for (int i = 0; i < K * M; i++)
//		{
//			std::cout << c[i] << " ";
//			if ((i + 1) % M == 0) {
//				std::cout << std::endl;
//			}
//		}
//	}
//	double endTime = MPI_Wtime();
//
//	double time = endTime - startTime;
//
//	if (rank == 0) {
//		std::cout << "Start time: " << startTime << std::endl;
//		std::cout << "Finish time: " << endTime << std::endl;
//		std::cout << "Time: " << time << std::endl;
//	}
//	MPI_Finalize();
//	return 0;
//}
