#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include <iostream>
#include <cmath>
#include "BlockSchemaAlgorithm.h"

double BlockSchemaAlgorithm::execute(int argc, char* argv[])
{
	double* pMatrix, * pVector, * pProcVector, * pResult, * pProcMatrix;
	int rank, proc_num, RowNum, ColumnNum;
	double* pProcRows, * pProcResult;

	int Size = atoi(argv[2]);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	ProcessInitialization(pMatrix, pVector, pProcVector, pResult, pProcRows, pProcResult, pProcMatrix,
		Size, RowNum, ColumnNum, proc_num, rank);

	double Start = MPI_Wtime();

	DataDistribution(pMatrix, pProcMatrix, pVector, pProcVector, Size, RowNum, ColumnNum, proc_num, rank);
	ParallelResultCalculation(pProcRows, pVector, pProcResult, Size, RowNum);
	ResultReplication(pProcResult, pResult, Size, RowNum, proc_num, rank);

	double Finish = MPI_Wtime();
	double Duration = Finish - Start;
	if (rank == 0)
		std::cout << "Duration = " << Duration << std::endl;

	ProcessTermination(pMatrix, pVector, pProcVector, pResult, pProcMatrix, pProcResult, rank);

	MPI_Finalize();
	return Duration;
}

void BlockSchemaAlgorithm::DataDistribution(double* pMatrix, double* pProcMatrix, double* pVector, double* pProcVector,
	int Size, int RowNum, int ColumnNum, int& ProcNum, int& ProcRank)
{
	int* pSendNum,* pSendInd,* pSendVecNum,* pSendVecInd;

	pSendInd = new int[ProcNum];
	pSendNum = new int[ProcNum];
	pSendVecInd = new int[ProcNum];
	pSendVecNum = new int[ProcNum];

	pSendNum[0] = RowNum * ColumnNum;
	pSendInd[0] = 0;
	pSendVecNum[0] = ColumnNum;
	pSendVecInd[0] = 0;
	for (int i = 1; i < ProcNum; i++) {
		pSendNum[i] = RowNum * ColumnNum;
		pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
		pSendVecNum[i] = ColumnNum;
		pSendVecInd[i] = (pSendVecInd[i - 1] + pSendVecNum[i - 1]) % Size;
	}

	MPI_Scatterv(pVector, pSendVecNum, pSendVecInd, MPI_DOUBLE, pProcVector, pSendVecNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcMatrix,
		pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

	delete[] pSendNum;
	delete[] pSendInd;
	delete[] pSendVecNum;
	delete[] pSendVecInd;
}

void BlockSchemaAlgorithm::ParallelResultCalculation(double* pProcMatrix, double* pProcVector, double* pProcResult, int RowNum, int ColumnNum) {
	int i, j;
	for (i = 0; i < RowNum; i++)
		pProcResult[i] = 0;
	for (i = 0; i < RowNum; i++) {
		for (j = 0; j < ColumnNum; j++)
			pProcResult[i] += pProcMatrix[i * ColumnNum + j] * pProcVector[j];
	}
}

void BlockSchemaAlgorithm::ProcessInitialization(double*& original, double*& pMatrix, double*& pVector, double*& pProcVector,
	double*& pResult, double*& pProcMatrix, double*& pProcResult,
	int& Size, int& RowNum, int& ColumnNum, int& ProcNum, int& ProcRank) {
	int RestColumns;
	int i;

	setvbuf(stdout, 0, _IONBF, 0);
	if (ProcRank == 0) {
		int s = determineMaxWidth(ProcNum), q = ProcNum / s;
		RowNum = Size / q, ColumnNum = Size / s;
	}

	MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&RowNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ColumnNum, 1, MPI_INT, 0, MPI_COMM_WORLD);

	pProcVector = new double[ColumnNum];
	pProcMatrix = new double[RowNum * ColumnNum];
	pProcResult = new double[RowNum];


	if (ProcRank == 0) {
		pMatrix = new double[Size * Size];
		pVector = new double[Size];
		pResult = new double[Size];
		RandomDataInitialization(pMatrix, pVector, Size);
	}
}

void BlockSchemaAlgorithm::ProcessTermination(double* pMatrix, double* pVector, double* pProcVector, double* pResult,
	double* pProcMatrix, double* pProcResult, int ProcRank) {
	if (ProcRank == 0) {
		delete[] pMatrix;
		delete[] pVector;
		delete[] pResult;
	}
	delete[] pProcVector;
	delete[] pProcResult;
	delete[] pProcMatrix;
}

void BlockSchemaAlgorithm::RandomDataInitialization(double*& pMatrix, double*& pVector, int& Size) {
	int i, j;
	srand(unsigned(clock()));
	for (i = 0; i < Size; i++) {
		pVector[i] = rand() / double(1000);
		for (j = 0; j < Size; j++) {
			double randVal = rand() / double(1000);
			pMatrix[i * Size + j] = randVal;
		}
	}
}

void BlockSchemaAlgorithm::ResultReplication(double* pProcResult, double* pResult, int Size, int RowNum, int ColumnNum, int ProcRank) {

	int posInRes = (ProcRank / (Size / ColumnNum)) * RowNum;
	double* allProcResult = new double[Size];
	for (int i = 0; i < Size; i++)
		allProcResult[i] = 0;
	for (int i = posInRes; i < posInRes + RowNum; i++) {
		allProcResult[i] = pProcResult[i - posInRes];
	}

	MPI_Reduce(allProcResult, pResult, Size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}
int BlockSchemaAlgorithm::determineMaxWidth(int Size)
{
	for (int i = sqrt(Size); i >= 1; i--)
		if (Size % i == 0)
			return i;
	return 0;
}
