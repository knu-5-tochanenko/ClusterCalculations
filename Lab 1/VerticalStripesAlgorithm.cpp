#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include <iostream>
#include "VerticalStripesAlgorithm.h"

double VerticalStripesAlgorithm::execute(int argc, char* argv[])
{
    double* pMatrix, * pVector, * pResult;
    int rank, proc_num, RowNum;
    double* pProcRows, * pProcResult;

    int Size = atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ProcessInitialization(pMatrix, pVector, pResult, pProcRows, pProcResult,  Size, RowNum, proc_num, rank);

    double Start = MPI_Wtime();
    DataDistribution(pMatrix, pProcRows, pVector, Size, RowNum, proc_num, rank);
    ParallelResultCalculation(pProcRows, pVector, pProcResult, Size, RowNum);
    ResultReplication(pProcResult, pResult, Size, RowNum, proc_num, rank);

    double Finish = MPI_Wtime();
    double Duration = Finish - Start;
    if (rank == 0)
        std::cout << "Duration = " << Duration << std::endl;

    ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcResult, rank);

    MPI_Finalize();
    return Duration;
}

void VerticalStripesAlgorithm::DataDistribution(double* pMatrix, double* pProcRows, double* pVector, int Size,
    int RowNum, int ProcNum, int ProcRank) {
    int* pSendNum,* pSendInd;
    int RestRows = Size;
    MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    pSendInd = new int[ProcNum];
    pSendNum = new int[ProcNum];
    RowNum = (Size / ProcNum);
    pSendNum[0] = RowNum * Size;
    pSendInd[0] = 0;
    for (int i = 1; i < ProcNum; i++) {
        RestRows -= RowNum;
        RowNum = RestRows / (ProcNum - i);
        pSendNum[i] = RowNum * Size;
        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
    }
    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows,
        pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    delete[] pSendNum;
    delete[] pSendInd;
}

void VerticalStripesAlgorithm::ParallelResultCalculation(double* pProcRows, double* pVector, double* pProcResult, int Size,
    int RowNum) {
    int i, j;
    for (i = 0; i < RowNum; i++) {
        pProcResult[i] = 0;
        for (j = 0; j < Size; j++)
            pProcResult[i] += pProcRows[i * Size + j] * pVector[j];
    }
}

void VerticalStripesAlgorithm::ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult,
    double*& pProcRows, double*& pProcResult, int& Size, int& RowNum, int& ProcNum, int& ProcRank) {
    int RestRows, i;
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    RestRows = Size;
    for (i = 0; i < ProcRank; i++)
        RestRows = RestRows - RestRows / (ProcNum - i);
    RowNum = RestRows / (ProcNum - ProcRank);
    pVector = new double[Size];
    pResult = new double[Size];
    pProcRows = new double[RowNum * Size];
    pProcResult = new double[RowNum];
    if (ProcRank == 0) {
        pMatrix = new double[Size * Size];
        RandomDataInitialization(pMatrix, pVector, Size);
    }
}

void VerticalStripesAlgorithm::ProcessTermination(double* pMatrix, double* pVector, double* pResult,
    double* pProcRows, double* pProcResult, int ProcRank) {
    if (ProcRank == 0)
        delete[] pMatrix;
    delete[] pVector;
    delete[] pResult;
    delete[] pProcRows;
    delete[] pProcResult;
}

void VerticalStripesAlgorithm::RandomDataInitialization(double*& pMatrix, double*& pVector, int& Size) {
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

void VerticalStripesAlgorithm::ResultReplication(double* pProcResult, double* pResult, int Size, int RowNum, int ProcNum, int ProcRank) {
    int i;
    int* pReceiveNum,* pReceiveInd;
    int RestRows = Size;
    pReceiveNum = new int[ProcNum];
    pReceiveInd = new int[ProcNum];
    pReceiveInd[0] = 0;
    pReceiveNum[0] = Size / ProcNum;
    for (i = 1; i < ProcNum; i++) {
        RestRows -= pReceiveNum[i - 1];
        pReceiveNum[i] = RestRows / (ProcNum - i);
        pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
    }
    MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_DOUBLE, pResult,
        pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
    delete[] pReceiveNum;
    delete[] pReceiveInd;
}
