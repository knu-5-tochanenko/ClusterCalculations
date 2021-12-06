#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpi.h>

class MatrixOperations {
public:
	static void DummyDataInitialization(double *pAMatrix, double *pBMatrix, int Size) {
		int i, j; // Loop variables
		for (i = 0; i < Size; i++)
			for (j = 0; j < Size; j++) {
				pAMatrix[i * Size + j] = 1;
				pBMatrix[i * Size + j] = 1;
			}
	}

	static void RandomDataInitialization(double *pAMatrix, double *pBMatrix,
										 int Size) {
		int i, j; // Loop variables
		srand(unsigned(clock()));
		for (i = 0; i < Size; i++)
			for (j = 0; j < Size; j++) {
				pAMatrix[i * Size + j] = rand() / double(1000);
				pBMatrix[i * Size + j] = rand() / double(1000);
			}
	}

	static void PrintMatrix(double *pMatrix, int RowCount, int ColCount) {
		int i, j; // Loop variables
		for (i = 0; i < RowCount; i++) {
			for (j = 0; j < ColCount; j++)
				printf("%7.4f ", pMatrix[i * ColCount + j]);
			printf("\n");
		}
	}
};

int ProcNum = 0;
int ProcRank = 0;
int GridSize;
int GridCoords[2];
MPI_Comm GridComm;
MPI_Comm ColComm;
MPI_Comm RowComm;

// Function for matrix multiplication
void SerialResultCalculation(const double *pAMatrix, const double *pBMatrix,
                             double *pCMatrix, int Size) {
    int i, j, k; // Loop variables
    for (i = 0; i < Size; i++) {
        for (j = 0; j < Size; j++)
            for (k = 0; k < Size; k++)
                pCMatrix[i * Size + j] += pAMatrix[i * Size + k] * pBMatrix[k * Size + j];
    }
}

// Function for block multiplication
void BlockMultiplication(double *pAblock, double *pBblock,
                         double *pCblock, int Size) {
    SerialResultCalculation(pAblock, pBblock, pCblock, Size);
}

// Creation of two-dimensional grid communicator
// and communicators for each row and each column of the grid
void CreateGridCommunicators() {
    int DimSize[2]; // Number of processes in each dimension of the grid
    int Periodic[2]; // =1, if the grid dimension should be periodic
    int Subdims[2]; // =1, if the grid dimension should be fixed

    DimSize[0] = GridSize;
    DimSize[1] = GridSize;
    Periodic[0] = 0;
    Periodic[1] = 0;
    // Creation of the Cartesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, 2, DimSize, Periodic, 1, &GridComm);
    // Determination of the cartesian coordinates for every process
    MPI_Cart_coords(GridComm, ProcRank, 2, GridCoords);

    // Creating communicators for rows
    Subdims[0] = 0; // Dimensionality fixing
    Subdims[1] = 1; // The presence of the given dimension in the subgrid
    MPI_Cart_sub(GridComm, Subdims, &RowComm);

    // Creating communicators for columns
    Subdims[0] = 1;
    Subdims[1] = 0;
    MPI_Cart_sub(GridComm, Subdims, &ColComm);
}

// Function for memory allocation and data initialization
void ProcessInitialization(double *&pAMatrix, double *&pBMatrix,
                           double *&pCMatrix, double *&pAblock, double *&pBblock, double *&pCblock,
                           double *&pTemporaryAblock, int &Size, int &BlockSize) {
    if (ProcRank == 0) {
        do {
            printf("\nEnter size of the initial objects: ");
            scanf("%d", &Size);
            if (Size % GridSize != 0) {
                printf("Size of matricies must be divisible by the grid size!\n");
            }
        } while (Size % GridSize != 0);
    }
    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    BlockSize = Size / GridSize;
    pAblock = new double[BlockSize * BlockSize];
    pBblock = new double[BlockSize * BlockSize];
    pCblock = new double[BlockSize * BlockSize];
    pTemporaryAblock = new double[BlockSize * BlockSize];
    for (int i = 0; i < BlockSize * BlockSize; i++) {
        pCblock[i] = 0;
    }
    if (ProcRank == 0) {
        pAMatrix = new double[Size * Size];
        pBMatrix = new double[Size * Size];
        pCMatrix = new double[Size * Size];
		MatrixOperations::RandomDataInitialization(pAMatrix, pBMatrix, Size);
    }
}

// Function for checkerboard matrix decomposition
void CheckerboardMatrixScatter(double *pMatrix, double *pMatrixBlock,
                               int Size, int BlockSize) {
    auto *MatrixRow = new double[BlockSize * Size];
    if (GridCoords[1] == 0) {
        MPI_Scatter(pMatrix, BlockSize * Size, MPI_DOUBLE, MatrixRow,
                    BlockSize * Size, MPI_DOUBLE, 0, ColComm);
    }
    for (int i = 0; i < BlockSize; i++) {
        MPI_Scatter(&MatrixRow[i * Size], BlockSize, MPI_DOUBLE,
                    &(pMatrixBlock[i * BlockSize]), BlockSize, MPI_DOUBLE, 0, RowComm);
    }
    delete[] MatrixRow;
}

// Data distribution among the processes
void DataDistribution(double *pAMatrix, double *pBMatrix, double *
pMatrixAblock, double *pBblock, int Size, int BlockSize) {
    // Scatter the matrix among the processes of the first grid column
    CheckerboardMatrixScatter(pAMatrix, pMatrixAblock, Size, BlockSize);
    CheckerboardMatrixScatter(pBMatrix, pBblock, Size, BlockSize);
}

// Function for gathering the result matrix
void ResultCollection(double *pCMatrix, double *pCblock, int Size, int BlockSize) {
    auto *pResultRow = new double[Size * BlockSize];
    for (int i = 0; i < BlockSize; i++) {
        MPI_Gather(&pCblock[i * BlockSize], BlockSize, MPI_DOUBLE,
                   &pResultRow[i * Size], BlockSize, MPI_DOUBLE, 0, RowComm);
    }
    if (GridCoords[1] == 0) {
        MPI_Gather(pResultRow, BlockSize * Size, MPI_DOUBLE, pCMatrix,
                   BlockSize * Size, MPI_DOUBLE, 0, ColComm);
    }
    delete[] pResultRow;
}

// Broadcasting matrix A blocks to process grid rows
void ABlockCommunication(int iter, double *pAblock, const double *pMatrixAblock,
                         int BlockSize) {
    // Defining the leading process of the process grid row
    int Pivot = (GridCoords[0] + iter) % GridSize;

    // Copying the transmitted block in a separate memory buffer
    if (GridCoords[1] == Pivot) {
        for (int i = 0; i < BlockSize * BlockSize; i++)
            pAblock[i] = pMatrixAblock[i];
    }

    // Block broadcasting
    MPI_Bcast(pAblock, BlockSize * BlockSize, MPI_DOUBLE, Pivot, RowComm);
}

// Cyclic shift of matrix B blocks in the process grid columns
void BblockCommunication(double *pBblock, int BlockSize) {
    MPI_Status Status;
    int NextProc = GridCoords[0] + 1;
    if (GridCoords[0] == GridSize - 1) NextProc = 0;
    int PrevProc = GridCoords[0] - 1;
    if (GridCoords[0] == 0) PrevProc = GridSize - 1;
    MPI_Sendrecv_replace(pBblock, BlockSize * BlockSize, MPI_DOUBLE,
                         PrevProc, 0, NextProc, 0, ColComm, &Status);
}

void ParallelResultCalculation(double *pAblock, double *pMatrixAblock,
                               double *pBblock, double *pCblock, int BlockSize) {
    for (int iter = 0; iter < GridSize; iter++) {
        // Sending blocks of matrix A to the process grid rows
        ABlockCommunication(iter, pAblock, pMatrixAblock, BlockSize);
        // Block multiplication
        BlockMultiplication(pAblock, pBblock, pCblock, BlockSize);
        // Cyclic shift of blocks of matrix B in process grid columns
        BblockCommunication(pBblock, BlockSize);
    }
}

// Test printing of the matrix block
void TestBlocks(double *pBlock, int BlockSize, char str[]) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (ProcRank == 0) {
        printf("%s \n", str);
    }
    for (int i = 0; i < ProcNum; i++) {
        if (ProcRank == i) {
            printf("ProcRank = %d \n", ProcRank);
			MatrixOperations::PrintMatrix(pBlock, BlockSize, BlockSize);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void TestResult(double *pAMatrix, double *pBMatrix, double *pCMatrix,
                int Size) {
    double StartSerial, FinishSerial, DurationSerial;

    double *pSerialResult; // Result matrix of serial multiplication
    double Accuracy = 0.5; // Comparison accuracy
    int equal = 0; // =1, if the matrices are not equal
    int i; // Loop variable
    if (ProcRank == 0) {
        pSerialResult = new double[Size * Size];
        for (i = 0; i < Size * Size; i++) {
            pSerialResult[i] = 0;
        }
        StartSerial = clock();
        BlockMultiplication(pAMatrix, pBMatrix, pSerialResult, Size);
        FinishSerial = clock();
        DurationSerial = (FinishSerial - StartSerial) / double(CLOCKS_PER_SEC);
        if (ProcRank == 0) {
            printf("Time of serial execution = % f\n", DurationSerial);
        }
        //PrintMatrix(pAMatrix, Size, Size);
        printf_s("\n");
        for (i = 0; i < Size * Size; i++) {
            if (fabs(pSerialResult[i] - pCMatrix[i]) >= Accuracy)
                equal = 1;
        }
        if (equal == 1)
            printf("The results of serial and parallel algorithms are NOT "
                   "identical. Check your code.\n");
        else
            printf("The results of serial and parallel algorithms are "
                   "identical.\n");
    }
}

// Function for computational process termination
void ProcessTermination(const double *pAMatrix, const double *pBMatrix,
                        const double *pCMatrix, const double *pAblock, const double *pBblock, const double *pCblock,
                        const double *pMatrixAblock) {
    if (ProcRank == 0) {
        delete[] pAMatrix;
        delete[] pBMatrix;
        delete[] pCMatrix;
    }
    delete[] pAblock;
    delete[] pBblock;
    delete[] pCblock;
    delete[] pMatrixAblock;
}

int main(int argc, char *argv[]) {
    double *pAMatrix; // The first argument of matrix multiplication
    double *pBMatrix; // The second argument of matrix multiplication
    double *pCMatrix; // The result matrix
    int Size; // Size of matricies
    int BlockSize; // Sizes of matrix blocks on current process
    double *pAblock; // Initial block of matrix A on current process
    double *pBblock; // Initial block of matrix B on current process
    double *pCblock; // Block of result matrix C on current process
    double *pMatrixAblock;
    double Start, Finish, Duration;
    setvbuf(stdout, 0, _IONBF, 0);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    GridSize = sqrt((double) ProcNum);
    if (ProcNum != GridSize * GridSize) {
        if (ProcRank == 0) {
            printf("Number of processes must be a perfect square \n");
        }
    } else {
        if (ProcRank == 0)
            printf("Parallel matrix multiplication program\n");
        // Creating the cartesian grid, row and column communcators
        CreateGridCommunicators();

        // Memory allocation and initialization of matrix elements
        ProcessInitialization(pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock,
                              pCblock, pMatrixAblock, Size, BlockSize);
        Start = MPI_Wtime();
        DataDistribution(pAMatrix, pBMatrix, pMatrixAblock, pBblock, Size,
                         BlockSize);
        // Execution of Fox method
        ParallelResultCalculation(pAblock, pMatrixAblock, pBblock,
                                  pCblock, BlockSize);
        ResultCollection(pCMatrix, pCblock, Size, BlockSize);
        Finish = MPI_Wtime();
        Duration = Finish - Start;
        TestResult(pAMatrix, pBMatrix, pCMatrix, Size);

        //PrintMatrix(pAMatrix, Size, Size);
        if (ProcRank == 0) {
            printf("Time of execution = % f\n", Duration);
        }
        // Process Termination
        ProcessTermination(pAMatrix, pBMatrix, pCMatrix, pAblock, pBblock,
                           pCblock, pMatrixAblock);
    }
    MPI_Finalize();
}