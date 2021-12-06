#pragma once
class BlockSchemaAlgorithm
{
public:
	double execute(int argc, char* argv[]);

private:
	void DataDistribution(double* pMatrix, double* pProcMatrix, double* pVector, double* pProcVector,
		int Size, int RowNum, int ColumnNum, int& ProcNum, int& ProcRank);

	void ParallelResultCalculation(double* pProcMatrix, double* pProcVector, double* pProcResult, int RowNum, int ColumnNum);

	void ProcessInitialization(double*& original, double*& pMatrix, double*& pVector, double*& pProcVector,
		double*& pResult, double*& pProcMatrix, double*& pProcResult,
		int& Size, int& RowNum, int& ColumnNum, int& ProcNum, int& ProcRank);

	void ProcessTermination(double* pMatrix, double* pVector, double* pProcVector, double* pResult,
		double* pProcMatrix, double* pProcResult, int ProcRank);

	void RandomDataInitialization(double*& pMatrix, double*& pVector, int& Size);

	void ResultReplication(double* pProcResult, double* pResult, int Size, int RowNum, int ProcNum, int ProcRank);

	int determineMaxWidth(int Size);
};

