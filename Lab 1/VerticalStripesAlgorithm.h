#pragma once
class VerticalStripesAlgorithm
{
public:
	double execute(int argc, char* argv[]);

private:
	void DataDistribution(double* pMatrix, double* pProcRows, double* pVector, int Size,
		int RowNum, int ProcNum, int ProcRank);

	void ParallelResultCalculation(double* pProcRows, double* pVector, double* pProcResult, int Size,
		int RowNum);

	void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult,
		double*& pProcRows, double*& pProcResult, int& Size, int& RowNum, int& ProcNum, int& ProcRank);

	void ProcessTermination(double* pMatrix, double* pVector, double* pResult,
		double* pProcRows, double* pProcResult, int ProcRank);

	void RandomDataInitialization(double*& pMatrix, double*& pVector, int& Size);

	void ResultReplication(double* pProcResult, double* pResult, int Size, int RowNum, int ProcNum, int ProcRank);
};

