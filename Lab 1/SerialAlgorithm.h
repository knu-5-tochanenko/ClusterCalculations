#pragma once
class SerialAlgorithm
{
public:
	double execute(int argc, char* argv[], int size);
private:
	void RandomDataInitialization(double*& pMatrix, double*& pVector, int& Size);

	void SequentialResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size);
};

