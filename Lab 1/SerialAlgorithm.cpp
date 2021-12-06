#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <ctime>  
#include "SerialAlgorithm.h"

double SerialAlgorithm::execute(int argc, char* argv[], int size)
{

	double* pMatrix, * pVector, * pResult;

	int Size = size;

	pVector = new double[Size];
	pResult = new double[Size];
	pMatrix = new double[Size * Size];
	RandomDataInitialization(pMatrix, pVector, Size);
	auto Start = std::chrono::system_clock::now();
	SequentialResultCalculation(pMatrix, pVector, pResult, Size);
	auto Finish = std::chrono::system_clock::now();
	double Duration = (Finish - Start).count();

	std::cout << "Duration = " << Duration << std::endl;

	return Duration;
}

void SerialAlgorithm::RandomDataInitialization(double*& pMatrix, double*& pVector, int& Size)
{
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

void SerialAlgorithm::SequentialResultCalculation(double* pMatrix, double* pVector, double* pResult, int Size)
{
	int i, j; 
	for (i = 0; i < Size; i++)
		pResult[i] = 0;
	for (j = 0; j < Size; j++) {
		for (i = 0; i < Size; i++)
			pResult[i] += pMatrix[j * Size + i] * pVector[j];
	}
}
