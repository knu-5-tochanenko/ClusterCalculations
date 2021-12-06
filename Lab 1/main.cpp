#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <cstring>
#include "BlockSchemaAlgorithm.h"
#include "HorisontalStripesAlgorithm.h"
#include "SerialAlgorithm.h"
#include "VerticalStripesAlgorithm.h"

int main(int argc, char* argv[]) {
	if (argv[1] == NULL || strlen(argv[1]) == 0) {
		int elements = 0;
		std::cout << "Number of elements:\n";
		std::cin >> elements;
		auto algo = new SerialAlgorithm();
		algo -> execute(argc, argv, elements);
	}
	else switch (atoi(argv[1])) {
		case 0: {
			printf("Block algorithm");
			auto algo = new BlockSchemaAlgorithm();
			algo -> execute(argc, argv);
			break;
		}
		case 1: {
			auto algo = new HorisontalStripesAlgorithm();
			algo -> execute(argc, argv);
			break;
		}
		case 2: {
			auto algo = new VerticalStripesAlgorithm();
			algo -> execute(argc, argv);
			break;
		}
		default: {
			printf("No algorithm found.");
		}
	}

    return 0;
}
