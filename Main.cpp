#include "Algorithms.h"
#include "Generator.h"
/*
Structure of config file:
0		#mode:
40		#size
0		#algorithm: 0-Prim, 1-Kruskal, 2-Dijkstra, 3-Ford
50 		#density: 0-100 in percent
*/
// FileData structure to hold the configuration data
struct fileData {
	int mode;
	int size;
	int algorithm;
	int density;
	int intervalAmount;
	int instanceAmount;
	int representation;
};
// resultsInfo structure to hold the results of the sorting algorithms
struct resultsInfo {
	double avgTime;
	std::vector<double> instanceTime;
};
int main() {
	Generator generator(1000);
	auto matrix = generator.generateWeightedIncidenceMatrix(10);
	return 0;
}