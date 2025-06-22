#include "Algorithms.h"
#include "Generator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <limits>
/*
Structure of config file:
0		#mode: 0-test, 1-simulation
40		#initialSize
7		#sizeAmount: amount of different sizes to test
50 		#density: 0-100 in percent
100		#maxRand: maximum weight for edges
40		#intervalAmount: difference in vertix size between instances
50		#instanceAmount: number of instances to generate
0		#task: 0-find MST, 1-find shortest paths
0		#displayGraph: 0-no, 1-yes
0		#algorithm in test mode: 0-Prim, 1-Kruskal, 2-Dijkstra, 3-Ford, Dijkstra and Ford only for directed graphs, Prim and Kruskal only for undirected graphs
0		#test graph source, 0-load from file 1-generate random
*/

// FileData structure to hold the configuration data
struct fileData {
	int mode;
	int initialSize;
	int sizeAmount;
	double density;
	int maxRand;
	int intervalAmount;
	int instanceAmount;
	int representation;
	int task;
	int displayGraph;
	int algorithm;
	int testGraphSource;
};
struct mstSimData {
	int size;
	int instanceAmount;
	double density;
	int maxRand;
};
struct mstResultsData {
	std::vector<double> primTimeAL;
	std::vector<double> kruskalTimeAL;
	std::vector<double> primTimeIM;
	std::vector<double> kruskalTimeIM;
};
struct spResultsData {
	std::vector<double> dijkstraTimeAL;
	std::vector<double> bellmanFordTimeAL;
	std::vector<double> dijkstraTimeIM;
	std::vector<double> bellmanFordTimeIM;
};
struct algorithmAvgTimes {
	double primAL;
	double kruskalAL;
	double primIM;
	double kruskalIM;
	double dijkstraAL;
	double bellmanFordAL;
	double dijkstraIM;
	double bellmanFordIM;
};
void performTest(const fileData fileData);
std::vector<mstResultsData> performMSTSimulation(mstSimData mstSimData);
std::vector<spResultsData> performSPSimulation(mstSimData mstSimData);
void saveToCSV(const std::vector<algorithmAvgTimes>& avgTimes) {
	std::cout << "Enter the path and filename to save results (e.g., C:/Users/userName/Desktop/results.csv): ";
	std::string fileName;
	std::cin >> fileName;
	std::ofstream outFile(fileName);

	if (!outFile.is_open()) {
		std::cerr << "Error opening file for writing: " << fileName << std::endl;
		return;
	}
	outFile << "PrimAL;KruskalAL;PrimIM;KruskalIM;DijkstraAL;BellmanFordAL;DijkstraIM;BellmanFordIM\n";
	for (const auto& times : avgTimes) {
		outFile << times.primAL << ";"
				<< times.kruskalAL << ";"
				<< times.primIM << ";"
				<< times.kruskalIM << ";"
				<< times.dijkstraAL << ";"
				<< times.bellmanFordAL << ";"
				<< times.dijkstraIM << ";"
				<< times.bellmanFordIM << "\n";
	}
	outFile.close();
}
void performTest(fileData fileData) {
	// Load graph from file and generate AL and IM or generate random graph based on parameters
	// Display graphs if required
	// Run algorithms and display results
	Generator generator(fileData.maxRand);
	if (fileData.task == 0) {
		if (fileData.testGraphSource == 0) {
			std::string filePath;
			std::cout << "Enter the path to the graph file: ";
			std::cin >> filePath;
			auto matrix = generator.generateWeightedIncidenceMatrixFromFile(filePath);
			auto adjacencyList = generator.generateWeightedAdjacencyListFromFile(filePath);
			if (fileData.displayGraph == 1) {
				generator.printIncidenceMatrix(matrix, fileData.initialSize);
				generator.printAdjacencyList(adjacencyList, fileData.initialSize);
			}
			Algorithms algorithms;
			if (fileData.algorithm == 0) {
				auto mstResult = algorithms.primAlgorithmIM(matrix, fileData.initialSize);
				algorithms.printMST(mstResult);
			} else if (fileData.algorithm == 1) {
				auto mstResult = algorithms.kruskalAlgorithmIM(matrix, fileData.initialSize);
				algorithms.printMST(mstResult);
			} else {
				std::cerr << "Invalid algorithm for MST task." << std::endl;
			}
		}
		else if (fileData.testGraphSource == 1) {
			auto adjacencyList = generator.generateWeightedAdjacencyList(fileData.initialSize);
			auto incidenceMatrix = generator.generateWeightedIncidenceMatrix(fileData.initialSize);
			// Reduce density if specified
			adjacencyList = generator.reduceAdjacencyListDensity(adjacencyList, fileData.density, fileData.initialSize);
			incidenceMatrix = generator.reduceIncidenceMatrixDensity(incidenceMatrix, fileData.density, fileData.initialSize);
			if (fileData.displayGraph == 1) {
				generator.printAdjacencyList(adjacencyList, fileData.initialSize);
				generator.printIncidenceMatrix(incidenceMatrix, fileData.initialSize);
			}
			Algorithms algorithms;
			if (fileData.algorithm == 0) {
				auto mstResult = algorithms.primAlgorithmAL(adjacencyList, fileData.initialSize);
				algorithms.printMST(mstResult);
			} else if (fileData.algorithm == 1) {
				auto mstResult = algorithms.kruskalAlgorithmAL(adjacencyList, fileData.initialSize);
				algorithms.printMST(mstResult);
			} else {
				std::cerr << "Invalid algorithm for MST task." << std::endl;
			}
			generator.deleteAdjacencyList(adjacencyList);
		}
	}
	else if (fileData.task == 1) {
		if (fileData.testGraphSource == 0) {
			std::string filePath;
			std::cout << "Enter the path to the graph file: ";
			std::cin >> filePath;
			auto matrix = generator.generateWeightedDirectedIncidenceMatrixFromFile(filePath);
			auto adjacencyList = generator.generateWeightedDirectedAdjacencyListFromFile(filePath);
			if (fileData.displayGraph == 1) {
				generator.printIncidenceMatrix(matrix, fileData.initialSize);
				generator.printAdjacencyList(adjacencyList, fileData.initialSize);
			}
			Algorithms algorithms;
			if (fileData.algorithm == 2) {
				auto pathResult = algorithms.dijkstraAlgorithmIM(matrix, fileData.initialSize);
				algorithms.printPaths(pathResult);
			} else if (fileData.algorithm == 3) {
				auto pathResult = algorithms.bellmanFordAlgorithmIM(matrix, fileData.initialSize);
				algorithms.printPaths(pathResult);
			} else {
				std::cerr << "Invalid algorithm for SP task." << std::endl;
			}
		}
		else if (fileData.testGraphSource == 1) {
			auto adjacencyList = generator.generateWeightedDirectedAdjacencyList(fileData.initialSize);
			auto incidenceMatrix = generator.generateWeightedDirectedIncidenceMatrix(fileData.initialSize);
			// Reduce density if specified
			adjacencyList = generator.reduceDirectedAdjacencyListDensity(adjacencyList, fileData.density, fileData.initialSize);
			incidenceMatrix = generator.reduceIncidenceMatrixDensity(incidenceMatrix, fileData.density, fileData.initialSize);
			if (fileData.displayGraph == 1) {
				generator.printAdjacencyList(adjacencyList, fileData.initialSize);
				generator.printIncidenceMatrix(incidenceMatrix, fileData.initialSize);
			}
			Algorithms algorithms;
			if (fileData.algorithm == 2) {
				auto pathResult = algorithms.dijkstraAlgorithmAL(adjacencyList, fileData.initialSize);
				algorithms.printPaths(pathResult);
			} else if (fileData.algorithm == 3) {
				auto pathResult = algorithms.bellmanFordAlgorithmAL(adjacencyList, fileData.initialSize);
				algorithms.printPaths(pathResult);
			} else {
				std::cerr << "Invalid algorithm for SP task." << std::endl;
			}
			generator.deleteAdjacencyList(adjacencyList);
		}
	} else {
		std::cerr << "Invalid task in config file." << std::endl;
		exit(EXIT_FAILURE);
	}
}
fileData loadFileData() {
	fileData fileData;
	std::string filePath;
	std::fstream fileStream;
	std::cout << "Enter the path to the config file: ";
	std::cin >> filePath;
	fileStream.open(filePath);
	if (!fileStream.is_open()) {
		std::cerr << "Error opening file: " << filePath << std::endl;
		exit(EXIT_FAILURE);
		return fileData;
	}
	std::string line;
	int lineIndex = 0;
	while (std::getline(fileStream, line)) {
		if (line.empty() || line[0] == '#') {
			continue;
		}
		std::istringstream iss(line);
		int value;
		if (iss >> value) {
			switch (lineIndex) {
				case 0: fileData.mode = value; break;
				case 1: fileData.initialSize = value; break;
				case 2: fileData.sizeAmount = value; break;
				case 3: fileData.density = value/100.0; break;
				case 4: fileData.maxRand = value; break;
				case 5: fileData.intervalAmount = value; break;
				case 6: fileData.instanceAmount = value; break;
				case 7: fileData.task = value; break;
				case 8: fileData.displayGraph = value; break;
				case 9: fileData.algorithm = value; break;
				case 10: fileData.testGraphSource = value; break;
				default:
					std::cerr << "Unexpected line in config file: " << line << std::endl;
					exit(EXIT_FAILURE);
			}
			lineIndex++;
		}
	}
	fileStream.close();
	return fileData;
}
double calculateAverageTime(const std::vector<double>& results) {
	double sum = 0.0;
	for (double time : results) {
		sum += time;
	}
	return std::round((sum / results.size()) * 100.0) / 100.0; //Round to 2 decimal places
}
std::vector<mstResultsData> performMSTSimulation(mstSimData simData) {
	// Create one result structure to store all times for this size
	std::vector<mstResultsData> results(1);  // Only one element needed
	mstResultsData& times = results[0];

	Generator generator(simData.maxRand);
	Algorithms algorithms;

	// Reserve space for better performance
	times.primTimeAL.reserve(simData.instanceAmount);
	times.kruskalTimeAL.reserve(simData.instanceAmount);
	times.primTimeIM.reserve(simData.instanceAmount);
	times.kruskalTimeIM.reserve(simData.instanceAmount);

	// Process each instance individually to save memory
	for (int i = 0; i < simData.instanceAmount; ++i) {
		// Generate and test adjacency list
		auto adjacencyList = generator.generateWeightedAdjacencyList(simData.size);
		adjacencyList = generator.reduceAdjacencyListDensity(adjacencyList, simData.density, simData.size);

		// Time Prim's algorithm on adjacency list
		auto start = std::chrono::high_resolution_clock::now();
		auto primResultAL = algorithms.primAlgorithmAL(adjacencyList, simData.size);
		auto end = std::chrono::high_resolution_clock::now();
		times.primTimeAL.push_back(std::chrono::duration<double, std::milli>(end - start).count());

		// Time Kruskal's algorithm on adjacency list
		start = std::chrono::high_resolution_clock::now();
		auto kruskalResultAL = algorithms.kruskalAlgorithmAL(adjacencyList, simData.size);
		end = std::chrono::high_resolution_clock::now();
		times.kruskalTimeAL.push_back(std::chrono::duration<double, std::milli>(end - start).count());

		// Clean up adjacency list
		generator.deleteAdjacencyList(adjacencyList);

		// Generate and test incidence matrix
		auto incidenceMatrix = generator.generateWeightedIncidenceMatrix(simData.size);
		incidenceMatrix = generator.reduceIncidenceMatrixDensity(incidenceMatrix, simData.density, simData.size);

		// Time Prim's algorithm on incidence matrix
		start = std::chrono::high_resolution_clock::now();
		auto primResultIM = algorithms.primAlgorithmIM(incidenceMatrix, simData.size);
		end = std::chrono::high_resolution_clock::now();
		times.primTimeIM.push_back(std::chrono::duration<double, std::milli>(end - start).count());

		// Time Kruskal's algorithm on incidence matrix
		start = std::chrono::high_resolution_clock::now();
		auto kruskalResultIM = algorithms.kruskalAlgorithmIM(incidenceMatrix, simData.size);
		end = std::chrono::high_resolution_clock::now();
		times.kruskalTimeIM.push_back(std::chrono::duration<double, std::milli>(end - start).count());

		// Incidence matrices are typically managed automatically (vector of vectors)
		// No explicit cleanup needed unless using dynamic allocation
	}

	return results;
}
std::vector<spResultsData> performSPSimulation(mstSimData simData) {
	// Create one result structure to store all times for this size
	std::vector<spResultsData> results(1);  // Only one element needed
	spResultsData& times = results[0];

	Generator generator(simData.maxRand);
	Algorithms algorithms;

	// Reserve space for better performance
	times.dijkstraTimeAL.reserve(simData.instanceAmount);
	times.bellmanFordTimeAL.reserve(simData.instanceAmount);
	times.dijkstraTimeIM.reserve(simData.instanceAmount);
	times.bellmanFordTimeIM.reserve(simData.instanceAmount);

	// Process each instance individually to save memory
	for (int i = 0; i < simData.instanceAmount; ++i) {
		// Generate and test directed adjacency list
		auto adjacencyList = generator.generateWeightedDirectedAdjacencyList(simData.size);
		adjacencyList = generator.reduceDirectedAdjacencyListDensity(adjacencyList, simData.density, simData.size);

		// Time Dijkstra's algorithm on adjacency list
		auto start = std::chrono::high_resolution_clock::now();
		auto dijkstraResultAL = algorithms.dijkstraAlgorithmAL(adjacencyList, simData.size);
		auto end = std::chrono::high_resolution_clock::now();
		times.dijkstraTimeAL.push_back(std::chrono::duration<double, std::milli>(end - start).count());

		// Time Bellman-Ford algorithm on adjacency list
		start = std::chrono::high_resolution_clock::now();
		auto bellmanFordResultAL = algorithms.bellmanFordAlgorithmAL(adjacencyList, simData.size);
		end = std::chrono::high_resolution_clock::now();
		times.bellmanFordTimeAL.push_back(std::chrono::duration<double, std::milli>(end - start).count());

		// Clean up adjacency list
		generator.deleteAdjacencyList(adjacencyList);

		// Generate and test directed incidence matrix
		auto incidenceMatrix = generator.generateWeightedDirectedIncidenceMatrix(simData.size);
		incidenceMatrix = generator.reduceIncidenceMatrixDensity(incidenceMatrix, simData.density, simData.size);

		// Time Dijkstra's algorithm on incidence matrix
		start = std::chrono::high_resolution_clock::now();
		auto dijkstraResultIM = algorithms.dijkstraAlgorithmIM(incidenceMatrix, simData.size);
		end = std::chrono::high_resolution_clock::now();
		times.dijkstraTimeIM.push_back(std::chrono::duration<double, std::milli>(end - start).count());

		// Time Bellman-Ford algorithm on incidence matrix
		start = std::chrono::high_resolution_clock::now();
		auto bellmanFordResultIM = algorithms.bellmanFordAlgorithmIM(incidenceMatrix, simData.size);
		end = std::chrono::high_resolution_clock::now();
		times.bellmanFordTimeIM.push_back(std::chrono::duration<double, std::milli>(end - start).count());

		// Incidence matrices are typically managed automatically
	}

	return results;
}
void performSim(const fileData fileData) {
	std::vector<algorithmAvgTimes> avgTimes;
	mstSimData simData;
	simData.size = fileData.initialSize;
	simData.instanceAmount = fileData.instanceAmount;
	simData.density = fileData.density;
	simData.maxRand = fileData.maxRand;
	algorithmAvgTimes avgTimesData;
	std::cout << "Starting simulation with the following parameters:\n";
	std::cout << "Initial size: " << simData.size << "\n";
	std::cout << "Size amount: " << fileData.sizeAmount << "\n";
	std::cout << "Density: " << simData.density*100 << "%\n";
	std::cout << "Max random weight: " << simData.maxRand << "\n";
	std::cout << "Interval amount: " << fileData.intervalAmount << "\n";
	std::cout << "Instance amount: " << simData.instanceAmount << "\n";
	for (int i = 0; i < fileData.sizeAmount; i++) {
		simData.size = fileData.initialSize + i * fileData.intervalAmount;
		std::cout << "Graph size: " << simData.size << "\n";
		//MST sim
		auto mstResults = performMSTSimulation(simData);
		for (const auto& result : mstResults) {
			std::cout << "Average Prim's time on adjacency list:         " << calculateAverageTime(result.primTimeAL) << " ms\n";
			std::cout << "Average Kruskal's time on adjacency list:      " << calculateAverageTime(result.kruskalTimeAL) << " ms\n";
			std::cout << "Average Prim's time on incidence matrix:       " << calculateAverageTime(result.primTimeIM) << " ms\n";
			std::cout << "Average Kruskal's time on incidence matrix:    " << calculateAverageTime(result.kruskalTimeIM) << " ms\n";
		}
		//SP sim
		auto spResults = performSPSimulation(simData);
		for (const auto& result : spResults) {
			std::cout << "Average Dijkstra's time on adjacency list:     " << calculateAverageTime(result.dijkstraTimeAL) << " ms\n";
			std::cout << "Average Bellman-Ford time on adjacency list:   " << calculateAverageTime(result.bellmanFordTimeAL) << " ms\n";
			std::cout << "Average Dijkstra's time on incidence matrix:   " << calculateAverageTime(result.dijkstraTimeIM) << " ms\n";
			std::cout << "Average Bellman-Ford time on incidence matrix: " << calculateAverageTime(result.bellmanFordTimeIM) << " ms\n";
		}
		//Calculate averages
		avgTimesData.primAL = calculateAverageTime(mstResults[0].primTimeAL);
		avgTimesData.kruskalAL = calculateAverageTime(mstResults[0].kruskalTimeAL);
		avgTimesData.primIM = calculateAverageTime(mstResults[0].primTimeIM);
		avgTimesData.kruskalIM = calculateAverageTime(mstResults[0].kruskalTimeIM);
		avgTimesData.dijkstraAL = calculateAverageTime(spResults[0].dijkstraTimeAL);
		avgTimesData.bellmanFordAL = calculateAverageTime(spResults[0].bellmanFordTimeAL);
		avgTimesData.dijkstraIM = calculateAverageTime(spResults[0].dijkstraTimeIM);
		avgTimesData.bellmanFordIM = calculateAverageTime(spResults[0].bellmanFordTimeIM);
		//Save averages from this size
		avgTimes.push_back(avgTimesData);
		std::cout << "----------------------------------------\n";
	}
	saveToCSV(avgTimes);
}
int main() {
	fileData fileData = loadFileData();
	if (fileData.mode == 1) {
		performSim(fileData);
	} else if (fileData.mode == 0) {
		performTest(fileData);
	} else {
		std::cerr << "Invalid mode in config file." << std::endl;
		return 1;
	}
	return 0;
}