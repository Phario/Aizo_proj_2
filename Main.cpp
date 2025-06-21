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
0		#representation: 0-incidence matrix, 1-adjacency list
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
	int density;
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
	int density;
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
std::vector<mstResultsData> performMSTSimulation(mstSimData mstSimData);
std::vector<spResultsData> performSPSimulation(mstSimData mstSimData);
void performTest();
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
				case 3: fileData.density = value; break;
				case 4: fileData.maxRand = value; break;
				case 5: fileData.intervalAmount = value; break;
				case 6: fileData.instanceAmount = value; break;
				case 7: fileData.representation = value; break;
				case 8: fileData.task = value; break;
				case 9: fileData.displayGraph = value; break;
				case 10: fileData.algorithm = value; break;
				case 11: fileData.testGraphSource = value; break;
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
	return sum / results.size();
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
void performTest() {
	int graphSize = 8;
	double density = 0.5;
	int maxRand = 20;
	Generator generator(maxRand);
	Algorithms algorithms;
	//Adjacency list tests (undirected graph)
	std::cout << "Weighted undirected adjacency list:\n";
	auto adjacencyList = generator.generateWeightedAdjacencyList(graphSize);
	generator.printAdjacencyList(adjacencyList, graphSize);
	std::cout << "Reduced weighted undirected adjacency list:\n";
	adjacencyList = generator.reduceAdjacencyListDensity(adjacencyList, density, graphSize);
	generator.printAdjacencyList(adjacencyList, graphSize);
	std::cout << "Prim's algorithm on adjacency list:\n";
	auto mstResult = algorithms.primAlgorithmAL(adjacencyList, graphSize);
	algorithms.printMST(mstResult);
	std::cout << "Kruskal's algorithm on adjacency list:\n";
	auto kruskalResult = algorithms.kruskalAlgorithmAL(adjacencyList, graphSize);
	algorithms.printMST(kruskalResult);
	std::cout << "Dijkstra's algorithm on adjacency list:\n";
	auto pathResult = algorithms.dijkstraAlgorithmAL(adjacencyList, graphSize);
	algorithms.printPaths(pathResult);
	std::cout << "Bellman-Ford algorithm on adjacency list:\n";
	pathResult = algorithms.bellmanFordAlgorithmAL(adjacencyList, graphSize);
	algorithms.printPaths(pathResult);

	//Adjacency list tests (directed graph)
	std::cout << "Weighted directed adjecency list:\n";
	adjacencyList = generator.generateWeightedDirectedAdjacencyList(graphSize);
	generator.printAdjacencyList(adjacencyList, graphSize);
	std::cout << "Reduced weighted directed adjacency list:\n";
	adjacencyList = generator.reduceDirectedAdjacencyListDensity(adjacencyList, density, graphSize);
	generator.printAdjacencyList(adjacencyList, graphSize);
	std::cout << "Dijkstra's algorithm on adjacency list:\n";
	pathResult = algorithms.dijkstraAlgorithmAL(adjacencyList, graphSize);
	algorithms.printPaths(pathResult);
	std::cout << "Bellman-Ford algorithm on adjacency list:\n";
	pathResult = algorithms.bellmanFordAlgorithmAL(adjacencyList, graphSize);
	algorithms.printPaths(pathResult);

	//Incidence matrix tests (undirected graph)
	std::cout << "Weighted incidence matrix:\n";
	auto incidenceMatrix = generator.generateWeightedIncidenceMatrix(graphSize);
	generator.printIncidenceMatrix(incidenceMatrix, graphSize);
	std::cout << "Reduced weighted incidence matrix:\n";
	incidenceMatrix = generator.reduceIncidenceMatrixDensity(incidenceMatrix, density, graphSize);
	generator.printIncidenceMatrix(incidenceMatrix, graphSize);
	std::cout << "Prim's algorithm on incidence matrix:\n";
	mstResult = algorithms.primAlgorithmIM(incidenceMatrix, graphSize);
	algorithms.printMST(mstResult);
	std::cout << "Kruskal's algorithm on incidence matrix:\n";
	kruskalResult = algorithms.kruskalAlgorithmIM(incidenceMatrix, graphSize);
	algorithms.printMST(kruskalResult);
	std::cout << "Dijkstra's algorithm on incidence matrix:\n";
	pathResult = algorithms.dijkstraAlgorithmIM(incidenceMatrix, graphSize);
	algorithms.printPaths(pathResult);
	std::cout << "Bellman-Ford algorithm on incidence matrix:\n";
	pathResult = algorithms.bellmanFordAlgorithmIM(incidenceMatrix, graphSize);
	algorithms.printPaths(pathResult);

	//Incidence matrix tests (directed graph)
	std::cout << "Weighted directed incidence matrix:\n";
	incidenceMatrix = generator.generateWeightedDirectedIncidenceMatrix(graphSize);
	generator.printIncidenceMatrix(incidenceMatrix, graphSize);
	std::cout << "Reduced weighted directed incidence matrix:\n";
	incidenceMatrix = generator.reduceIncidenceMatrixDensity(incidenceMatrix, density, graphSize);
	generator.printIncidenceMatrix(incidenceMatrix, graphSize);
	std::cout << "Dijkstra's algorithm on incidence matrix:\n";
	pathResult = algorithms.dijkstraAlgorithmIM(incidenceMatrix, graphSize);
	algorithms.printPaths(pathResult);
	std::cout << "Bellman-Ford algorithm on incidence matrix:\n";
	pathResult = algorithms.bellmanFordAlgorithmIM(incidenceMatrix, graphSize);
	algorithms.printPaths(pathResult);

	// Clean up dynamically allocated memory in adjacency list
	generator.deleteAdjacencyList(adjacencyList);
}
int main() {
	std::vector<double> avgTimes;
	fileData fileData = loadFileData();
	if (fileData.mode == 1) {
		// Perform simulation mode
		mstSimData simData;
		simData.size = fileData.initialSize;
		simData.instanceAmount = fileData.instanceAmount;
		simData.density = fileData.density;
		simData.maxRand = fileData.maxRand;
		std::cout << "Starting simulation with the following parameters:\n";
		std::cout << "Initial size: " << simData.size << "\n";
		std::cout << "Size amount: " << fileData.sizeAmount << "\n";
		std::cout << "Density: " << simData.density << "%\n";
		std::cout << "Max random weight: " << simData.maxRand << "\n";
		std::cout << "Interval amount: " << fileData.intervalAmount << "\n";
		std::cout << "Instance amount: " << simData.instanceAmount << "\n";
		for (int i = 0; i < fileData.sizeAmount; i++) {
			simData.size = fileData.initialSize + i * fileData.intervalAmount;
			std::cout << "Graph size: " << simData.size << "\n";
			// MST simulation
			auto mstResults = performMSTSimulation(simData);
			for (const auto& result : mstResults) {
				std::cout << "Average Prim's time on adjacency list:         " << calculateAverageTime(result.primTimeAL) << " ms\n";
				std::cout << "Average Kruskal's time on adjacency list:      " << calculateAverageTime(result.kruskalTimeAL) << " ms\n";
				std::cout << "Average Prim's time on incidence matrix:       " << calculateAverageTime(result.primTimeIM) << " ms\n";
				std::cout << "Average Kruskal's time on incidence matrix:    " << calculateAverageTime(result.kruskalTimeIM) << " ms\n";
			}
			// SP simulation
			auto spResults = performSPSimulation(simData);
			for (const auto& result : spResults) {
				std::cout << "Average Dijkstra's time on adjacency list:     " << calculateAverageTime(result.dijkstraTimeAL) << " ms\n";
				std::cout << "Average Bellman-Ford time on adjacency list:   " << calculateAverageTime(result.bellmanFordTimeAL) << " ms\n";
				std::cout << "Average Dijkstra's time on incidence matrix:   " << calculateAverageTime(result.dijkstraTimeIM) << " ms\n";
				std::cout << "Average Bellman-Ford time on incidence matrix: " << calculateAverageTime(result.bellmanFordTimeIM) << " ms\n";
			}
		}
	} else if (fileData.mode == 0) {
		// Perform test mode
		performTest();
	} else {
		std::cerr << "Invalid mode in config file." << std::endl;
		return 1;
	}
	return 0;
}