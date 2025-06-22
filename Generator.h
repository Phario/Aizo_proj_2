#pragma once
#ifndef GENERATOR_H
#define GENERATOR_H
#include <ctime>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <tuple>
#include <random>
#include <unordered_set> 
#include <fstream>

/*
Structure of a source file for a graph:
E V						# E - number of edges, V - number of vertices separated by space, the edges are numbered from 0 to V-1
start end weight		# start and end vertices of the edge and its weight, one edge per line
...						# MST - edges are undirected, SP - edges are directed
*/
// Structure for a neighbour in the adjacency list
struct neighbour {
	neighbour* next;
	int vertex;
	int weight;
};

class Generator {

public:
	int maxRand; // Maximum weight for edges
	Generator(int maxRand) : maxRand(maxRand){// Maximum weight for edges
		srand(time(NULL)); 
	}
	// Generating graphs from scratch
	// Generate a weighted incidence matrix for an undirected graph
	std::vector<std::vector<int>> generateWeightedIncidenceMatrix(int n) {
		std::vector<std::vector<int>> incidenceMatrix(n, std::vector<int>(n * (n - 1) / 2, 0));
		int edgeIndex = 0;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				int weight = rand() % maxRand + 1;
				incidenceMatrix[i][edgeIndex] = weight;
				incidenceMatrix[j][edgeIndex] = weight;
				edgeIndex++;
			}
		}
		return incidenceMatrix;
	}
	// Generate a weighted incidence matrix for a directed graph
	std::vector<std::vector<int>> generateWeightedDirectedIncidenceMatrix(int n) {
		std::vector<std::vector<int>> incidenceMatrix(n, std::vector<int>(n * (n - 1), 0));
		int edgeIndex = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					int weight = rand() % maxRand + 1;
					incidenceMatrix[i][edgeIndex] = weight;
					incidenceMatrix[j][edgeIndex] = -weight;
					edgeIndex++;
				}
			}
		}
		return incidenceMatrix;
	}
	// Generate a weighted adjacency list for an undirected graph
	std::vector<neighbour*> generateWeightedAdjacencyList(int n) {
		std::vector<neighbour*> adjacencyList(n, nullptr);

		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				int weight = rand() % maxRand + 1;
				// Add j as neighbor of i
				neighbour* newNeighbour1 = new neighbour{ adjacencyList[i], j, weight };
				adjacencyList[i] = newNeighbour1;
				// Add i as neighbor of j (undirected graph)
				neighbour* newNeighbour2 = new neighbour{ adjacencyList[j], i, weight };
				adjacencyList[j] = newNeighbour2;
			}
		}

		return adjacencyList;
	}
	// Generate a weighted incidence matrix from a file
	std::vector<neighbour*> generateWeightedDirectedAdjacencyList(int n) {
		std::vector<neighbour*> adjacencyList(n, nullptr);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					int weight = rand() % maxRand + 1;
					// Add j as neighbor of i (directed edge i -> j)
					neighbour* newNeighbour = new neighbour{ adjacencyList[i], j, weight };
					adjacencyList[i] = newNeighbour;
				}
			}
		}
		return adjacencyList;
	}
	// Generate a weighted incidence matrix from a file
	std::vector<std::vector<int>> reduceIncidenceMatrixDensity(const std::vector<std::vector<int>>& matrix, double density, int vertices) {
		int originalCols = matrix[0].size();
		int newEdges = static_cast<int>(originalCols * density);
		// Safety check
		if (newEdges > originalCols) {
			newEdges = originalCols;
		}

		// Find all valid edge columns
		std::vector<int> validColumns;
		for (int col = 0; col < originalCols; col++) {
			int nonZeroCount = 0;
			for (int row = 0; row < vertices; row++) {
				if (matrix[row][col] != 0) {
					nonZeroCount++;
				}
			}
			if (nonZeroCount == 2) {
				validColumns.push_back(col);
			}
		}
		// Randomly shuffle the valid columns
		std::random_device rd;
		std::mt19937 gen(rd());
		std::shuffle(validColumns.begin(), validColumns.end(), gen);
		// Take the first of the newEdges columns from the shuffled list
		int edgesToKeep = std::min(newEdges, static_cast<int>(validColumns.size()));
		std::vector<std::vector<int>> reducedMatrix(vertices, std::vector<int>(edgesToKeep, 0));
		// Fill the reduced matrix with the selected edges
		for (int i = 0; i < edgesToKeep; i++) {
			int originalCol = validColumns[i];
			for (int row = 0; row < vertices; row++) {
				reducedMatrix[row][i] = matrix[row][originalCol];
			}
		}
		return reducedMatrix;
	}
	// Reduce the density of an undirected adjacency list
	std::vector<neighbour*> reduceAdjacencyListDensity(const std::vector<neighbour*>& list, double density, int vertices) {
		// Collect all unique undirected edges (i < j)
		std::vector<std::tuple<int, int, int>> uniqueEdges;
		std::unordered_set<int> originallyConnected;
		for (int i = 0; i < vertices; ++i) {
			neighbour* current = list[i];
			while (current) {
				int j = current->vertex;
				int weight = current->weight;
				originallyConnected.insert(i);
				originallyConnected.insert(j);
				if (i < j) {
					uniqueEdges.emplace_back(i, j, weight);
				}
				current = current->next;
			}
		}
		// Shuffle edges to randomize
		std::random_device rd;
		std::mt19937 gen(rd());
		std::shuffle(uniqueEdges.begin(), uniqueEdges.end(), gen);
		// Calculate number of edges to keep based on density
		int totalEdges = static_cast<int>(uniqueEdges.size());
		int edgesToKeep = static_cast<int>(totalEdges * density);
		// Ensure minimum connectivity
		int minEdgesForConnectivity = std::max(1, static_cast<int>(originallyConnected.size()) - 1);
		edgesToKeep = std::max(edgesToKeep, minEdgesForConnectivity);
		std::vector<neighbour*> reducedList(vertices, nullptr);
		std::unordered_set<int> connectedVertices;
		int edgesAdded = 0;
		// First pass: prioritize connecting isolated vertices
		for (const auto& edge : uniqueEdges) {
			if (edgesAdded >= edgesToKeep) break;
			int u = std::get<0>(edge);
			int v = std::get<1>(edge);
			int w = std::get<2>(edge);
			// Add edge if it connects a new vertex or limit isn't reached
			bool connectsNewVertex = (connectedVertices.count(u) == 0 || connectedVertices.count(v) == 0);
			if (connectsNewVertex || edgesAdded < edgesToKeep) {
				// Add u -> v
				reducedList[u] = new neighbour{ reducedList[u], v, w };
				// Add v -> u  
				reducedList[v] = new neighbour{ reducedList[v], u, w };
				connectedVertices.insert(u);
				connectedVertices.insert(v);
				edgesAdded++;
			}
		}
		return reducedList;
	}
	// Reduce the density of a directed adjacency list
	std::vector<neighbour*> reduceDirectedAdjacencyListDensity(const std::vector<neighbour*>& list, double density, int vertices) {
		// Collect all directed edges
		std::vector<std::tuple<int, int, int>> allEdges; // (from, to, weight)
		for (int i = 0; i < vertices; i++) {
			neighbour* current = list[i];
			while (current) {
				allEdges.emplace_back(i, current->vertex, current->weight);
				current = current->next;
			}
		}
		// Shuffle edges randomly
		std::random_device rd;
		std::mt19937 gen(rd());
		std::shuffle(allEdges.begin(), allEdges.end(), gen);
		// Calculate how many edges to keep
		int totalEdges = static_cast<int>(allEdges.size());
		int edgesToKeep = static_cast<int>(totalEdges * density);
		edgesToKeep = std::max(edgesToKeep, 1); // Keep at least 1 edge if any exist
		// Create new adjacency list with selected edges
		std::vector<neighbour*> reducedList(vertices, nullptr);
		for (int i = 0; i < edgesToKeep && i < totalEdges; i++) {
			int from = std::get<0>(allEdges[i]);
			int to = std::get<1>(allEdges[i]);
			int weight = std::get<2>(allEdges[i]);
			// Add edge to adjacency list (insert at head for simplicity)
			reducedList[from] = new neighbour{ reducedList[from], to, weight };
		}
		return reducedList;
	}	//print adjacency list for debugging
	// Generating graphs from file
	// Generate a weighted incidence matrix from a file
	std::vector<std::vector<int>> generateWeightedIncidenceMatrixFromFile(const std::string& filename) {
		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cerr << "Error opening file: " << filename << std::endl;
			return {};
		}
		int edges, vertices;
		file >> edges >> vertices;
		std::vector<std::vector<int>> incidenceMatrix(vertices, std::vector<int>(edges, 0));
		for (int i = 0; i < edges; i++) {
			int start, end, weight;
			file >> start >> end >> weight;
			incidenceMatrix[start][i] = weight;
			incidenceMatrix[end][i] = weight; // For undirected graph
		}
		file.close();
		return incidenceMatrix;
	}
	// Generate a weighted adjacency list from a file
	std::vector<neighbour*> generateWeightedAdjacencyListFromFile(const std::string& filename) {
		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cerr << "Error opening file: " << filename << std::endl;
			return {};
		}
		int edges, vertices;
		file >> edges >> vertices;
		std::vector<neighbour*> adjacencyList(vertices, nullptr);
		for (int i = 0; i < edges; i++) {
			int start, end, weight;
			file >> start >> end >> weight;
			// Add edge to adjacency list (undirected)
			adjacencyList[start] = new neighbour{ adjacencyList[start], end, weight };
			adjacencyList[end] = new neighbour{ adjacencyList[end], start, weight };
		}
		file.close();
		return adjacencyList;
	}
	// Generate a weighted directed adjacency list from a file
	std::vector<neighbour*> generateWeightedDirectedAdjacencyListFromFile(const std::string& filename) {
		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cerr << "Error opening file: " << filename << std::endl;
			return {};
		}
		int edges, vertices;
		file >> edges >> vertices;
		std::vector<neighbour*> adjacencyList(vertices, nullptr);
		for (int i = 0; i < edges; i++) {
			int start, end, weight;
			file >> start >> end >> weight;
			// Add directed edge to adjacency list
			adjacencyList[start] = new neighbour{ adjacencyList[start], end, weight };
		}
		file.close();
		return adjacencyList;
	}
	// Generate a weighted directed incidence matrix from a file
	std::vector<std::vector<int>> generateWeightedDirectedIncidenceMatrixFromFile(const std::string& filename) {
		std::ifstream file(filename);
		if (!file.is_open()) {
			std::cerr << "Error opening file: " << filename << std::endl;
			return {};
		}
		int edges, vertices;
		file >> edges >> vertices;
		std::vector<std::vector<int>> incidenceMatrix(vertices, std::vector<int>(edges, 0));
		for (int i = 0; i < edges; i++) {
			int start, end, weight;
			file >> start >> end >> weight;
			incidenceMatrix[start][i] = weight; // Directed edge
			incidenceMatrix[end][i] = -weight; // Indicating direction
		}
		file.close();
		return incidenceMatrix;
	}
	// Utils
	// Print adjacency list in formatted table format for debugging
	void printAdjacencyList(const std::vector<neighbour*>& list, int vertices) {
		for (int i = 0; i < vertices; i++) {
			std::cout << "Vertex " << i << ": ";
			neighbour* current = list[i];
			while (current) {
				std::cout << "(" << current->vertex << ", " << current->weight << ") ";
				current = current->next;
			}
			std::cout << std::endl;
		}
	}
	// Print incidence matrix in formatted table format for debugging
	void printIncidenceMatrix(const std::vector<std::vector<int>>& matrix, int vertices) {
		int edges = matrix[0].size();
		std::cout << "Incidence Matrix (" << vertices << " vertices, " << edges << " edges):" << std::endl;

		for (int i = 0; i < vertices; i++) {
			for (int j = 0; j < edges; j++) {
				std::cout << std::setw(4) << matrix[i][j];
			}
			std::cout << std::endl;
		}
	}
	// Deallocate adjacency list
	void deleteAdjacencyList(std::vector<neighbour*>& list) {
		for (auto& head : list) {
			while (head) {
				neighbour* temp = head;
				head = head->next;
				delete temp;
			}
		}
		list.clear();
	}
};
#endif