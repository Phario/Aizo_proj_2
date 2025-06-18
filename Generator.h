#pragma once
#ifndef GENERATOR_H
#define GENERATOR_H
#include<ctime>
#include<cstdlib>
#include<vector>
#include<string>
#include <algorithm>
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
	std::vector<std::vector<int>> reduceIncidenceMatrixDensity(const std::vector<std::vector<int>>& matrix, double density, int vertices) {

		int originalCols = matrix[0].size();
		int newEdges = static_cast<int>(originalCols * density);
		//safety check
		if (newEdges > originalCols) {
			newEdges = originalCols;
		}
		std::vector<std::vector<int>> reducedMatrix(vertices, std::vector<int>(newEdges, 0));
		int reducedEdgeIndex = 0;
		for (int col = 0; col < originalCols && reducedEdgeIndex < newEdges; col++) {
			int nonZeroCount = 0;
			for (int row = 0; row < vertices; row++) {
				if (matrix[row][col] != 0) {
					nonZeroCount++;
				}
			}
			if (nonZeroCount == 2) {
				for (int row = 0; row < vertices; row++) {
					reducedMatrix[row][reducedEdgeIndex] = matrix[row][col];
				}
				reducedEdgeIndex++;
			}
		}

		return reducedMatrix;
	}
	std::vector<neighbour*> reduceAdjacencyListDensity(const std::vector<neighbour*>& list, double density, int vertices) {
		//count edges
		int totalEdgeCount = 0;
		for (int i = 0; i < vertices; i++) {
			neighbour* current = list[i];
			while (current) {
				totalEdgeCount++;
				current = current->next;
			}
		}

		// For undirected graphs, divide by 2 to get actual edge count
		int actualEdges = totalEdgeCount / 2;
		int newEdges = static_cast<int>(actualEdges * density);

		std::vector<neighbour*> reducedList(vertices, nullptr);

		int edgesAdded = 0;

		//iterate through vertices and their neighbors
		for (int i = 0; i < vertices && edgesAdded < newEdges; i++) {
			neighbour* current = list[i];
			while (current && edgesAdded < newEdges) {
				int j = current->vertex;
				int weight = current->weight;

				//For undirected graphs, only process each edge once (i < j)
				//This avoids duplicates since each edge appears in both adjacency lists
				if (i < j) {
					// Add edge i-j to vertex i's list
					neighbour* newNeighbour1 = new neighbour{ reducedList[i], j, weight };
					reducedList[i] = newNeighbour1;

					// Add edge j-i to vertex j's list
					neighbour* newNeighbour2 = new neighbour{ reducedList[j], i, weight };
					reducedList[j] = newNeighbour2;

					edgesAdded++;
				}

				current = current->next;
			}
		}

		return reducedList;
	}
	std::vector<neighbour*> reduceDirectedAdjacencyListDensity(const std::vector<neighbour*>& list, double density, int vertices) {
		// Count total directed edges
		int totalEdges = 0;
		for (int i = 0; i < vertices; i++) {
			neighbour* current = list[i];
			while (current) {
				totalEdges++;
				current = current->next;
			}
		}

		int newEdges = static_cast<int>(totalEdges * density);

		// Create reduced adjacency list
		std::vector<neighbour*> reducedList(vertices, nullptr);

		int edgesAdded = 0;

		// Iterate through all vertices and their neighbors
		for (int i = 0; i < vertices && edgesAdded < newEdges; i++) {
			neighbour* current = list[i];
			while (current && edgesAdded < newEdges) {
				int j = current->vertex;
				int weight = current->weight;

				// For directed graphs, add each edge as-is
				neighbour* newNeighbour = new neighbour{ reducedList[i], j, weight };
				reducedList[i] = newNeighbour;

				edgesAdded++;
				current = current->next;
			}
		}

		return reducedList;
	}
};
#endif