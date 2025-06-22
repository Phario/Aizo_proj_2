#pragma once
#ifndef ALGORITHMS_H
#define ALGORITHMS_H
#include<ctime>
#include<cstdlib>
#include<vector>
#include<string>
#include <algorithm>
#include "Generator.h"
#include <iostream>
#include <queue>
// structure to hold the results of pathfinding algorithms
struct PathResult {
	std::vector<int> cost;           // cost[i] = shortest distance from source to vertex i
	std::vector<int> predecessors;        // predecessors[i] = previous vertex in shortest path to i
	std::vector<std::vector<int>> paths;  // paths[i] = complete path from source to vertex i
	int source;
};
//structure to hold the results of MST algorithms
struct MSTResult {
	std::vector<std::string> edges;
	int cost;
};
// Utility functions for Union-Find data structure
struct Edge {
	int u;
	int v;
	int weight;
	bool operator<(const Edge& other) const {
		return weight < other.weight;
	}
};
struct DijkstraNode {
	int vertex;
	int distance;
	// For min heap
	bool operator>(const DijkstraNode& other) const {
		return distance > other.distance;
	}
};
class Algorithms {
public:
	// Prim's algorithm using adjacency list
	MSTResult primAlgorithmAL(const std::vector<neighbour*>& adjacencyList, int vertices) {
		std::vector<bool> inMST(vertices, false);
		std::vector<int> minEdgeWeight(vertices, INT_MAX);
		std::vector<int> parent(vertices, -1);
		minEdgeWeight[0] = 0; // Start from the first vertex
		for (int i = 0; i < vertices - 1; i++) {
			int u = -1;
			for (int j = 0; j < vertices; j++) {
				if (!inMST[j] && (u == -1 || minEdgeWeight[j] < minEdgeWeight[u])) {
					u = j;
				}
			}
			inMST[u] = true;
			for (neighbour* v = adjacencyList[u]; v != nullptr; v = v->next) {
				if (!inMST[v->vertex] && v->weight < minEdgeWeight[v->vertex]) {
					minEdgeWeight[v->vertex] = v->weight;
					parent[v->vertex] = u;
				}
			}
		}
		MSTResult result;
		result.cost = 0;
		for (int i = 1; i < vertices; i++) {
			result.edges.push_back(std::to_string(parent[i]) + " - " + std::to_string(i) + " (" + std::to_string(minEdgeWeight[i]) + ")");
			result.cost += minEdgeWeight[i];
		}
		return result;
	}
	// Prim's algorithm using incidence matrix
	MSTResult primAlgorithmIM(const std::vector<std::vector<int>>& incidenceMatrix, int vertices) {
		std::vector<bool> inMST(vertices, false);
		std::vector<int> minEdgeWeight(vertices, INT_MAX);
		std::vector<int> parent(vertices, -1);
		minEdgeWeight[0] = 0;
		int numEdges = incidenceMatrix[0].size();
		for (int i = 0; i < vertices - 1; i++) {
			int u = -1;
			// Find the vertex with minimum edge weight that's not in MST
			for (int j = 0; j < vertices; j++) {
				if (!inMST[j] && (u == -1 || minEdgeWeight[j] < minEdgeWeight[u])) {
					u = j;
				}
			}
			inMST[u] = true;
			// Check all edges to find adjacent vertices
			for (int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++) {
				// Check if vertex u participates in this edge
				if (incidenceMatrix[u][edgeIdx] != 0) {
					// Find the other vertex in this edge
					for (int v = 0; v < vertices; v++) {
						if (v != u && incidenceMatrix[v][edgeIdx] != 0) {
							// Found adjacent vertex v
							if (!inMST[v]) {
								// For weighted incidence matrix, the weight is the absolute value
								int weight = abs(incidenceMatrix[u][edgeIdx]);
								// Check if this edge is better than the current minimum for v
								if (weight < minEdgeWeight[v]) {
									minEdgeWeight[v] = weight;
									parent[v] = u;
								}
							}
							break; // Each edge connects exactly 2 vertices
						}
					}
				}
			}
		}
		// Build result
		MSTResult result;
		result.cost = 0;
		for (int i = 1; i < vertices; i++) {
			result.edges.push_back(std::to_string(parent[i]) + " - " + std::to_string(i) + " (" + std::to_string(minEdgeWeight[i]) + ")");
			result.cost += minEdgeWeight[i];
		}
		return result;
	}
	// Kruskal's algorithm using adjacency list
	MSTResult kruskalAlgorithmAL(const std::vector<neighbour*>& adjacencyList, int vertices) {
		std::vector<Edge> edges;
		// Extract all edges from adjacency list
		for (int u = 0; u < vertices; u++) {
			for (neighbour* current = adjacencyList[u]; current != nullptr; current = current->next) {
				int v = current->vertex;
				int weight = current->weight;
				// Only add each edge once (u < v to avoid duplicates)
				if (u < v) {
					edges.push_back({ u, v, weight });
				}
			}
		}
		// Sort edges by weight
		std::sort(edges.begin(), edges.end());
		// Initialize Union find data structure
		std::vector<int> parent(vertices);
		std::vector<int> rank(vertices, 0);
		for (int i = 0; i < vertices; i++) {
			parent[i] = i;
		}
		MSTResult result;
		result.cost = 0;
		// Process edges in order of increasing weight
		for (const Edge& edge : edges) {
			int rootU = find(parent, edge.u);
			int rootV = find(parent, edge.v);
			// If adding this edge doesn't create a cycle
			if (rootU != rootV) {
				unionSets(parent, rank, edge.u, edge.v);
				// Add edge to result
				result.edges.push_back(std::to_string(edge.u) + " - " + std::to_string(edge.v) + " (" + std::to_string(edge.weight) + ")");
				result.cost += edge.weight;
				// MST is complete when we have vertices-1 edges
				if (result.edges.size() == vertices - 1) {
					break;
				}
			}
		}

		return result;
	}
	// Kruskal's algorithm using incidence matrix
	MSTResult kruskalAlgorithmIM(const std::vector<std::vector<int>>& incidenceMatrix, int vertices) {
		std::vector<Edge> edges;
		int numEdges = incidenceMatrix[0].size();
		// Extract all edges from incidence matrix
		for (int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++) {
			int u = -1, v = -1;
			int weight = 0;
			// Find the two vertices connected by this edge
			for (int vertex = 0; vertex < vertices; vertex++) {
				if (incidenceMatrix[vertex][edgeIdx] != 0) {
					weight = abs(incidenceMatrix[vertex][edgeIdx]);
					if (u == -1) {
						u = vertex;
					}
					else {
						v = vertex;
						break; // Found both vertices
					}
				}
			}
			// Add edge if we found exactly two vertices
			if (u != -1 && v != -1) {
				// Ensure consistent ordering (u < v) to match adjacency list behavior
				if (u > v) {
					std::swap(u, v);
				}
				edges.push_back({ u, v, weight });
			}
		}
		// Sort edges by weight
		std::sort(edges.begin(), edges.end());
		// Initialize Union-Find data structure
		std::vector<int> parent(vertices);
		std::vector<int> rank(vertices, 0);
		for (int i = 0; i < vertices; i++) {
			parent[i] = i;
		}

		MSTResult result;
		result.cost = 0;
		// Process edges in order of increasing weight
		for (const Edge& edge : edges) {
			int rootU = find(parent, edge.u);
			int rootV = find(parent, edge.v);
			// If adding this edge doesn't create a cycle
			if (rootU != rootV) {
				unionSets(parent, rank, edge.u, edge.v);
				result.edges.push_back(std::to_string(edge.u) + " - " + std::to_string(edge.v) + " (" + std::to_string(edge.weight) + ")");
				result.cost += edge.weight;
				// MST is complete when we have vertices-1 edges
				if (result.edges.size() == vertices - 1) {
					break;
				}
			}
		}

		return result;
	}
	// Dijkstra's algorithm using adjacency list
	PathResult dijkstraAlgorithmAL(const std::vector<neighbour*>& adjacencyList, int vertices) {
		int source = 0;
		std::vector<int> cost(vertices, INT_MAX);
		std::vector<int> prev(vertices, -1);
		std::vector<bool> visited(vertices, false);
		// Priority queue
		std::priority_queue<DijkstraNode, std::vector<DijkstraNode>, std::greater<DijkstraNode>> pq;
		cost[source] = 0;
		pq.push({ source, 0 });
		// Process the priority queue
		while (!pq.empty()) {
			// Get the node with the smallest distance
			DijkstraNode current = pq.top();
			pq.pop();
			int u = current.vertex;
			// Skip if already visited (handles duplicate entries in pq)
			if (visited[u]) continue;
			visited[u] = true;
			// Relaxation
			for (neighbour* v = adjacencyList[u]; v != nullptr; v = v->next) {
				int vertex = v->vertex;
				int weight = v->weight;
				// Only process unvisited vertices and check for overflow
				if (!visited[vertex] && cost[u] != INT_MAX && cost[u] + weight < cost[vertex]) {
					cost[vertex] = cost[u] + weight;
					prev[vertex] = u;
					pq.push({ vertex, cost[vertex] });
				}
			}
		}
		// Build result structure
		PathResult result;
		result.source = source;
		result.cost = std::move(cost);
		result.predecessors = std::move(prev);
		result.paths.resize(vertices);
		// Reconstruct all paths
		for (int target = 0; target < vertices; target++) {
			if (result.cost[target] == INT_MAX) {
				// No path to this vertex - paths[target] is already empty
				continue;
			}
			// Reconstruct path from source to target
			std::vector<int> path;
			int current = target;
			while (current != -1) {
				path.push_back(current);
				current = result.predecessors[current];
			}
			std::reverse(path.begin(), path.end());
			result.paths[target] = std::move(path);
		}
		return result;
	}
	// Dijkstra's algorithm using incidence matrix
	PathResult dijkstraAlgorithmIM(const std::vector<std::vector<int>>& incidenceMatrix, int vertices) {
		int source = 0; // Always start from vertex 0
		std::vector<int> cost(vertices, INT_MAX);
		std::vector<int> prev(vertices, -1);
		std::vector<bool> visited(vertices, false);
		int numEdges = incidenceMatrix[0].size();
		// Convert O(V*E) edge lookups to O(V) adjacency lookups
		std::vector<std::vector<std::pair<int, int>>> adjacencyFromIM(vertices);
		for (int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++) {
			int u = -1, v = -1;
			int weight = 0;
			// Find the two vertices connected by this edge
			for (int vertex = 0; vertex < vertices; vertex++) {
				if (incidenceMatrix[vertex][edgeIdx] != 0) {
					if (u == -1) {
						u = vertex;
						weight = abs(incidenceMatrix[vertex][edgeIdx]);
					}
					else {
						v = vertex;
						break;
					}
				}
			}
			// Add bidirectional edges to adjacency representation
			if (u != -1 && v != -1) {
				adjacencyFromIM[u].emplace_back(v, weight);
				adjacencyFromIM[v].emplace_back(u, weight);
			}
		}
		std::priority_queue<DijkstraNode, std::vector<DijkstraNode>, std::greater<DijkstraNode>> pq;
		cost[source] = 0;
		pq.push({ source, 0 });
		while (!pq.empty()) {
			DijkstraNode current = pq.top();
			pq.pop();
			int u = current.vertex;
			// Skip if already visited
			if (visited[u]) continue;
			visited[u] = true;
			// Process all adjacent vertices using pre-computed adjacency
			for (const auto& edge : adjacencyFromIM[u]) {
				int vertex = edge.first;
				int weight = edge.second;
				if (!visited[vertex] && cost[u] != INT_MAX && cost[u] + weight < cost[vertex]) {
					cost[vertex] = cost[u] + weight;
					prev[vertex] = u;
					pq.push({ vertex, cost[vertex] });
				}
			}
		}
		// Build result structure
		PathResult result;
		result.source = source;
		result.cost = std::move(cost);
		result.predecessors = std::move(prev);
		result.paths.resize(vertices);
		// Reconstruct all paths
		for (int target = 0; target < vertices; target++) {
			if (result.cost[target] == INT_MAX) {
				// No path to this vertex - paths[target] is already empty
				continue;
			}
			// Reconstruct path from source to target
			std::vector<int> path;
			int current = target;
			while (current != -1) {
				path.push_back(current);
				current = result.predecessors[current];
			}
			std::reverse(path.begin(), path.end());
			result.paths[target] = std::move(path);
		}
		return result;
	}	
	// Bellman-Ford algorithm using adjacency list
	PathResult bellmanFordAlgorithmAL(const std::vector<neighbour*>& adjacencyList, int vertices) {
		int source = 0;
		std::vector<int> cost(vertices, INT_MAX);
		std::vector<int> prev(vertices, -1);
		cost[source] = 0;
		// Relax edges v-1 times
		for (int i = 0; i < vertices - 1; i++) {
			bool updated = false;
			for (int u = 0; u < vertices; u++) {
				if (cost[u] != INT_MAX) {
					for (neighbour* v = adjacencyList[u]; v != nullptr; v = v->next) {
						if (cost[u] + v->weight < cost[v->vertex]) {
							cost[v->vertex] = cost[u] + v->weight;
							prev[v->vertex] = u;
							updated = true;
						}
					}
				}
			}
			// Early termination
			if (!updated) break;
		}
		// Build result structure
		PathResult result;
		result.source = source;
		result.cost = cost;
		result.predecessors = prev;
		result.paths.resize(vertices);
		// Reconstruct all paths
		for (int target = 0; target < vertices; target++) {
			if (cost[target] == INT_MAX) {
				// No path to this vertex
				result.paths[target] = {}; // Empty path
			}
			else {
				// Reconstruct path from source to target
				std::vector<int> path;
				int current = target;
				while (current != -1) {
					path.push_back(current);
					current = prev[current];
				}
				std::reverse(path.begin(), path.end());
				result.paths[target] = path;
			}
		}
		return result;
	}
	// Bellman-Ford algorithm using incidence matrix
	PathResult bellmanFordAlgorithmIM(const std::vector<std::vector<int>>& incidenceMatrix, int vertices) {
		int source = 0;
		std::vector<int> cost(vertices, INT_MAX);
		std::vector<int> prev(vertices, -1);
		int numEdges = incidenceMatrix[0].size();
		cost[source] = 0;
		// Relax edges vertices-1 times
		for (int i = 0; i < vertices - 1; i++) {
			bool updated = false;
			// Check all edges
			for (int edgeIdx = 0; edgeIdx < numEdges; edgeIdx++) {
				int u = -1, v = -1;
				int weight = 0;
				// Find the two vertices connected by this edge
				for (int vertex = 0; vertex < vertices; vertex++) {
					if (incidenceMatrix[vertex][edgeIdx] != 0) {
						weight = abs(incidenceMatrix[vertex][edgeIdx]);
						if (u == -1) {
							u = vertex;
						}
						else {
							v = vertex;
							break; // Found both vertices
						}
					}
				}
				// Relax edge in both directions
				if (u != -1 && v != -1) {
					// Edge u -> v
					if (cost[u] != INT_MAX && cost[u] + weight < cost[v]) {
						cost[v] = cost[u] + weight;
						prev[v] = u;
						updated = true;
					}
					// Edge v -> u
					if (cost[v] != INT_MAX && cost[v] + weight < cost[u]) {
						cost[u] = cost[v] + weight;
						prev[u] = v;
						updated = true;
					}
				}
			}
			// Early termination if no updates occurred
			if (!updated) break;
		}
		// Build result structure
		PathResult result;
		result.source = source;
		result.cost = cost;
		result.predecessors = prev;
		result.paths.resize(vertices);
		// Reconstruct all paths
		for (int target = 0; target < vertices; target++) {
			if (cost[target] == INT_MAX) {
				// No path to this vertex
				result.paths[target] = {}; // Empty path
			}
			else {
				// Reconstruct path from source to target
				std::vector<int> path;
				int current = target;
				while (current != -1) {
					path.push_back(current);
					current = prev[current];
				}
				// Reverse the path to get it from source to target
				std::reverse(path.begin(), path.end());
				result.paths[target] = path;
			}
		}
		return result;
	}
	// Utility functions to print results
	void printMST(const MSTResult& result) {
		std::cout << "Minimum Spanning Tree Edges:\n";
		for (const auto& edge : result.edges) {
			std::cout << edge << "\n";
		}
		std::cout << "Total Cost: " << result.cost << "\n";
	}
	void printPaths(const PathResult& result) {
		std::cout << "Shortest Paths from Source " << result.source << ":\n";
		for (int i = 0; i < result.paths.size(); i++) {
			std::cout << "To Vertex " << i << ": Cost = " << result.cost[i] << ", Path = ";
			if (result.paths[i].empty()) {
				std::cout << "No path\n";
			} else {
				for (int j = 0; j < result.paths[i].size(); j++) {
					std::cout << result.paths[i][j];
					if (j < result.paths[i].size() - 1) std::cout << " -> ";
				}
				std::cout << "\n";
			}
		}
	}
private:
		// Union-Find (Disjoint Set) helper functions for Kruskal's algorithm
		int find(std::vector<int>& parent, int x) {
			if (parent[x] != x) {
				parent[x] = find(parent, parent[x]); // Path compression
			}
			return parent[x];
		}
		void unionSets(std::vector<int>& parent, std::vector<int>& rank, int x, int y) {
			int rootX = find(parent, x);
			int rootY = find(parent, y);
			if (rootX != rootY) {
				// Union by rank
				if (rank[rootX] < rank[rootY]) {
					parent[rootX] = rootY;
				}
				else if (rank[rootX] > rank[rootY]) {
					parent[rootY] = rootX;
				}
				else {
					parent[rootY] = rootX;
					rank[rootX]++;
				}
			}
		}

};
#endif