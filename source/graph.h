#ifndef GRAPH_H
#define GRAPH_H

#include "vertex.h"
class SimplexCollection;

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <queue>



// This class represents an undirected simplicial graph.
// Simplicial means that there are no self-loops or multiple edges.
class Graph{
    
    std::map<Vertex,std::set<Vertex>> adj;
    
public:
    
    // Constructs a graph with zero vertices
    Graph();
    
    // Constructs a graph with given adjacency lists
    // Does not check whether adj is symmetric, this must be ensured by the caller!
    Graph(const std::map<Vertex,std::set<Vertex>> &adj);
    
    // Reads a graph from an input stream
    // Format: {x_1 y_1 x_2 y_2 ... x_M y_M}
    // where x_i are vertices and y_i are either vertices or sets of vertices,
    // meaning that, for every i, x_i is adjacent to every vertex in y_i.
    // If y_i is a vertex (not a set), then an optional ":" between x_i and y_i is allowed.
    // The delimiters (first and last chars) can be {}, [] or ().
    // As separators between x_i, y_i commas and colons are allowed.
    friend std::istream& operator >> (std::istream &s, Graph &G);
    
    // Writes a graph on an output stream
    friend std::ostream& operator << (std::ostream &s, const Graph &G);
    
    int NumVertices() const;
    int NumEdges() const;
    
    // Does the graph have vertex v?
    bool IsVertex(const Vertex &v) const;
    
    // Adds an edge between u and v (if not present already, and u != v).
    // It is not required that the vertices be already part of the graph:
    // if they are missing, they are added.
    void Connect(const Vertex &v, const Vertex &u);
    
    void AddVertex(const Vertex &v); // Has no effect if v is already a vertex.
    void RemoveVertex(const Vertex &v); // Has no effect if v does not appear in the graph.
    
    // Returns true if u and v are in the graph and are equal or connected by an edge.
    bool Adjacent(const Vertex &v, const Vertex &u) const;
    
    // Returns the set of vertices of G.
    std::set<Vertex> Vertices() const;
    
    // Returns the set of vertices adjacent to v.
    // v must be a vertex of the graph.
    const std::set<Vertex>& Neighbors(Vertex v) const;
    
    // Edges; only one among {u,v} and {v,u} is returned, for u and v distinct adjacent vertices.
    std::vector<std::pair<Vertex,Vertex>> Edges() const;
    
    // Returns the set of vertices x such that x is adjacent (and distinct) to every vertex in V.
    std::set<Vertex> CommonNeighbors(const std::set<Vertex> &V) const;
    
    // Every vertex adjacent to y becomes adjacent to x, and y is deleted (if != x).
    // x and y must be vertices of the graph.
    void IdentifyVertices(const Vertex &x, const Vertex &y);
    
    // Adds a new vertex subdividing the edge [x,y]. It becomes adjacent to common neighs. of x and y.
    // x and y must be distinc adjacent vertices.
    void SubdivideEdge(const Vertex &x, const Vertex &y);
    
    // Determines whether x and y are distinct vertices which are part of a square (an induced subgraph isomorphic to a circle of length 4).
    bool InSquare(const Vertex &x, const Vertex &y) const;
    
    // A circular graph with n vertices, numbered from 0 to n-1 in order
    static Graph Circle(int n);
    
    // Returns the join with another graph.
    // Vertices of the other graph are relabelled as needed.
    Graph Join(const Graph &other) const;
    
    // Takes two copies of G-v, identifying the neighbous of the two copies of v.
    Graph VertexDouble(const Vertex &v) const;
    
    // The induced subgraph whose vertices are in V.
    Graph InducedSubgraph(const std::set<Vertex> &V) const;
    
    // Returns a maximal simple path, meaning that it cannot be extended without repeating vertices.
    std::vector<Vertex> MaximalSimplePath() const;
    
private:
    
    // Input:
    // - s: a set of vertices of G forming a clique
    // - cone_v: a set of vertices of G, each of them adjacent to every vertex of s, but not in s
    // Output: the collection of cliques in G containing s and contained in the union of s and cone_v 
    SimplexCollection CliqueCollection(std::set<Vertex> &s, const std::set<Vertex> &cone_v) const;
    
    // Input:
    // - s: a set of vertices of G forming a clique.
    // Output: the collection of cliques in G containing s.
    SimplexCollection CliqueCollection(std::set<Vertex> &s) const;
    
public:
    
    // Returns the collection of cliques of G.
    SimplexCollection CliqueCollection() const;
    
    // Returns a graph whose vertices are relabelled as prescribed by f.
    // f is not required to preserve the ordering.
    Graph Relabelled(const std::map<Vertex,Vertex> &f) const;

    // Returns a graph whose vertices are relabelled starting from 0, preserving the vertex ordering.
    Graph Relabelled() const;
    
    // Computes gamma_2 = 16 - 8V + 4E - 2F + T, where:
    // - V = number of vertices;
    // - E = number of edges;
    // - F = number of triangles (cliques of size 3);
    // - T = number of tetrahedrons (cliques of size 4);
    int Gamma2() const;
    
    // Returns true if f sends adjacent vertices of G to adjacent (or equal) vertices of H.
    static bool IsHomomorphism(const Graph &G, const Graph &H, const std::map<Vertex,Vertex> f);
};

#endif
