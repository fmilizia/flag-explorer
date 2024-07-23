#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "vertex.h"

#include <vector>
#include <set>
#include <map>
#include <string>



// A Simplex encodes a set of vertices.
// This class exposes methods to perform usual operations on simplices.
class Simplex{
    
    std::vector<Vertex> vertex_vector; // A vector containing the vertices of the simplex, in increasing order.
    unsigned int rank; // Rank of the simplex, i.e., the cardinality of its vertex set.
    
public:
    
    Simplex();
    Simplex(const std::vector<Vertex> &vertex_vector);
    Simplex(const std::set<Vertex> &vertex_set);
    
    bool operator <(const Simplex &other) const;
    bool operator ==(const Simplex &other) const;
    
    unsigned int Rank() const;
    
    const std::vector<Vertex>& GetVertices() const;
    std::vector<Vertex> Vertices() const;
    
    // Returns a string describing the simplex.
    // Format: [v_1 v_2 ... ]
    std::string ToString() const;
    
    // Returns the codimension-1 face opposite to the given vertex.
    // If v is not a vertex, returns a simplex identical to the current one.
    Simplex Facet(const Vertex &v) const;
    
    // Returns a vector containing all the codimension-1 faces of the simplex.
    std::vector<Simplex> Facets() const;
    
    // Returns the subsimplex whose vertices are those NOT contained in s.
    Simplex OppositeFace(const Simplex &s) const;
    
    // Returns a vector storing all the faces of the simplex, of any dimension.
    std::vector<Simplex> Faces() const;
    
    // Returns a new simplex with vertices relabelled as prescribed by f.
    // f is not required to preserve order.
    Simplex Relabel(const std::map<Vertex,Vertex> &f) const;
    
    // Returns a simplex by merging the current one with another.
    Simplex Merge(const Simplex &other) const;
};

#endif
