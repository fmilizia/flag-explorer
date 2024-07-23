#ifndef SIMPLEX_COLLECTION_H
#define SIMPLEX_COLLECTION_H

#include "vertex.h"
class Simplex;
class Graph;

#include <vector>
#include <set>
#include <map>
#include <string>



// This class encodes a set of Simplex objects.
// The set does not have to be a simplicial complex (i.e., subset-closed).
// It exposes methods to merge, relabel, take join of collections, among others.
class SimplexCollection{
    
    std::set<Simplex> simplex_set;
    
public:
    
    typedef std::set<Simplex>::iterator ptr;
    
    SimplexCollection();
    SimplexCollection(const std::set<Simplex> &simplex_set);
    
    const std::set<Simplex>& GetSimplices() const;
    std::set<Simplex> Simplices() const;
    
    // Returns a string describing the collection.
    // Format: {S_1 S_2 ...}
    //         where the S_i's are the simplices in the collection.
    std::string ToString() const;
    
    // Returns the set of vertices appearing in the collection.
    std::set<Vertex> Vertices() const;
    
    int Rank() const;
    
    // Returns true if the collection of simplices is subset-closed, false otherwise.
    bool IsSimplicialComplex() const;
    
    // Adds a simplex and returns:
    // - A pointer to the newly inserted simplex (or an equal one already there);
    // - Whether a new simplex was actually inserted.
    std::pair<ptr,bool> Insert(const Simplex &s);
    
    // Inserts all the simplices of the provided collection.
    void Merge(const SimplexCollection &other);
    
    // Returns a collection whose vertices are relabelled as prescribed by f.
    // f is not required to preserve the ordering.
    SimplexCollection Relabelled(const std::map<Vertex,Vertex> &f) const;
    
    // Returns a collection whose vertices are relabelled starting from 0, preserving the vertex ordering.
    SimplexCollection Relabelled() const;
    
    // Returns the join with another collection C.
    // Vertices of C are relabelled if needed.
    SimplexCollection Join(const SimplexCollection &C) const;
    
    // Returns a collection in which the vertices U and V have been identified.
    // U is kept and V is discarded.
    SimplexCollection Identify(const Vertex &U, const Vertex &V);
    
    // Triangulation of the circle, with n vertices, numbered from 0 to n-1 in order.
    static SimplexCollection Circle(int n);

    // Returns the one-skeleton, as a graph.
    Graph OneSkeleton() const;
};

#endif
