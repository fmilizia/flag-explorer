#ifndef SIMPLICIAL_COMPLEX_H
#define SIMPLICIAL_COMPLEX_H

#include "vertex.h"
class Simplex;
class SimplexCollection;
class Graph;

#include <vector>
#include <string>
#include <map>
#include <set>



// This class encodes a simplicial complex.
// It does not allow manipulations after it has been created; this is to allow
// performing some precomputations in order to make computations, e.g.,  of
// links and stars more efficient.
// If you want to manipulate simplicial complexes and then perform some 
// computations on them, you should first use a SimplexCollection object, to 
// perform manipulations; then, build a SimplicialComplex object from it. 
class SimplicialComplex{
    
private:
    
    std::vector<Simplex> simplices; // Stores the simplices of the complex, in increasing order.
    // The id of a simplex is defined as its position in the vector.
    // The empty simplex is always included, its id is 0.
    // The id of a singleton {v} is always equal to v+1.
    
    std::vector<unsigned int> simplex_count; // i -> The number of simplices with rank < i, for i <= rank + 1.
    // As a consequence, the simplices of rank k are those with index from simplex_count[k] included to simplex_count[k+1] excluded.
    
    unsigned int rank; // max{rk(s): s simplex in the complex}, i.e., dimension + 1.
    
    std::vector<Vertex> label; // Stores the original vertices of the complex.
    // When the SimplicialComplex is created, the vertices are relabelled from 0; label[i] is the original label of vertex i.
    
    std::vector<std::vector<unsigned int>> link; // Stores the link of each simplex of the complex, as (ordered) vectors of id's.
    std::vector<std::vector<unsigned int>> open_star; // Stores the open star of each simplex of the complex, as (ordered) vectors of id's.
    std::vector<std::vector<unsigned int>> facet; // facet[s_id][j] = id of the j-th facet of the simplex with id s_id (obtained by removing its j-th vertex).
    
public:
    
    SimplicialComplex(const SimplexCollection &collection);
    
private:
    
    void Precomputations();
    void ComputeFacets();
    void ComputeLinksAndStars();
    
public:
    
    std::string ToString() const;
    
    SimplexCollection Collection() const;
    
    unsigned int Rank() const;
    
    int Dimension() const;
    
    int NumVertices() const;
    
    // Returns the SimplexId corresponding to the given simplex s, if the latter is part of the complex;
    // otherwise, returns an invalid SimplexId.
    unsigned int GetSimplexId(const Simplex &s) const;
    
    bool ValidSimplexId(const unsigned int &id) const;
    
    // Returns the link of a given simplex.
    const std::vector<unsigned int>& GetSimplexLink(const unsigned int &id) const;
    
    // Given a simplex and one of its facets (the i-th one), gives the simplex sharing this facet, assuming it is unique.
    unsigned int OppositeFacet(const unsigned int &id, const unsigned int &i) const;
    
    friend bool IsSimplicial(const std::map<Vertex,Vertex> &f, const SimplicialComplex &S, const SimplicialComplex &T);
    friend int DegreeMod2(const std::map<Vertex,Vertex> &f, const SimplicialComplex &S, const SimplicialComplex &T);
    
    // Given a simplex and an order on its vertices  (as a permutation of {0,...,rank-1}),
    // it assigns a number f[v] to every vertex v, 
    // telling the order of "discovery" during a "visit" of the complex (starting from 0).
    // Assumption: the simplicial complex is a strongly-connected pseudomanifold.
    std::vector<int> VisitOrder(const unsigned int &s_id, const std::vector<int> &s_order) const;
    
    
    // This subclass furnishes a way to "encode" a simplicial complex which is
    // assumed to be a flag strongly-connected pseudomanifold.
    // The encodings of two complexes are equal if and only if the complexes are
    // isomorphic; moreover, the set of encodings is totally ordered. 
    class Encoding{
        std::vector<int> enc;
    
    public:
        Encoding(std::vector<int> enc);
        Encoding();
        bool operator == (const Encoding &other) const;
        bool operator < (const Encoding &other) const;
        friend std::ostream& operator << (std::ostream &o, const Encoding &E);
        friend std::istream& operator >> (std::istream &s, Encoding &E);
        Graph OneSkeleton() const; // The 1-skeleton of a triangulation whith the given encoding.
    };
    
    // Computes the encoding of the simplicial complex.
    // Assumption: the complex is a flag strongly-connected pseudomanifold.
    Encoding StandardEncoding() const;
    
    
    
    // The following subclass can be used to determine if two simplicial complexes X and Y are isomorphic.
    // It assumes that X and Y are strongly-connected pseudomanifolds of the same dimension.
    // Example of how to use this class:
    //   X = SimplicialComplex(...); Y = SimplicialComplex(...);
    //   SimplicialComplex::IsomorphismFinder IF(X,Y);
    //   if(IF.Isomorphic()){ std::vector<Vertex> v = IF.Isomorphism(); ... }
    class IsomorphismFinder{
        const SimplicialComplex &X;
        const SimplicialComplex &Y;
        unsigned int rank;
        bool computed;
        bool isomorphic;
        std::vector<Vertex> f; // Contains the found isomorphism, if it exists and has been computed. {Vertex -> Vertex}
        
    public:
        
        IsomorphismFinder(const SimplicialComplex &X, const SimplicialComplex &Y);
        
    private:
        
        void BasicChecks();
        bool AttemptFromFacets(const unsigned int &sX_id, const unsigned int &sY_id, const std::vector<Vertex> &sY_order);
        void Compute();
        
    public:
        
        bool Isomorphic();
        std::map<Vertex,Vertex> Isomorphism();
    };
    
    // The following subclass can be used to determine if a simplicial complexes X dominates another complex Y,
    // meaning that there is a simplicial map from X to Y of nonzero degree.
    // It assumes that X and Y are orientable strongly-connected pseudomanifolds of the same dimension.
    // Example of how to use this class:
    //  X = SimplicialComplex(...); Y = SimplicialComplex(...);
    //  SimplicialComplex::NZMapFinder F(X,Y);
    //  if(F.Exists()){ std::vector<Vertex> v = F.NZMap(); ... }
    class NZMapFinder{
        const SimplicialComplex &X;
        const SimplicialComplex &Y;
        unsigned int rank;
        bool computed;
        bool found;
        std::vector<Vertex> f; // Contains the found map, if it exists and has been computed. {Vertex -> Vertex}
        
    public:
        
        NZMapFinder(const SimplicialComplex &X, const SimplicialComplex &Y);
        
    private:
        
        void BasicChecks();
        void Compute();
        void Compute(std::vector<std::set<Vertex>> candidates);
        int DegreeMod2(); // Should be replaced by the degree over Z, but modulo 2 it is easier to implement.
        
    public:
        
        bool Exists();
        std::map<Vertex,Vertex> NZMap();
    };
};

#endif
