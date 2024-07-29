#ifndef LOCALPICTURE_H
#define LOCALPICTURE_H

#include "vertex.h"
class Graph;

#include <vector>
#include <set>
#include <queue>
#include <map>


// The class LocalPicture encodes a structure consisting of:
// - A graph;
// - A partition of the vertices of the graph into three subsets: equator and southern/nothern hemispheres;
// - A subset of vertices which are considered "marked",
// such that the subgraph induced by the equatorial vertices is a cycle; moreover, the structure includes:
// - An orientation of the equator.
class LocalPicture{
    
    // Internally, vertices are numbered with integers, starting from 0.
    std::vector<std::vector<int>> graph; // v -> the neighbours of v, in increasing order.
    int equator_length = 0; // Vertices from 0 to equator_length-1 are part of the "equator".
    int north_size = 0; // From equator_length to equator_length+north_size-1 they are part of the "northern hemisphere".
    int south_size = 0; // From equator_length+north_size to equator_length+north_size+south_size-1 they are part of the "southern hemisphere".
    std::vector<bool> marked;
    
    bool reflected = false; // Exchange north / south.
    bool inverted = false; // Change the orientation of the equator.
    
    std::vector<Vertex> labels; // Original labels of the vertices.
    std::map<Vertex,int> idx; // Original vertex -> its internal number.
    
    void Init(const Graph &G,
              const std::vector<Vertex> &equator,
              const std::set<Vertex> &north,
              const std::set<Vertex> &south,
              const std::set<Vertex> &marked_vertices);
    
    // Do x and y correspond to adjacent (or equal) vertices?
    bool AdjacentIdx(int x, int y) const;
    
    // The i-th index in the north/south.
    int IthNorthIdx(int i) const;
    int IthSouthIdx(int i) const;
    
    // Is i a north/south/equator index? (The equator is neither north nor south)
    bool IsNorthIdx(int i) const;
    bool IsSouthIdx(int i) const;
    bool IsEquatorIdx(int i) const;
    
    int EquatorSuccessorIdx(int i) const; // The index that comes after i in the equator, following the orientation.
    int EquatorPredecessorIdx(int i) const;
    
public:
    
    LocalPicture(const Graph &G,
            const std::vector<Vertex> &equator,
            const std::set<Vertex> &north,
            const std::set<Vertex> &south,
            const std::set<Vertex> &marked_vertices);
    
    // The following constructor takes a graph G and two vertices x,y.
    // Assumptions:
    // - x and y are distinct and adjacent in G;
    // - The common neighbours of x and y induce a closed path in G.
    // The resulting LocalPicture has:
    // - The common neighbours of x and y as equator;
    // - The neighbours of x at distance >= 2 from y, as northern hemisphere;
    // - The neighbours of y at distance >= 2 from x, as southern hemisphere;
    // - Graph: the one induced in G by the set of vertices of the three types;
    // - Marked vertices are those which are adjacent, in G, to vertices which
    //     are not in the LocalPicture (that is, not adjacent to x nor y).
    LocalPicture(const Graph &G, const Vertex &x, const Vertex &y);
    
    int NumVertices() const;
    int EquatorLength() const;
    
    bool Adjacent(const Vertex &v, const Vertex &w) const;
    int NorthSize() const;
    int SouthSize() const;
    Vertex IthNorth(int i) const;
    Vertex IthSouth(int i) const;
    Vertex IthEquator(int i) const; // Consecutive along the equator, but start and direction do not mean anything.
    Vertex IthVertex(int i) const;
    Vertex EquatorSuccessor(const Vertex &v) const;
    Vertex EquatorPredecessor(const Vertex &v) const;
    bool IsVertex(const Vertex &v) const;
    bool IsNorth(const Vertex &v) const;
    bool IsSouth(const Vertex &v) const;
    bool IsEquator(const Vertex &v) const;
    bool IsMarked(const Vertex &v) const;
    
    void Reflect(); // Exchange north and south.
    void Invert();  // Reverse the equator.
    
    // Check if f describes a homomorphism from P to Q:
    // - Its domain is the set of vertices of P;
    // - It sends adjacent vertices to adjacent or equal vertices;
    // - It sends equatorial vertices to equatorial vertices;
    // - It sends southern vertices to southern or equatorial vertices;
    // - Restricted to equators, it is a nondecreasing map of degree 1 between oriented circles.
    static bool IsHomomorphism(const LocalPicture &P, const LocalPicture &Q, const std::map<Vertex,Vertex> &f);
    
    // Check if f describes a homomorphism from P to Q',
    // where Q' is obtain from Q by reflecting and/or inverting as dictated by the parameters.
    static bool IsHomomorphism(const LocalPicture &P, const LocalPicture &Q, const std::map<Vertex,Vertex> &f, bool reflect, bool invert);
    
    // The group C_2 x C_2 acts on localpictures by reflections / inversions.
    // This function computes a set of representatives for the cosets w.r.t.
    // the subset acting trivially on the locapicture.
    std::vector<std::bitset<2>> ReflectInvertRepresentatives() const;
    
    class HomFinder;
};



// A helper class to find homomorphisms between two LocalPicture objects.
// Example of how to use this class, assuming A and B are two existing LocalPicture objects:
//
// LocalPicture::Homfinder hf(A, B);
// if(hf.HomExists()){
//     std::map<Vertex,Vertex> hom = hf.Homomorphism();
//     -- do something with hom --
// }
class LocalPicture::HomFinder{
    
    const LocalPicture P;
    const LocalPicture Q;
    std::vector<std::vector<int>> rQ; // The complement of Q.graph.
    
    bool inverted;
    bool reflected;
    
    bool computed;
    bool hom_found;
    std::vector<int> f; // The homomorphism. -1 means: not assigned.
    
    std::vector<std::vector<bool>> is_allowed; // is_allowed[i][j] == false means that vert. i (in P) cannot be sent to j (in Q).
    std::vector<int> count_alternatives; // how many possible images for each vertex in P.
    std::priority_queue<std::pair<int,int>> pq; // p-queue of (-#alternatives, vertex in P), to select the one with less alternatives.
    std::vector<int> eqP;
    
    void RemoveAlternative(int v, int u, std::vector<std::pair<int,int>> &deleted_alternatives, std::vector<int> &changed_vertices);
    void ComputeRecursive();
    void Compute();
    
public:
    
    HomFinder(const LocalPicture &P, const LocalPicture &Q);
    HomFinder(const LocalPicture &P, const LocalPicture &Q, bool reflect, bool invert);
    
    bool HomExists();
    std::map<Vertex,Vertex> Homomorphism();
};

#endif
