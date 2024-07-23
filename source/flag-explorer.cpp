#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <bitset>
#include <cassert>
#include <random>

#include "vertex.cpp"
#include "simplicial_complex.cpp"
#include "graph.cpp"
#include "localpicture.cpp"
#include "localpicture_homfinder.cpp"



const std::string localpictures_folder = "localpictures/";
const std::array<std::string, 3> localpictures_files{
    "univ_5_1_T10.txt",
    "univ_5_1_T12_B.txt",
    "univ_5_1_T12_A.txt",
};
constexpr int n_localpicture = localpictures_files.size();



// Device for pseudo-random number generation
std::mt19937 random_box;
void PrepareRandomBox(int seed){
    std::cerr << "SEED: " << seed << std::endl;
    random_box = std::mt19937(seed);
}
void PrepareRandomBox(){
    std::random_device rd;
    int seed = rd();
    PrepareRandomBox(seed);
}



// Read a graph from a file.
// If the same file has been already used to read a graph, the file is not 
// reopened, because already read graphs are kept in memory.
Graph GraphFromFile(std::string filename){
    static std::map<std::string, Graph> memo;
    if(memo.count(filename)) return memo.at(filename);
    
    std::ifstream fin(filename);
    if(fin.fail()) throw std::runtime_error("File " + filename + " has some problems");
    
    Graph G;
    DiscardComments(fin);
    fin >> G;
    
    memo.emplace(filename, G);
    return memo.at(filename);
}



// A struct containing a graph, two vertices forming an edge in it, and the 
// LocalPicture constructed from them.
struct GraphLocalPicture{
    const Graph graph;
    const LocalPicture locpic;
    const Vertex u, v;
    
    std::map<std::string, std::map<Vertex,Vertex>> reductions;
    
    GraphLocalPicture(const Graph &G, const Vertex &u, const Vertex &v):
        graph(G),
        locpic(G, u, v),
        u(u),
        v(v),
        reductions()
    {}
};



// Extracts a GraphLocalPicture from a file containing a graph and two vertices
// forming an edge.
GraphLocalPicture LocalPictureFromFile(std::string s){
    static std::map<std::string, GraphLocalPicture> memo;
    if(memo.count(s)) return memo.at(s);
    
    std::ifstream fin(s);
    if(fin.fail()) throw std::runtime_error("File " + s + " has some problems");
    
    Graph G;
    DiscardComments(fin);
    fin >> G;
    
    Vertex x, y;
    DiscardComments(fin);
    fin >> x >> y;
    assert(x != y);
    
    GraphLocalPicture PG(G, x, y);
    DiscardComments(fin);
    fin >> PG.reductions;
    
    memo.emplace(s, PG);
    return memo.at(s);
}



// Changes the graph G by collapsing edges not contained in squares,
// until every edge is contained in a square.
// Selects randomly the edges to collapse.
void CollapseToSquares(Graph &G){
    auto e = G.Edges();
    std::vector<std::pair<Vertex,Vertex>> es;
    
    for(auto p: e) if(not G.InSquare(p.first, p.second)) es.emplace_back(p.first, p.second);
    
    if(es.size() == 0) return;
    
    std::uniform_int_distribution<> distrib(0, es.size()-1);
    int u = distrib(random_box);
    
    G.IdentifyVertices(es[u].first, es[u].second);
    CollapseToSquares(G);
}



// Changes a graph G by performing, in order:
// - sub edge-subdivisions;
// - col edge-collapses (of edges not contained in squares).
// The edges on which to apply these operations are chosen randomly.
// If bl == true (stands for BigLinks) the edge to collapse is (at every step)
// chosen among those having the biggest link size.
// If at some point all edges are contained in squares, the procedure ends
// (without errors) even if the number of performed collapses is less than col.
void NextRandomTriangulation(Graph &G, int sub, int col, bool bl){
    for(int j=0; j<sub; j++){ //Perform edge subdivisions
        auto e = G.Edges();
        std::uniform_int_distribution<> distrib(0, e.size()-1);
        int u = distrib(random_box);
        G.SubdivideEdge(e[u].first, e[u].second);
    }
    for(int j=0; j<col; j++){ //Perform edge collapses
        auto e = G.Edges();
        int cn = 0;
        std::vector<std::pair<Vertex,Vertex>> es;
        for(auto p: e) if(not G.InSquare(p.first, p.second)){
            int ncn = G.CommonNeighbors({p.first,p.second}).size();
            if(bl and ncn > cn){
                cn = ncn;
                es.clear();
            }
            if(ncn >= cn or not bl) es.emplace_back(p.first, p.second);
        }
        if(es.empty()) break;
        std::uniform_int_distribution<> distrib(0, es.size()-1);
        int u = distrib(random_box);
        G.IdentifyVertices(es[u].first, es[u].second);
    }
}



// Checks if there is a map between local pictures coming from two graphs.
// The target local picture (B) is assumed to be defined around an almost-omniscent edge.
// If there is a map, a map between the GRAPHS is returned.
// Otherwise, an empty map is returned.
std::map<Vertex,Vertex> SpecialHomomorphism(const GraphLocalPicture &A, const GraphLocalPicture &B){
    for(int k=0; k<2; k++){ //invert A?
        for(int j=0; j<2; j++){ //reflect A?
            LocalPicture::HomFinder F(A.locpic, B.locpic, j > 0, k > 0);
            if(F.HomExists()){
                auto FF = F.Homomorphism();
                if(j){
                    FF[A.u] = B.v;
                    FF[A.v] = B.u;
                }else{
                    FF[A.u] = B.u;
                    FF[A.v] = B.v;
                }
                
                Vertex out = -1;
                for(const Vertex &x: B.graph.Vertices()) if(x != B.u and x != B.v and not B.locpic.IsVertex(x)){
                    assert(out == -1);
                    out = x;
                }
                for(const Vertex &x: A.graph.Vertices()) if(x != A.u and x != A.v and not A.locpic.IsVertex(x)){
                    assert(out >= 0);
                    FF[x] = out;
                }
                
                for(auto &p: FF) assert(A.graph.IsVertex(p.first));
                assert(Graph::IsHomomorphism(A.graph, B.graph, FF));
                return FF;
            }
        }
    }
    return {};
}



// The following declarations and methods allow to store/read a history
// of analyzed triangulations.
// The history is a vector of objects of type ProcessedTriangulations.
// Each of these objects is meant to represent the event of finding a
// triangulation (which might be isomorphic to a previously encountered one) and
// contains:
// - A number 'first_isomorphic', which is the index (in the history vector) of 
//   the first encountered triangulation which is isomorphic to the current one;
// - A number 'isomorphic_count' which tells how many triangulations of the past
//   are isomorphic to the current one;
// - The encoding of the triangulation, if isomorphic_count == 0 (which means
//   that the current triangulation had never appeared before);
// - A bitset 'dominance' telling to which of the "target localpictures" a
//   homomorphism from the current triangulation has been found (again, only if
//   isomorphic_count == 0).
struct ProcessedTriangulation{
    SimplicialComplex::Encoding E;
    int first_isomorphic;
    int isomorphic_count;
    std::bitset<n_localpicture> dominance;
    
    ProcessedTriangulation(const SimplicialComplex::Encoding &E, int first_iso, const std::bitset<n_localpicture> &dom):
        E(E),
        first_isomorphic(first_iso),
        isomorphic_count(0),
        dominance(dom) {}
    ProcessedTriangulation(int first_iso, int iso_c):
        E({}),
        first_isomorphic(first_iso),
        isomorphic_count(iso_c) {}
};



std::ostream& operator << (std::ostream &o, const ProcessedTriangulation &P){
    o << P.first_isomorphic << ' ' << P.isomorphic_count;
    if(P.isomorphic_count == 0) o << " " << P.dominance << " " << P.E;
    return o;
}



std::vector<ProcessedTriangulation> history;
std::map<SimplicialComplex::Encoding, int> triangulation_collection;



void ReadHistory(std::istream &in){
    std::map<int,int> f_iso_update;
    DiscardComments(in);
    int f_iso, iso_c;
    while(in >> f_iso >> iso_c){
        if(iso_c == 0){
            std::bitset<n_localpicture> dom;
            SimplicialComplex::Encoding E({});
            in >> dom >> E;
            
            const auto [it_e, novel] = triangulation_collection.emplace(E, history.size());
            if(novel){
                history.emplace_back(E, history.size(), dom);
            }else{
                history.emplace_back(history[it_e->second].first_isomorphic, history[it_e->second].isomorphic_count + 1);
                it_e->second = history.size() - 1;
            }
            
            f_iso_update[f_iso] = history.back().first_isomorphic;
        }else{
            auto it_e = triangulation_collection.find(history[f_iso_update[f_iso]].E);
            history.emplace_back(history[it_e->second].first_isomorphic, history[it_e->second].isomorphic_count + 1);
            it_e->second = history.size() - 1;
        }
    }
}



void ReadHistoryFromFile(std::string s){
    std::ifstream fin(s);
    if(fin.fail()) throw std::runtime_error("File " + s + " has some problems");
    ReadHistory(fin);
}



// Check if the clique complex of G dominates one of the target
// triangulations via a homomorphism of local pictures.
// Assumption: G is the one-skeleton of a flag triangulation of S^3.
std::bitset<n_localpicture> Check(const Graph &G){
    const int gamma2 = 16 - 5 * G.NumVertices() + G.NumEdges();
    std::bitset<n_localpicture> result{};
    if(gamma2 == 0) assert(G.NumVertices() == 8);
    if(gamma2 > 0){
        for(int i=0; i<n_localpicture; i++){
            const GraphLocalPicture TP = LocalPictureFromFile(localpictures_folder + localpictures_files[i]);
            
            bool dominates_s = false;
            for(auto e: G.Edges()){
                const GraphLocalPicture PG(G, e.first, e.second);
                if(not SpecialHomomorphism(PG, TP).empty()){
                    dominates_s = true;
                    break;
                }
            }
            if(dominates_s){
                result[i] = true;
            }
        }
    }
    return result;
}



// Process the triangulation given by the clique complex of G.
// Assumption: G is the one-skeleton of a flag triangulation of S^3.
// This adds an entry in the history vector and prints it on the stdout stream.
void ProcessTriangulation(const Graph &G){
    SimplicialComplex::Encoding E = SimplicialComplex(G.CliqueCollection()).StandardEncoding();
    const auto [it_e, novel] = triangulation_collection.emplace(E, history.size());
    if(novel){
        auto res = Check(G);
        history.emplace_back(E, history.size(), res);
    }else{
        history.emplace_back(history[it_e->second].first_isomorphic, history[it_e->second].isomorphic_count + 1);
        it_e->second = history.size() - 1;
    }
    std::cout << history.back() << std::endl;
}



// This is the procedure that implements a "random walk" on the space of 
// flag triangulations of S^3.
void Explore(Graph G, int iterations){
    int sub;
    int col;
    bool bl = false;
    int num_v = G.NumVertices();
    
    for(int i=0; i<iterations; i++){
        std::clog << i << ")\t";
        
        sub = 30; // How many subdivisions in one iteration of the process.
        
        col = sub-3; // How many collapses in one iteration of the process.
        if(num_v > 12) col += 3;
        if(num_v > 30) col += 2;
        col++;
        
        if(i > 0) NextRandomTriangulation(G, sub, col, bl);
        
        Graph H = G;
        CollapseToSquares(H);
        int gamma2 = 16 - 5 * H.NumVertices() + H.NumEdges();
        num_v = H.NumVertices();
        
        if(num_v > 25) bl = true;
        if(num_v <= 15) bl = false;
        
        std::clog << G.NumVertices() << "\tV = " << H.NumVertices() << "\tg_2 = " << gamma2 << std::endl;
        
        ProcessTriangulation(H);
    }
    
    int count_tot = triangulation_collection.size();
    std::clog << "Checked " << count_tot << " pairwise non-isomorphic triangulations" << std::endl;
}



int main(){
    PrepareRandomBox(7978679); // Put any seed for pseudo-random number generation
    
    // Possibly, read a history of previously found triangulations
    // ReadHistoryFromFile("flag_history.txt");
    
    Graph G = Graph::Circle(4).Join(Graph::Circle(4));
    Explore(G, 1000);
    
    return 0;
}
