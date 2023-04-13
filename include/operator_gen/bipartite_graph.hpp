#ifndef UTILS_BIPARTITE_GRAPH_HPP
#define UTILS_BIPARTITE_GRAPH_HPP

#include <vector>
#include <list>
#include <stack>
#include <iostream>
#include <queue>

namespace utils
{
template <typename T, typename V> class bipartite_graph;
template <typename T, typename V> std::ostream& operator<<(std::ostream& os, const utils::bipartite_graph<T, V>& op);




//a structure for storing a bipartite graph that allows for information to be stored at each edge
template <typename T, typename Val>
class bipartite_graph
{
public:
    static void generate_connected_subgraphs(const bipartite_graph<T, Val>& G, std::list<bipartite_graph<T, Val>>& SG)
    {
        size_t nodes_filled = 0;
        std::vector<bool> Uc(G.N());     std::fill(Uc.begin(), Uc.end(), false);
        std::vector<bool> Vc(G.M());      std::fill(Vc.begin(), Vc.end(), false);
        while(nodes_filled < G.N() + G.M())
        {
            size_t Ui = 0;
            bool Uifound = false;
            //find the first U node that has not yet been coloured
            for(size_t i = 0; i < Uc.size() && !Uifound; ++i)
            {
                if(!Uc[i]){Ui = i;  Uifound = true;}
            }

            //if we have found a Ui index, then we proceed to flood fill to identify all of its connected nodes
            if(Uifound)
            {
                std::stack<size_t> Us, Vs;
                bipartite_graph<T, Val> bcg;    
                //insert the first value into Us 
                Us.push(Ui);   
                Uc[Ui] = true;
                bcg.add_U(G.U(Ui));

                //while there is a node that is still in either of our stacks we go through and do as is needed to empty everything
                while(!Us.empty() || !Vs.empty())
                {
                    //we take the top U and add all of its Vs that have not already been added and mark them as coloured
                    if(!Us.empty())
                    {
                        size_t ui = Us.top();  Us.pop();
                        ++nodes_filled;

                        const std::list<size_t>& connected_nodes = G.U_edges(ui);
                        const std::list<Val>& connected_nodes_d = G.U_edges_data(ui);
                        for(auto z : zip(connected_nodes, connected_nodes_d))
                        {
                            size_t vi = std::get<0>(z);     
                            if(!Vc[vi])
                            {
                                bcg.add_V(G.V(vi));
                                Vs.push(vi);
                                Vc[vi] = true;
                            }

                            //connect this u to all of the v's it is connected to
                            bcg.add_edge_v(G.U(ui), G.V(vi), std::get<1>(z));
                        }
                    }

                    if(!Vs.empty())
                    {
                        size_t vi = Vs.top();  Vs.pop();
                        ++nodes_filled;
                        const std::list<size_t>& connected_nodes = G.V_edges(vi);
                        const std::list<Val>& connected_nodes_d = G.V_edges_data(vi);
                        for(auto z : zip(connected_nodes, connected_nodes_d))
                        {
                            size_t ui = std::get<0>(z);
                            if(!Uc[ui])
                            {
                                bcg.add_U(G.U(ui));
                                Us.push(ui);
                                Uc[ui] = true;
                            }

                            //as we are adding all of the connections from u to v we don't need to worry about adding them this side, as the U code will iterate over every allowed u and as a consequence every allowed edge
                            //bcg.add_edge_v(G.U(ui), G.V(vi), std::get<1>(z));
                        }
                    }
                }

                SG.push_back(bcg);
            }
            //if we have already coloured all Uc nodes then search for any uncoloured Vi nodes and add them to the SG list
            else
            {
                bool Vifound = false;
                for(size_t i = 0; i < Vc.size(); ++i)
                {
                    if(!Vc[i])
                    {
                        Vifound = true;
                        bipartite_graph<T, Val> bcg;    bcg.add_V(G.V(i));
                        SG.push_back(bcg);
                        ++nodes_filled;
                    }
                }

                ASSERT(Vifound, "Something went horribly wrong.  We haven't coloured all nodes, but no nodes were found to have not been coloured.");
            }
        }
    }


public:
    bipartite_graph() : _nedge(0) {}
    bipartite_graph(size_t n, size_t m) : _nedge(0)
    {
        CALL_AND_HANDLE(resize(n, m), "Failed to compute bipartite graph");
    }

    bipartite_graph(const bipartite_graph& o) = default;
    bipartite_graph(bipartite_graph&& o) = default;

    bipartite_graph& operator=(const bipartite_graph& o) = default;
    bipartite_graph& operator=(bipartite_graph&& o) = default;
    
    void clear()
    {
        _U.clear();
        _V.clear();
        _Eu.clear();
        _Ev.clear();
        _Eud.clear();
        _Evd.clear();
        _nedge = 0;
    }

    void resize(size_t n, size_t m)
    {
        _U.reserve(n);
        _V.reserve(m);
        _Eu.reserve(n);
        _Ev.reserve(m);
        _Eud.reserve(n);
        _Evd.reserve(m);
    }

    size_t add_U(const T& u)
    {
        size_t ind = _U.size();
        _U.push_back(u);
        _Eu.push_back(std::list<size_t>());
        _Eud.push_back(std::list<Val>());
        return ind;
    }

    size_t add_V(const T& v)
    {
        size_t ind = _V.size();
        _V.push_back(v);
        _Ev.push_back(std::list<size_t>());
        _Evd.push_back(std::list<Val>());
        return ind;
    }

    void add_edge(size_t u, size_t v)
    {
        ASSERT(u < _U.size() && v < _V.size(), "Failed to add edge.  Node indices out of bounds.");
        _Eu[u].push_back(v);
        _Ev[v].push_back(u);

        _Eud[u].push_back(Val());
        _Evd[v].push_back(Val());

        ++_nedge;
    }

    void add_edge(size_t u, size_t v, const Val& d)
    {
        ASSERT(u < _U.size() && v < _V.size(), "Failed to add edge.  Node indices out of bounds.");
        _Eu[u].push_back(v);
        _Ev[v].push_back(u);

        _Eud[u].push_back(d);
        _Evd[v].push_back(d);

        ++_nedge;
    }

    void add_edge_v(const T& u, const T& v)
    {
        size_t ui, vi;
        CALL_AND_HANDLE(ui = U_ind(u), "Cannot add edge.  Nodes not found.");
        CALL_AND_HANDLE(vi = V_ind(v), "Cannot add edge.  Nodes not found.");
        add_edge(ui, vi);
    }

    void add_edge_v(const T& u, const T& v, const Val& d)
    {
        size_t ui, vi;
        CALL_AND_HANDLE(ui = U_ind(u), "Cannot add edge.  Nodes not found.");
        CALL_AND_HANDLE(vi = V_ind(v), "Cannot add edge.  Nodes not found.");
        add_edge(ui, vi, d);
    }

    bool U_contains(const T& u) const
    {
        return std::find(_U.begin(), _U.end(), u) != _U.end();
    }

    bool V_contains(const T& v) const
    {
        return std::find(_V.begin(), _V.end(), v) != _V.end();
    }


    size_t U_ind(const T& u)
    {
        auto ui = std::find(_U.begin(), _U.end(), u);
        ASSERT(ui != _U.end(), "Failed to get index of U.  It is not in the set.");
        return (ui - _U.begin());
    }

    size_t V_ind(const T& v)
    {
        auto vi = std::find(_V.begin(), _V.end(), v);
        ASSERT(vi != _V.end(), "Failed to get index of V.  It is not in the set.");
        return (vi - _V.begin());
    }


    const T& U(size_t i) const
    {
        ASSERT(i<_U.size(), "Failed to access U element.");
        return _U[i];
    }

    T& U(size_t i)
    {
        ASSERT(i<_U.size(), "Failed to access U element.");
        return _U[i];
    }

    const T& V(size_t i) const
    {
        ASSERT(i<_V.size(), "Failed to access V element.");
        return _V[i];
    }

    T& V(size_t i)
    {
        ASSERT(i<_V.size(), "Failed to access V element.");
        return _V[i];
    }

    const std::list<size_t>& U_edges(size_t i) const
    {
        ASSERT(i<_Eu.size(), "Failed to access U element.");
        return _Eu[i];
    }

    const std::list<size_t>& V_edges(size_t i) const
    {
        ASSERT(i<_Ev.size(), "Failed to access U element.");
        return _Ev[i];
    }

    const std::list<Val>& U_edges_data(size_t i) const
    {
        ASSERT(i<_Eud.size(), "Failed to access U element.");
        return _Eud[i];
    }

    const std::list<Val>& V_edges_data(size_t i) const
    {
        ASSERT(i<_Evd.size(), "Failed to access U element.");
        return _Evd[i];
    }

    size_t N() const{return _U.size();}
    size_t M() const{return _V.size();}
    size_t nedge() const{return _nedge;}

    friend std::ostream& operator<< <T, Val>(std::ostream& os, const bipartite_graph<T, Val>& op);
protected:
    std::vector<T> _U;
    std::vector<T> _V;
    std::vector<std::list<size_t> > _Eu;
    std::vector<std::list<size_t> > _Ev;

    std::vector<std::list<Val> > _Eud;
    std::vector<std::list<Val> > _Evd;
    size_t _nedge;

};



template <typename T, typename V> 
std::ostream& operator<<(std::ostream& os, const utils::bipartite_graph<T, V>& op)
{
    os << "{\"U\" : [" << std::endl;
    std::string sep("");
    for(size_t i = 0; i < op.N(); ++i)
    {
        os << sep << "\"" << op.U(i) << "\"";
        sep = std::string(", ");
    }
    os << std::endl << "]," << std::endl;

    os << "\"V\" : [" << std::endl;
    sep = std::string("");
    for(size_t i = 0; i < op.M(); ++i)
    {
        os << sep << "\"" << op.V(i) << "\"";
        sep = std::string(", ");
    }
    os << std::endl << "]," << std::endl;

    os << "\"E\" : [" << std::endl;
    sep = std::string("");
    for(size_t i = 0; i < op.N(); ++i)
    {
        const auto& Ue = op.U_edges(i);
        const auto& Ued = op.U_edges_data(i);
        for(auto z : zip(Ue, Ued))
        {
            os << sep << "[ \"" << op.U(i) << "\", \"" << op.V(std::get<0>(z)) << "\", " << std::get<1>(z) << "] ";
            sep = std::string(", ");
        }
    }
    os << std::endl << "]}" << std::endl;
    return os;
}


class bipartite_matching
{
public:
    bipartite_matching() : _nU(0), _nV(0) {}
    template <typename T, typename V>
    bipartite_matching(const bipartite_graph<T, V>& o)
    {
        compute_maximal_matching(o);
    }

    /*
     * This seems to work but it will be necessary to perform more testing before I am happy to use it in a production environment
     */
    template <typename T, typename V>
    size_t compute_maximal_matching(const bipartite_graph<T, V>& o)
    {
        _nU = o.N();     _nV=o.M();
        match_U.resize(_nU);     std::fill(match_U.begin(), match_U.end(), -1);
        match_V.resize(_nV);   std::fill(match_V.begin(), match_V.end(), -1);
        dist.resize(_nU);
        

        bool may_contain_augmenting_path = true;
        //while there is a chance that there may be augmenting paths.  
        size_t nterms = 0;
        while(may_contain_augmenting_path)
        {
            //perform the breadth first search and separate the graph into layers of matching and not matching edges
            bfs(o);

            int n_augmenting_paths = 0;

            //now iterate over all edges in set U and if they are unpaired see if they are the root of an augmenting path
            for(size_t u = 0; u < _nU; ++u)
            {
                if( (match_U[u] == -1) && dfs(o, u)) 
                {
                    ++n_augmenting_paths;
                    ++nterms;
                }
                
                //if we didn't find any augmenting paths then we are done
                if(n_augmenting_paths == 0)
                {
                    may_contain_augmenting_path = false;
                }
            }
        }
        return nterms;
    }

    std::vector<std::pair<size_t, size_t>> edges() const
    {
        std::vector<std::pair<size_t, size_t>> edg;
        for(size_t u = 0; u < _nU; ++u)
        {
            if(match_U[u] != -1)
            {
                edg.emplace_back(u, static_cast<size_t>(match_U[u]));
            }
        }
        return edg;
    }


    template <typename T, typename Val>
    std::pair<std::vector<size_t>, std::vector<size_t>> minimum_vertex_cover(const bipartite_graph<T, Val>& o) const
    {
        return _minimum_vertex_cover_v(o);
    }


protected:
    /*
     * This appears to work but it will be necessary to perform more testing before I am happy to use it in a production environment
     */
    template <typename T, typename Val>
    std::pair<std::vector<size_t>, std::vector<size_t>> _minimum_vertex_cover_u(const bipartite_graph<T, Val>& o) const
    {
        std::vector<size_t> U, V;
        std::vector<bool> U_visited(_nU);     std::fill(U_visited.begin(), U_visited.end(), false);
        std::vector<bool> V_visited(_nV);     std::fill(V_visited.begin(), V_visited.end(), false);

        //iterate over all indices in U and if the element is not included in the matching perform 
        //a dfs from the current node to identify and alternating path
        for(size_t u = 0; u < _nU; ++u)
        {
            dfs_visited_u(o, u, U_visited, V_visited);
        }

        //now that we have done the search we can actually work out the covering.  If the v index was visited then we add it into
        //the set.  Otherwise we check the u value that was matched with it and add it to the vertex cover if it wasn't visited
        for(size_t v = 0; v < _nV; ++v)
        {
            if(match_V[v] != -1)
            {
                int u = match_V[v];
                if(V_visited[v])
                {
                    V.push_back(v);
                }
                else if(!U_visited[u])
                {
                    U.push_back(u);
                }
            }
        }

        //get the set of unmatched vertices in th right set
        return {U, V};
    }

    template <typename T, typename Val>
    std::pair<std::vector<size_t>, std::vector<size_t>> _minimum_vertex_cover_v(const bipartite_graph<T, Val>& o) const
    {
        std::vector<size_t> U, V;
        std::vector<bool> U_visited(_nU);     std::fill(U_visited.begin(), U_visited.end(), false);
        std::vector<bool> V_visited(_nV);     std::fill(V_visited.begin(), V_visited.end(), false);

        //iterate over all indices in U and if the element is not included in the matching perform 
        //a dfs from the current node to identify and alternating path
        for(size_t v = 0; v < _nV; ++v)
        {
            dfs_visited_v(o, v, U_visited, V_visited);
        }

        //now that we have done the search we can actually work out the covering.  If the v index was visited then we add it into
        //the set.  Otherwise we check the u value that was matched with it and add it to the vertex cover if it wasn't visited
        for(size_t u = 0; u < _nU; ++u)
        {
            if(match_U[u] != -1)
            {
                int v = match_U[u];
                if(U_visited[u])
                {
                    U.push_back(u);
                }
                else if(!V_visited[v])
                {
                    V.push_back(v);
                }
            }
        }

        //get the set of unmatched vertices in th right set
        return {U, V};
    }

protected:
    //perform the breadth first search to determine to construct the augmenting paths  
    template <typename T, typename V>
    void bfs(const bipartite_graph<T, V>& o)
    {
        std::queue<int> q;

        //add any free vertices to the queue - this is the first layer of vertices
        for(size_t u = 0;  u < _nU; ++u)
        {
            //if we find a free vertex then we add it to the queue
            if(match_U[u] == -1)
            {
                q.push(u);  dist[u] = 0;
            }
            //otherwise set its distance to -1 so that it will be considered in the next iteration
            else
            {
                dist[u] = -1;
            }
        }

        //now iterate over the queue and determine if it can provide a shorter path
        while(!q.empty())
        {
            //get the first free vertex and remove it from the queue
            int u = q.front();  q.pop();

            //if the distance associated with this node is not -1 (that is it can provide a shorter path)
            if(dist[u]  != -1)
            {
                //iterate over the vertices of v that are connected to this node u
                for(auto v : o.U_edges(u))
                {
                    //add the v index to the augmenting path  
                    if(match_V[v] != -1 && dist[match_V[v]] == -1)
                    {   
                        dist[match_V[v]] = dist[u]+1;

                        //and add the node connected to v back into the queue
                        q.push(match_V[v]);
                    }
                }
            }
        }
    }

    //perform the depth first search to add augmenting paths to the current matching
    template <typename T, typename V>
    bool dfs(const bipartite_graph<T, V>& o, int u)
    {
        //check the base case of the recursion
        for(auto v : o.U_edges(u))
        {
            //if any of the nodes aren't in the matching then we add it and we are done.
            if(match_V[v] == -1)
            {
                match_U[u] = v;     match_V[v] = u;     return true;
            }
        }

        //dist[u] = -1;

        //now the recursive case of the dfs
        for(auto v : o.U_edges(u))
        {
            if(dist[match_V[v]] == dist[u]  + 1 && dfs(o, match_V[v]))
            {
                match_U[u] = v; match_V[v] = u;
                return true;
            }
        }
        //dist[u] = -1;

        //and if we get to the end there is no augmenting path beginning with u
        return false;
    }


    template <typename T, typename V>
    void dfs_visited_u(const bipartite_graph<T, V>& o, int u, std::vector<bool>& U_visited, std::vector<bool>& V_visited) const
    {
        if(match_U[u] == -1 && !U_visited[u])
        {
            U_visited[u] = true;

            for(auto v : o.U_edges(u))
            {
                //this node hasn't been visited
                if(!V_visited[v])
                {
                    //mark it as being visited
                    V_visited[v] = true;

                    if(match_V[v] != -1)
                    {
                        dfs_visited_u(o, match_V[v], U_visited, V_visited);
                    }
                }
            }
        }
    }

    template <typename T, typename V>
    void dfs_visited_v(const bipartite_graph<T, V>& o, int v, std::vector<bool>& U_visited, std::vector<bool>& V_visited) const
    {
        if(match_V[v] == -1 && !V_visited[v])
        {
            V_visited[v] = true;

            for(auto u : o.V_edges(v))
            {
                //this node hasn't been visited
                if(!U_visited[u])
                {
                    //mark it as being visited
                    U_visited[u] = true;

                    if(match_U[u] != -1)
                    {
                        dfs_visited_v(o, match_U[u], U_visited, V_visited);
                    }
                }
            }
        }
    }

    size_t _nU, _nV;
    std::vector<int> match_U, match_V, dist;

};


}

#endif

