#ifndef TTNS_LIB_OPERATOR_GEN_AUTO_SUM_OF_PRODUCT_HPP
#define TTNS_LIB_OPERATOR_GEN_AUTO_SUM_OF_PRODUCT_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>

#include <operator_gen/sSOP.hpp>

#include "compressedSOP.hpp"
#include "../utils/bipartite_graph.hpp"

namespace ttns
{


class bipartitionSOP
{
public:
    //a function for taking a bipartitioning of the Hamiltonain.  Here the modes U store the 
    template <typename T>
    static inline bipartite_graph<sPOP, T> form_bipartite_graph(const sSOP<T>& sop, const std::vector<size_t>& modes)
    {
        //sort the modes vector into smodes
        std::vector<size_t> smodes(modes.begin(), modes.end());
        std::sort(smodes.begin(), smodes.end());

        //start by getting the bipartitioned sets of operators
        utils::bipartite_graph<sPOP, T> bpg;
        {
            for(const auto& NBO : sop)
            {
                const auto& pop = NBO.pop();
                sPOP _u, _v;

                //extract tdb
                bipartition_term(pop, smodes, _u, _v);

                size_t _ui, _vi;
                if(!bpg.U_contains(_u))
                {
                    _ui = bpg.add_U(_u);
                }
                else
                {
                    _ui = bpg.U_ind(_u);
                }


                if(!bpg.V_contains(_v))
                {
                    _vi = bpg.add_V(_v);
                }
                else
                {
                    _vi = bpg.V_ind(_v);
                }

                bpg.add_edge(_ui, _vi, NBO.coeff());

            }
        }

        return bpg;
    }


    //take a sum of product form for a Hamiltonian and determine which terms can be combined together to simplify the representation of the Hamiltonian
    template <typename T>
    static inline void simplify_bipartitioning(const sSOP<T>& sop, const std::vector<size_t>& modes)
    {
        //bipartition the graph
        auto bpg = bipartitionSOP::form_bipartite_graph(sop, modes);

        //get a list of all connected subgraphs
        std::list<decltype(bpg)> sbpg;
        decltype(bpg)::generate_connected_subgraphs(bpg, sbpg);

        //now iterate over the list of connected subgraphs and form the minimal vertex cover
        
        std::list<std::list<size_t>> sum_ops;
        size_t counter = 0;
        for(const auto& g : sbpg)
        {
            bipartite_matching bpm(g);
            auto m = bpm.edges();
            //std::vector<size_t> U,V;

            //std::cerr << "Matching: " << std::endl;
            //for(const auto& z : m)
            //{
            //    std::cerr << g.U(std::get<0>(z)) << " <-> " << g.V(std::get<1>(z)) << std::endl;
            //}

            auto [_U, _V] = bpm.minimum_vertex_cover(g);

            //now for all of the nodes in _V we can take all of the connected nodes and construct from these a sum operator 
            for(auto v : _V)
            {
                const auto& ve = g.V_edges(v);
                std::list<size_t> connected_indices;
                for(const auto& _u : ve)
                {
                    if(std::find(_U.begin(), _U.end(), _u) == _U.end())
                    {
                        connected_indices.push_back(_u);
                    }
                }
                sum_ops.push_back(connected_indices);

                std::cout << g.V(v) << std::endl;
                for(auto u : connected_indices)
                {
                    std::cout << "\t: " << g.U(u) << std::endl;
                }
            }

            //increment the index counter based on the number of nodes we will have already had to have treated.
            counter += g.N();
        }

    }

protected:
    static inline void bipartition_term(const sPOP& in, const std::vector<size_t>& ld, sPOP& l, sPOP& r)
    {
        l.clear();  r.clear();
        for(const auto& op : in)
        {
            if(std::binary_search(ld.begin(), ld.end(), op.mode()))
            {
                l*=op;
            }
            else
            {
                r*=op;
            }
        }
    }
};





class AutoSOP
{
protected:
    template <typename T>
    static inline cSSOP<T> compress_sop_simple(const sSOP<T>& sop, size_t nmodes)
    {

        std::vector< std::list<std::pair<sOP, std::list<size_t> > >> sops(nmodes);

        for(const auto& pop : sop)
        {
            
        }
    }

public:
    template <typename T>
    static inline compressedSOP<T> simple(const sSOP<T>& _sop, const system_modes& sys_info)
    {
        try
        {
            sSOP<T> sop;

            CALL_AND_HANDLE(sop = SOPUtils::primitive(_sop, sys_info), "Failed to simplify the SOP object.");

            size_t _nmodes = get_nmodes(sop);
            size_t nmodes = sys_info.nmodes();
            ASSERT(nmodes >= _nmodes, "Failed to construct sum_of_product_operator object with the specified number of modes.  The input object has more modes than have been requested.");

            CALL_AND_HANDLE(return compress_sop_simple(sop, nmodes), "Failed to compute the compressed sop object using the simple strategy.");
        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to compute the simple form of the SOP operator.");
        }
    }

    template <typename T>
    static inline compressedSOP<T> simple(const sSOP<T>& _sop, size_t nmodes = 0)
    {
        try
        {

            sSOP<T> sop;
            CALL_AND_HANDLE(sop = SOPUtils::primitive(_sop, nmodes), "Failed to simplify the SOP object.");

            size_t _nmodes = get_nmodes(sop);
            if(nmodes == 0){nmodes  = _nmodes;}
            else{ASSERT(nmodes >= _nmodes, "Failed to construct sum_of_product_operator object with the specified number of modes.  The input object has more modes than have been requested.");}

            CALL_AND_HANDLE(return compress_sop_simple(sop, nmodes), "Failed to compute the compressed sop object using the simple strategy.");

        }
        catch(const std::exception& ex)
        {
            std::cerr << ex.what() << std::endl;
            RAISE_EXCEPTION("Failed to compute the simple form of the SOP operator.");
        }
    }
};


}   //namespace utils

#endif

