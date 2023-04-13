#ifndef UTILS_OPERATOR_GEN_AUTO_SUM_OF_PRODUCT_HPP
#define UTILS_OPERATOR_GEN_AUTO_SUM_OF_PRODUCT_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>

#include "sSOP.hpp"
#include "compressedSOP.hpp"
#include "system_information.hpp"

#include "bipartite_graph.hpp"

namespace utils
{


class bipartitionSOP
{
public:
    //a function for taking a bipartitioning of the Hamiltonain 
    template <typename T>
    static inline bipartite_graph<sPOP, T> form_bipartite_graph(const sSOP<T>& sop, const std::vector<size_t>& modes)
    {
        //sort the modes vector into smodes
        std::vector<size_t> smodes(modes.begin(), modes.end());
        std::sort(smodes.begin(), smodes.end());

        //start by getting the bipartitioned sets of operators
        bipartite_graph<sPOP, T> bpg;
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

            //now for all of the nodes in _V we push all connected nodes 
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





class SOPUtils
{
public:
    template <typename T>
    static inline sSOP<T> primitive(const sSOP<T>& sop, const system_modes& sys_info)
    {
        //check that the size of the system implied by the SOP is consistent with the system information
        size_t _nmodes = get_nmodes(sop);
        size_t nmodes = sys_info.nmodes();
        ASSERT(nmodes >= _nmodes, "Failed to construct sum_of_product_operator object with the specified number of modes.  The input object has more modes than have been requested.");

        //determine the modes that have been specified as fermionic due to the operator
        std::vector<bool> is_fermion_mode(nmodes);      std::fill(is_fermion_mode.begin(), is_fermion_mode.end(), false);
        ASSERT(set_is_fermionic_mode(sop, is_fermion_mode), "Failed to determine if operators were fermionic from SOP.  Inconsistent definitions throughout.");

        //now check that the fermionic modes specified by the operator is consistent with the system specification.  Here we only check 
        bool contains_fermionic_operator = false;
        for(size_t i = 0; i < _nmodes; ++i)
        {
            ASSERT(is_fermion_mode[i] == sys_info[i].fermionic(), "The system information and SOP information about which modes are fermionic are inconsistent.");
            if(is_fermion_mode[i]){contains_fermionic_operator = true;}
        }
        for(size_t i = _nmodes; i < nmodes; ++i)
        {
            //and set the remaining modes that the sop does not provide any information about to be fermionic 
            is_fermion_mode[i] = sys_info[i].fermionic();
        }

        //now we build a new sop operator that does not contain any fermionic operators - e.g. perform the required Jordan wigner string mappings to remove all traces of operators
        sSOP<T> mapped_sop;
        CALL_AND_HANDLE(map_operator(sop, is_fermion_mode, mapped_sop), "Failed to map sum of product operator");

        //now that we have mapped all fermionic operators to qubit operators we no longer have to worry about the anticommutation relations and we can begin to simplify the 
        sSOP<T> simplified_sop;
        CALL_AND_HANDLE(simplify_operator(mapped_sop, is_fermion_mode, simplified_sop), "Failed to simplify the form of the sum of product operator.");

        return simplified_sop;
    }

    template <typename T>
    static inline sSOP<T> primitive(const sSOP<T>& sop, size_t nmodes = 0)
    {
        //check that the size of the system implied by the SOP is consistent with the specified number of modes
        size_t _nmodes = get_nmodes(sop);
        if(nmodes == 0){nmodes  = _nmodes;}
        else{ASSERT(nmodes >= _nmodes, "Failed to construct sum_of_product_operator object with the specified number of modes.  The input object has more modes than have been requested.");}

        //determine the modes that have been specified as fermionic due to the operator
        std::vector<bool> is_fermion_mode(nmodes);      std::fill(is_fermion_mode.begin(), is_fermion_mode.end(), false);
        ASSERT(set_is_fermionic_mode(sop, is_fermion_mode), "Failed to determine if operators were fermionic from SOP.  Inconsistent definitions throughout.");


        //now we build a new sop operator that does not contain any fermionic operators - e.g. perform the required Jordan wigner string mappings to remove all traces of operators
        sSOP<T> mapped_sop;
        CALL_AND_HANDLE(map_operator(sop, is_fermion_mode, mapped_sop), "Failed to map sum of product operator");

        //now we go through and we map each 
        sSOP<T> simplified_sop;
        CALL_AND_HANDLE(simplify_operator(mapped_sop, is_fermion_mode, simplified_sop), "Failed to simplify the form of the sum of product operator.");

        return simplified_sop;
    }

protected:
    template <typename T> 
    static inline void simplify_operator(const sSOP<T>& sop, const std::vector<bool>& is_fermi, sSOP<T>& _op)
    {
        //now we perform the JW mapping for this operator
        for(const auto& pop : sop)
        {
            sNBO<T> term;    term.coeff() = pop.coeff();
            //only add in the term if the coefficient is non-zero
            if(term.coeff() != T(0))
            {
                //check if the mode contains a fermionic operator at all.  If it doesn't then we just directly insert it into the new operator
                if(simplify_string_operator(pop.pop(), is_fermi, term.pop())){term.coeff() = -pop.coeff();}
                _op += term;
            }
        }

    }

    static inline bool simplify_string_operator(const sPOP& in, const std::vector<bool>& is_fermi, sPOP& out)
    {
        size_t nmodes = is_fermi.size();
        std::vector<std::list<std::string>> terms(nmodes);

        //now we are free to reorder the operators acting on different modes but we need to preserve the ordering of operators acting on a single mode.
        for(const auto& op : in)
        {
            terms[op.mode()].push_back(op.op());
        }


        bool flip_coeff = false;
        std::vector<size_t> njw(nmodes);    std::fill(njw.begin(), njw.end(), 0);
        //now we attempt to simplify the JW string operators acting on a given mode
        for(size_t i = 0; i < nmodes; ++i)
        {
            if(is_fermi[i])
            {
                njw[i] = 0;
                bool attempt_simplify = true;
                for(const auto& str : terms[i])
                {
                    if(str == std::string("jw"))
                    {
                        ++njw[i];
                    }
                    else if( (str != std::string("n") && str != std::string("a") && str != std::string("adag")))
                    {
                        attempt_simplify = false;
                    }
                }

                //if the number of jw strings is greater than 1 then we pull the JW strings to the front
                //and cancel in pairs
                if(attempt_simplify && njw[i] != 0)
                {
                    size_t nc = 0;
                    for(const auto& str : terms[i])
                    {
                        if(str == std::string("a") || str == std::string("adag"))
                        {
                            ++nc;
                        }
                        else if(str == std::string("jw"))
                        {
                            //if there are an odd number of creation and annihilation operators then we use the fact that bringing the JW string through a creation or annihilation operator flips the sign
                            if(nc % 2 == 1)
                            {
                                flip_coeff = !flip_coeff;
                            }
                        }
                    }
                }
            }
        }


        //and finally we create the new string product operator to return
        for(size_t i = 0; i < nmodes; ++i)
        {
            if(is_fermi[i])
            {
                bool njwadd = false;
                for(const auto& str : terms[i])
                {
                    if(str != std::string("jw"))
                    {
                        //if this is the first non-jw term we are adding to this mode
                        if(!njwadd)
                        {
                            if(njw[i] % 2 == 1)
                            {
                                //now we check the first operator
                                if(str == std::string("n") || str == std::string("adag")){flip_coeff = !flip_coeff;}
                                else if(str == std::string("a")){}
                                else
                                {
                                    out*=sOP(std::string("jw"), i);
                                }
                            }
                        }
                        njwadd = true;
                        out *= sOP(str, i);
                    }
                }
                if(!njwadd)
                {
                    if(njw[i] % 2 == 1)
                    {
                        out*=sOP(std::string("jw"), i);
                    }
                }
            }
            else
            {
                for(const auto& str : terms[i])
                {
                    out *= sOP(str, i);
                }
            }
        }
        return flip_coeff;
    }


    template <typename T> 
    static inline void map_operator(const sSOP<T>& sop, const std::vector<bool>& is_fermi, sSOP<T>& _op)
    {

        bool contains_fermionic_operator = false;
        size_t nfermionic = 0;
        for(const auto& isf : is_fermi)
        {
            if(isf){contains_fermionic_operator = true; ++nfermionic;}
        }

        std::vector<size_t> fermion_modes(nfermionic);

        size_t counter = 0;
        for(size_t i = 0; i < is_fermi.size(); ++i)
        {
            if(is_fermi[i]){fermion_modes[counter] = i; ++counter;}
        }


        sSOP<T> _sop;

        //here we are just doing the processing required to handle the fermionic modes
        if(contains_fermionic_operator)
        {
            //now we perform the JW mapping for this operator
            for(const auto& pop : sop)
            {
                sNBO<T> term;    term.coeff() = pop.coeff();
                if(pop.coeff() != T(0))
                {
                    //check if the mode contains a fermionic operator at all.  If it doesn't then we just directly insert it into the new operator
                    bool term_has_fermion_operator = false;
                    for(const auto& op : pop)
                    {
                        if(op.fermionic())
                        {
                            term_has_fermion_operator = true;
                            break;
                        }
                    }

                    if(term_has_fermion_operator)
                    {
                        //now iterate over the terms in this operator and if the 
                        for(const auto& op : pop)
                        {
                            if(op.fermionic())
                            {
                                sOP oterm;
                                int include_JW_string = jordan_wigner(op, oterm);

                                if(include_JW_string == -1)
                                {
                                    //add in the sigma_z terms for all fermion modes before the current mode
                                    for(size_t i = 0; i < op.mode(); ++i)
                                    {
                                        if(is_fermi[i])
                                        {
                                            term*=sOP({"jw", i});
                                        }
                                    }
                                    term*=oterm;
                                }
                                else if(include_JW_string == 1)
                                {
                                    term*=oterm;
                                    for(size_t i = 0; i < op.mode(); ++i)
                                    {
                                        size_t l = op.mode()-(i+1);
                                        if(is_fermi[l])
                                        {
                                            term*=sOP({"jw", l});
                                        }
                                    }
                                }
                                else
                                {
                                    term*=oterm;
                                }
                            }
                            else
                            {
                                term *= op;
                            }
                        }
                        _sop += term;
                    }
                    else
                    {
                        _sop += pop;
                    }
                }
            }
        }
        else
        {
            for(const auto& pop : sop)
            {
                if(pop.coeff() != T(0))
                {
                    _sop += pop;
                }
            }
        }
    
        _op = _sop;
    }

    static inline int jordan_wigner(const sOP& in, sOP& out)
    {
        std::string label(in.op());
        io::remove_whitespace_and_to_lower(label);

        if(label == "c" || label == "a"|| label == "f")
        {
            out = sOP({"a", in.mode()});
            return -1;
        }
        else if(label == "cdag" || label == "c^\\dagger"  || label == "adag" || label == "a^\\dagger"  || label == "fdag" || label == "f^\\dagger")
        {
            out = sOP({"adag", in.mode()});
            return 1;
        }
        else if(label == "n")
        {
            out = sOP({"n", in.mode()});
            return 0;
        }
        else
        {
            //if it isn't a standard fermionic string then we just assume that it transforms as the number operator 
            out = sOP({label, in.mode()});
            return 0;
        }
    }


    template <typename T> 
    static inline bool set_is_fermionic_mode(const sSOP<T>& sop, std::vector<bool>& is_fermion_mode)
    {
        std::fill(is_fermion_mode.begin(), is_fermion_mode.end(), false);
        std::vector<bool> mode_properties_set(is_fermion_mode.size());  std::fill(mode_properties_set.begin(), mode_properties_set.end(), false);

        //first verify that the input sSOP is valid
        for(const auto& pop : sop)
        {
            for(const auto& op : pop)
            {
                //if we have set the mode properties of this mode already ensure that the new properties we have obtained from this operator are consistent.
                if(mode_properties_set[op.mode()])
                {
                    if(op.fermionic() != is_fermion_mode[op.mode()]){return false;}
                }
                else
                {
                    if(op.fermionic()){is_fermion_mode[op.mode()] = true;}
                    mode_properties_set[op.mode()] = true;
                }
            }
        }
        return true;
    }


    template <typename T> 
    static inline size_t get_nmodes(const sSOP<T>& sop)
    {
        size_t _nmodes = 0;
        for(const auto& pop : sop)
        {
            for(const auto& op : pop)
            {
                if(op.mode()+1 > _nmodes)
                {
                    _nmodes = op.mode()+1;
                }
            }
        }
        return _nmodes;
    }

};


class AutoSOP
{
protected:
    template <typename T>
    static inline compressedSOP<T> compress_sop_simple(const sSOP<T>& sop, size_t nmodes)
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

