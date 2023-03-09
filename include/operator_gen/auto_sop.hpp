#ifndef UTILS_OPERATOR_GEN_AUTO_SUM_OF_PRODUCT_HPP
#define UTILS_OPERATOR_GEN_AUTO_SUM_OF_PRODUCT_HPPR

#include <vector>
#include <algorithm>
#include <unordered_map>

#include "sSOP.hpp"
#include "system_information.hpp"

namespace utils
{


class bipartitionSOP
{
public:
    //a function for taking a bipartitioning of the Hamiltonain 
    template <typename T>
    static inline void form_bipartite_graph(const sSOP<T>& sop, const std::vector<size_t>& modes)
    {
        //sort the modes vector into smodes
        std::vector<size_t> smodes(modes.begin(), modes.end());
        std::sort(smodes.begin(), smodes.end());

        //start by getting the bipartitioned sets of operators
        std::vector<sPOP> U, V;
        std::list<std::tuple<size_t, size_t, T>> E;
        {
            size_t nu = 0, nv = 0;
            std::unordered_map<std::string, std::pair<sPOP, size_t>> _U, _V;
            for(const auto& NBO : sop)
            {
                const auto& pop = NBO.pop();
                sPOP _u, _v;

                //extract tdb
                bipartition_term(pop, smodes, _u, _v);

                std::string ukey = get_key(_u);

                size_t _ui, _vi;
                //if the key isn't present in the map
                auto Ui = _U.find(ukey);
                if(Ui == _U.end())
                {
                    _U.insert(std::make_pair(ukey, std::make_pair(_u, nu)));
                    _ui = nu;
                    ++nu;
                }
                else
                {
                    _ui = Ui->second.second;
                }

                std::string vkey = get_key(_v);
                auto Vi = _V.find(vkey);
                if(Vi == _V.end())
                {
                    _V.insert(std::make_pair(vkey, std::make_pair(_v, nv)));
                    _vi = nv;
                    ++nv;
                }
                else
                {
                    _vi = Vi->second.second;
                }
                E.push_back(std::make_tuple(_ui, _vi, NBO.coeff()));
            }

            U.resize(_U.size());
            for(const auto& t : _U)
            {
                U[t.second.second] = t.second.first;
            }
            V.resize(_V.size());
            for(const auto& t : _V)
            {
                V[t.second.second] = t.second.first;
            }

        }

        std::cout << "U" << std::endl;
        for(size_t i = 0; i < U.size(); ++i)
        {
            std::cout << i << " " << U[i] <<std::endl;
        }

        std::cout << "V" << std::endl;
        for(size_t i = 0; i < V.size(); ++i)
        {
            std::cout << i << " " << V[i] <<std::endl;
        }

        std::cout << "E" << std::endl;
        for(const auto& t : E)
        {
            std::cout << std::get<0>(t) << " " << std::get<1>(t) << " " << std::get<2>(t) << std::endl;
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

    static inline std::string get_key(const sPOP& l)
    {
        std::ostringstream ss;
        ss << l;
        return ss.str();
    }
};

class AutoSOP
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

        //now we go through and we amp each 
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
            //check if the mode contains a fermionic operator at all.  If it doesn't then we just directly insert it into the new operator
            if(simplify_string_operator(pop.pop(), is_fermi, term.pop())){term.coeff() = -pop.coeff();}
            _op += term;
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
        else
        {
            _sop = sop;
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


}   //namespace utils

#endif

