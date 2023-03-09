#ifndef UTILS_SYSTEM_INFORMATION_HPP
#define UTILS_SYSTEM_INFORMATION_HPP

namespace utils
{

class mode_data
{
public:
    mode_data() : m_lhd(1), m_is_fermionic(false){}
    mode_data(size_t d) : m_lhd(d), m_is_fermionic(false){}
    mode_data(size_t d, bool is_fermionic) : m_lhd(d), m_is_fermionic(is_fermionic){}

    mode_data(const mode_data& o) = default;
    mode_data(mode_data&& o) = default;

    mode_data& operator=(const mode_data& o) = default;
    mode_data& operator=(mode_data&& o) = default;

    void make_fermionic(){m_is_fermionic = true;    m_lhd = 2;}

    const bool& fermionic() const{return m_is_fermionic;}
    bool& fermionic(){return m_is_fermionic;}

    const size_t& lhd() const{return m_lhd;}
    size_t& lhd(){return m_lhd;}
protected:
    size_t m_lhd;           //local hilbert space dimension
    bool m_is_fermionic;    //whether or not the mode is fermionic
};


class composite_mode
{   
public:
    composite_mode(){}
    composite_mode(const composite_mode& o) = default;
    composite_mode(composite_mode&& o) = default;

    composite_mode& operator=(const composite_mode& o) = default;
    composite_mode& operator=(composite_mode&& o) = default;
protected:
    std::vector<size_t> m_primitive_mode_indices;
    std::shared_ptr<occupation_number_basis> m_composite_basis;
};

class system_modes
{
public:
    using iterator = typename std::vector<mode_data>::iterator;
    using const_iterator = typename std::vector<mode_data>::const_iterator;
    using reverse_iterator = typename std::vector<mode_data>::reverse_iterator;
    using const_reverse_iterator = typename std::vector<mode_data>::const_reverse_iterator;

public:
    system_modes(size_t N) : m_primitive_modes(N) {}
    system_modes(size_t N, size_t d) : m_primitive_modes(N)
    {
        for(auto& mode : m_primitive_modes){mode.lhd() = d;}
    }
    system_modes(const system_modes& o) = default;
    system_modes(system_modes&& o) = default;

    system_modes& operator=(const system_modes& o) = default;
    system_modes& operator=(system_modes&& o) = default;

    size_t nmodes() const{return m_primitive_modes.size();}

    void resize(size_t N)
    {
        if(N >= nmodes())
        {
            m_primitive_modes.resize(N);
        }
        else
        {
            clear();
            m_primitive_modes.resize(N);
        }
    }


    mode_data& operator[](size_t i)
    {
        ASSERT(i < m_primitive_modes.size(), "Index out of bounds.");
        return m_primitive_modes[i];
    }

    const mode_data& operator[](size_t i) const
    {
        ASSERT(i < m_primitive_modes.size(), "Index out of bounds.");
        return m_primitive_modes[i];
    }

    mode_data& mode(size_t i)
    {
        ASSERT(i < m_primitive_modes.size(), "Index out of bounds.");
        return m_primitive_modes[i];
    }

    const mode_data& mode(size_t i) const
    {
        ASSERT(i < m_primitive_modes.size(), "Index out of bounds.");
        return m_primitive_modes[i];
    }

    void clear() noexcept
    {
        m_primitive_modes.clear();
        m_composite_modes.clear();
        m_bound_primitive_modes.clear();
    }

public:
    iterator begin() {  return iterator(m_primitive_modes.begin());  }
    iterator end() {  return iterator(m_primitive_modes.end());  }
    const_iterator begin() const {  return const_iterator(m_primitive_modes.begin());  }
    const_iterator end() const {  return const_iterator(m_primitive_modes.end());  }

    reverse_iterator rbegin() {  return reverse_iterator(m_primitive_modes.rbegin());  }
    reverse_iterator rend() {  return reverse_iterator(m_primitive_modes.rend());  }
    const_reverse_iterator rbegin() const {  return const_reverse_iterator(m_primitive_modes.rbegin());  }
    const_reverse_iterator rend() const {  return const_reverse_iterator(m_primitive_modes.rend());  }

protected:
    std::vector<mode_data> m_primitive_modes;
    std::vector<composite_mode> m_composite_modes;
    std::list<size_t> m_bound_primitive_modes;
};

inline system_modes fermionic_system(size_t N)
{
    system_modes sys(N);
    for(auto& mode : sys)
    {
        mode.make_fermionic();
    }
    return sys;
}



}

#endif

