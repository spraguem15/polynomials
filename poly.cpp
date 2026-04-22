#include <poly.h>
#include <vector>
#include <utility>
#include <cstddef>
#include <map>
#include <algorithm>

using power = size_t;
using coeff = int;

polynomial::polynomial()
{
    terms.push_back({0,0});
}

polynomial::polynomial(const polynomial &other)
{
    terms = other.terms;
}

void polynomial::print() const
{
    
}

polynomial &polynomial::operator=(const polynomial &other)
{
    if (this != &other)
    {
        terms = other.terms;
    }
    return *this;
}

size_t polynomial::find_degree_of()
{
    size_t max = 0;
    for (const auto &term : terms)
    {
        if (term.first > max)
        {
            max = term.first;
        }
    }
    return max;
}

std::vector<std::pair<power, coeff>> polynomial::canonical_form() const
{

    std::map<power, coeff> acc;
    for (const auto &term : terms)
    {
        acc[term.first] += term.second;
    }

    std::vector<std::pair<power, coeff>> result;

    for (const auto &place : acc)
    {
        if (place.second != 0)
        {
            result.push_back({place.first, place.second});
        }
    }

    if (result.empty())
    {
        result.push_back({0,0});
    }

    std::reverse(result.begin(), result.end());

    return result;
}


polynomial polynomial::operator*(const polynomial &other)
{
    polynomial product;

    for (const auto &i : this->terms)
    {
        for (const auto &j : other.terms)
        {
            if (i.first == 0 && i.second == 0)
            {
                product.terms.push_back({0,0});
                break;
            }
            product.terms.push_back({i.first + j.first, i.second * j.second});
        }
    }
    return product;
}
