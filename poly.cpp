#include <poly.h>
#include <vector>
#include <utility>
#include <cstddef>
#include <map>
#include <algorithm>
#include <iostream>
#include <stdexcept>

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
    auto cf = canonical_form();
    for (const auto &term : cf) {
        std::cout << term.second << "x^" << term.first << " ";
    }
    std::cout << std::endl;
}

polynomial &polynomial::operator=(const polynomial &other)
{
    if (this != &other)
    {
        terms = other.terms;
    }
    return *this;
}

size_t polynomial::find_degree_of() const
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
        return {{0,0}};
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


polynomial polynomial::operator+(const polynomial &other) const {
    polynomial result;
    result.terms.clear();
    result.terms.insert(result.terms.end(), terms.begin(), terms.end());
    result.terms.insert(result.terms.end(), other.terms.begin(), other.terms.end());

    return result;
}

polynomial polynomial::operator+(int value) const {
    polynomial result = *this;
    result.terms.push_back({0, value});
    return result;
}

polynomial operator+(int value, const polynomial &poly) {
    return poly + value;
}

polynomial polynomial::operator%(const polynomial &divisor) const {
    auto r = this->canonical_form();
    auto d = divisor.canonical_form();

    if(d.size() == 1 && d[0].second == 0) {
        throw std::invalid_argument("Modulo by a zero polynomial");
    }

    while (!r.empty() && r[0].first >= d[0].first) {
        power power_diff = r[0].first - d[0].first;
        coeff coeff_ratio = r[0].second / d[0].second;

        std::vector<std::pair<power, coeff>> temp;

        for (const auto &term : d) {
            temp.push_back({term.first + power_diff, term.second * coeff_ratio});
        }
        std::map<power, coeff> acc;

        for (const auto &term : r) {
            acc[term.first] += term.second;
        }

        r.clear();

        for (const auto &p : acc) {
            if (p.second != 0) {
                r.push_back({p.first, p.second});
            }
        }

        std::sort(r.begin(), r.end(), [](auto &a, auto &b)
        {
            return a.first > b.first;
        });
    }
    if (r.empty()) {
        return polynomial();
    }
    polynomial result;
    result.terms = r;
    return result;
}