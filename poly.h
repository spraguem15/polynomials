#ifndef POLY_H
#define POLY_H

#include <vector>
#include <utility>
#include <cstddef>

using power = size_t;
using coeff = int;

class polynomial
{
private:
    std::vector<std::pair<power, coeff>> terms;
    void canonicalize();

public:
    polynomial();

    template <typename Iter>
    polynomial(Iter begin, Iter end)
    {
        for (Iter it = begin; it != end; ++it) {
            terms.push_back(*it);
        }
        canonicalize();
    }

    polynomial(const polynomial &other);
    polynomial &operator=(const polynomial &other);

    void print() const;

    polynomial operator+(const polynomial &other) const;
    polynomial operator+(int value) const;

    polynomial operator*(const polynomial &other) const;
    polynomial operator*(int value) const;

    polynomial operator%(const polynomial &other) const;

    size_t find_degree_of() const;
    std::vector<std::pair<power, coeff>> canonical_form() const;
};

polynomial operator+(int value, const polynomial &poly);
polynomial operator*(int value, const polynomial &poly);

#endif
