#include "poly.h"
#include <vector>
#include <utility>
#include <cstddef>
#include <map>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <pthread.h>
#include <cmath>
#include <unordered_map>

using power = size_t;
using coeff = int;

static constexpr int NUM_THREADS = 8;

// struct ThreadData {
//     const polynomial* subset_data;
//     const polynomial* full_data;
//     int start;
//     int end;
//     std::vector<std::pair<power, coeff>>* result;
// };

void polynomial::canonicalize() {
    if(terms.empty()) {
        terms.push_back({0,0});
        return;
    }

    std::sort(terms.begin(), terms.end(), [](const std::pair<power, coeff> &a, const std::pair<power, coeff> &b)
                {
                    return a.first > b.first;
                });
    
    std::vector<std::pair<power, coeff>> combined;
    combined.reserve(terms.size());

    for (const auto &t : terms) {
        //
    }
}

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

void* multiplication_thread(void* arg)
{
    ThreadData* data = static_cast<ThreadData*>(arg);
    const auto& sub_terms = data->subset_data->terms;
    const auto& full_terms = data->full_data->terms;

    std::map<power, coeff> local;

    for (int i = data->start; i < data->end; i++)
    {
        for (auto& j : full_terms)
        {
            power exp = sub_terms[i].first + j.first;
            coeff val = sub_terms[i].second + j.second;

            local[exp] += val;
        }
    }

    for (const auto& terms : local)
    {
        data->result->push_back({terms.first, terms.second});
    }

    return nullptr;
}

polynomial polynomial::operator*(const polynomial &other) const
{
    polynomial product;
    int terms = this->terms.size();
    int other_terms = other.terms.size();

    if ((terms * other_terms) >= 10000)
    {
        int num_threads = 12;
        int thread_workload = terms / num_threads;
        int last_workload = thread_workload + (terms % num_threads);

        std::vector<pthread_t> threads(num_threads);
        std::vector<ThreadData> thread_data(num_threads);
        std::vector<std::vector<std::pair<power, coeff>>> temp(num_threads);

        int start = 0;
        int end, remainder;

        for (int i = 0; i < num_threads; i++) 
        {

            remainder = thread_workload;
            if (i == num_threads - 1)
            {
                remainder += last_workload;
            }
            end = start + remainder;
            temp[i].reserve((end - start) * other_terms);
            thread_data[i] = {this, &other, start, end, &temp[i]};

            pthread_create(&threads[i], nullptr, multiplication_thread, &thread_data[i]);

            start = end;

        }

        for (int i = 0; i < num_threads; i++)
        {
            pthread_join(threads[i], nullptr);
        }

        for (int i = 0; i < num_threads; i++)
        {
            auto& thread_results = *(thread_data[i].result);
            for (const auto& term : thread_results)
            {
                product.terms.push_back(term);  
            }
        }
        return product;
    }
    else
    {
        for (const auto &i : this->terms)
        {
            for (const auto &j : other.terms)
            {
                product.terms.push_back({i.first + j.first, i.second * j.second});
            }
        }
        return product;
    }

}

polynomial polynomial::operator*(int scalar) const
{

    polynomial product;
    for (const auto &i : this->terms)
    {
        product.terms.push_back({i.first, scalar * i.second});
    }
    return product;
}

polynomial operator*(int scalar, const polynomial &other)
{
    return other * scalar;
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

        for (const auto &term : d)
        {
            acc[term.first + power_diff] -= term.second * coeff_ratio;
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