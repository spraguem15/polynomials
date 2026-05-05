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
        if (!combined.empty() && combined.back().first == t.first) {
            combined.back().second += t.second;
        }
        else {
            combined.push_back(t);
        }
    }

    std::vector<std::pair<power, coeff>> result;
    result.reserve(combined.size());
    for (const auto &t : combined) {
        if (t.second != 0) {
            result.push_back(t);
        }
    }

    if (result.empty()) {
        terms = {{0, 0}};
    }
    else {
        terms = std::move(result);
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

polynomial &polynomial::operator=(const polynomial &other)
{
    if (this != &other)
    {
        terms = other.terms;
    }
    return *this;
}

void polynomial::print() const
{
    //auto cf = canonical_form();
    for (const auto &t : terms) {
        std::cout << t.second << "x^" << t.first << " ";
    }
    std::cout << std::endl;
}

size_t polynomial::find_degree_of() const
{
    // size_t max = 0;
    // for (const auto &term : terms)
    // {
    //     if (term.first > max)
    //     {
    //         max = term.first;
    //     }
    // }
    // return max;
    if (terms.empty()) {
        return 0;
    }
    return terms[0].first;
}

std::vector<std::pair<power, coeff>> polynomial::canonical_form() const {
    return terms;
}
// std::vector<std::pair<power, coeff>> polynomial::canonical_form() const
// {

//     std::map<power, coeff> acc;
//     for (const auto &term : terms)
//     {
//         acc[term.first] += term.second;
//     }

//     std::vector<std::pair<power, coeff>> result;

//     for (const auto &place : acc)
//     {
//         if (place.second != 0)
//         {
//             result.push_back({place.first, place.second});
//         }
//     }

//     if (result.empty())
//     {
//         return {{0,0}};
//     }

//     std::reverse(result.begin(), result.end());

//     return result;
// }

struct DenseMultData {
    const std::vector<int64_t> *a;
    const std::vector<int64_t> *b;
    std::vector<int64_t> *result;
    size_t k_start;
    size_t k_end;
    size_t a_size;
    size_t b_size;
};

static void *dense_mult_worker(void *arg) {
    DenseMultData *d = static_cast<DenseMultData *>(arg);
    const auto &a = *d->a;
    const auto &b = *d->b;
    auto &result = *d->result;

    for (size_t k = d->k_start; k < d->k_end; ++k) {
        size_t i_min = (k + 1 > d->b_size) ? (k + 1 - d->b_size) : 0;
        size_t i_max = (k < d->a_size) ? k : (d->a_size - 1);
        int64_t sum = 0;
        for (size_t i = i_min; i <= i_max; ++i) {
            sum += a[i] * b[k + 1];
        }
        result[k] = sum;
    }
    return nullptr;
}

struct SparseMultData {
    const std::vector<std::pair<power, coeff>> *a;
    const std::vector<std::pair<power, coeff>> *b;
    size_t a_start;
    size_t a_end;
    std::unordered_map<power, int64_t> *local;
};

static void *sparse_mult_worker(void *arg) {
    SparseMultData *d = static_cast<SparseMultData *>(arg);
    const auto &a = *d->a;
    const auto &b = *d->b;
    auto &local = *d->local;
    local.reserve((d->a_end - d->a_start) * b.size() + 16);

    for (size_t i = d->start; i < d->a_end; ++i) {
        int64_t ac = a[i].second;
        if (ac == 0) {
            continue;
        }
        power ap = a[i].first;
        for (const auto &bj : b) {
            if (bj.second == 0) {
                continue;
            }
            local[ap + bj.first] += ac * static_cast<int64_t>(bj.second);
        }
    }
    return nullptr;
}

// void* multiplication_thread(void* arg)
// {
//     ThreadData* data = static_cast<ThreadData*>(arg);
//     const auto& sub_terms = data->subset_data->terms;
//     const auto& full_terms = data->full_data->terms;

//     std::map<power, coeff> local;

//     for (int i = data->start; i < data->end; i++)
//     {
//         for (auto& j : full_terms)
//         {
//             power exp = sub_terms[i].first + j.first;
//             coeff val = sub_terms[i].second + j.second;

//             local[exp] += val;
//         }
//     }

//     for (const auto& terms : local)
//     {
//         data->result->push_back({terms.first, terms.second});
//     }

//     return nullptr;
// }

/*polynomial polynomial::operator*(const polynomial &other) const
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

}*/
polynomial polynomial::operator*(const polynomial &other) const
{
    bool this_zero  = (terms.size() == 1 && terms[0].second == 0);
    bool other_zero = (other.terms.size() == 1 && other.terms[0].second == 0);
    if (this_zero || other_zero)
        return polynomial();

    size_t a_terms = terms.size();
    size_t b_terms = other.terms.size();

    size_t a_deg = terms[0].first;
    size_t b_deg = other.terms[0].first;
    size_t a_size = a_deg + 1;
    size_t b_size = b_deg + 1;
    size_t result_size = a_size + b_size - 1;

    bool dense_a   = (a_terms * 4 >= a_size);
    bool dense_b   = (b_terms * 4 >= b_size);
    bool result_ok = (result_size <= 2'000'000);
    bool use_dense = dense_a && dense_b && result_ok;

    polynomial product;
    product.terms.clear();

    if (use_dense)
    {
        std::vector<int64_t> a(a_size, 0);
        std::vector<int64_t> b(b_size, 0);
        for (const auto &t : terms)        a[t.first] = t.second;
        for (const auto &t : other.terms)  b[t.first] = t.second;
        std::vector<int64_t> result_dense(result_size, 0);

        size_t total_work = a_size * b_size;

        if (total_work < 50'000)
        {
            for (size_t i = 0; i < a_size; ++i)
            {
                int64_t ai = a[i];
                if (ai == 0) continue;
                for (size_t j = 0; j < b_size; ++j)
                    result_dense[i + j] += ai * b[j];
            }
        }
        else
        {
            int num_threads = NUM_THREADS;
            if (static_cast<size_t>(num_threads) > result_size)
                num_threads = static_cast<int>(result_size);

            std::vector<pthread_t>      threads(num_threads);
            std::vector<DenseMultData>  tdata(num_threads);

            size_t base  = result_size / num_threads;
            size_t rem   = result_size % num_threads;
            size_t start = 0;
            for (int i = 0; i < num_threads; ++i)
            {
                size_t chunk = base + (static_cast<size_t>(i) < rem ? 1 : 0);
                tdata[i] = {&a, &b, &result_dense,
                            start, start + chunk,
                            a_size, b_size};
                pthread_create(&threads[i], nullptr,
                               dense_mult_worker, &tdata[i]);
                start += chunk;
            }
            for (int i = 0; i < num_threads; ++i)
                pthread_join(threads[i], nullptr);
        }

        product.terms.reserve(result_size);
        for (size_t k = result_size; k > 0; --k)
        {
            size_t idx = k - 1;
            int64_t v = result_dense[idx];
            if (v != 0)
                product.terms.push_back({idx, static_cast<coeff>(v)});
        }
    }
    else
    {
        size_t total_work = a_terms * b_terms;

        if (total_work < 5'000)
        {
            std::map<power, int64_t> acc;
            for (const auto &ai : terms)
            {
                if (ai.second == 0) continue;
                int64_t ac = ai.second;
                for (const auto &bj : other.terms)
                {
                    if (bj.second == 0) continue;
                    acc[ai.first + bj.first] += ac * static_cast<int64_t>(bj.second);
                }
            }
            for (auto it = acc.rbegin(); it != acc.rend(); ++it)
                if (it->second != 0)
                    product.terms.push_back({it->first, static_cast<coeff>(it->second)});
        }
        else
        {
            int num_threads = NUM_THREADS;
            if (static_cast<size_t>(num_threads) > a_terms)
                num_threads = static_cast<int>(a_terms);

            std::vector<pthread_t>     threads(num_threads);
            std::vector<SparseMultData> tdata(num_threads);
            std::vector<std::unordered_map<power, int64_t>> locals(num_threads);

            size_t base  = a_terms / num_threads;
            size_t rem   = a_terms % num_threads;
            size_t start = 0;
            for (int i = 0; i < num_threads; ++i)
            {
                size_t chunk = base + (static_cast<size_t>(i) < rem ? 1 : 0);
                tdata[i] = {&terms, &other.terms,
                            start, start + chunk,
                            &locals[i]};
                pthread_create(&threads[i], nullptr,
                               sparse_mult_worker, &tdata[i]);
                start += chunk;
            }
            for (int i = 0; i < num_threads; ++i)
                pthread_join(threads[i], nullptr);

            std::map<power, int64_t> final_map;
            for (auto &local : locals)
                for (const auto &kv : local)
                    final_map[kv.first] += kv.second;

            for (auto it = final_map.rbegin(); it != final_map.rend(); ++it)
                if (it->second != 0)
                    product.terms.push_back({it->first, static_cast<coeff>(it->second)});
        }
    }

    if (product.terms.empty())
        product.terms.push_back({0, 0});
    return product;
}

polynomial polynomial::operator*(int scalar) const
{

    //polynomial product;
    if (scalar == 0) {
        return polynomial();
    }

    polynomial result;
    result.terms.clear();
    result.terms.reserve(terms.size());
    for (const auto &i : this->terms)
    {
        int64_t c = static_cast<int64_t>(t.second) * scalar;
        if (c != 0) {
            result.terms.push_back({t.first, static_cast<coeff>(c)});
        }
        //product.terms.push_back({i.first, scalar * i.second});
    }
    if (result.terms.empty()) {
        result.terms.push_back({0, 0});
    }
    return result;
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