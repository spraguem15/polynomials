#include "poly.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <pthread.h>
#include <stdexcept>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

static constexpr int MAX_THREADS = 8;
static constexpr size_t NTT_THRESHOLD = 4096;
static constexpr size_t SPARSE_THREAD_THRESHOLD = 20000;

static int detect_thread_count()
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n < 1) n = 1;
    if (n > MAX_THREADS) n = MAX_THREADS;
    return static_cast<int>(n);
}

static bool is_zero_terms(const std::vector<std::pair<power, coeff>>& v)
{
    return v.size() == 1 && v[0].first == 0 && v[0].second == 0;
}


template <int MOD, int G>
static int mod_pow_ll(long long a, long long e)
{
    long long r = 1;
    while (e > 0) {
        if (e & 1) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1;
    }
    return static_cast<int>(r);
}

template <int MOD, int G>
static void ntt(std::vector<int>& a, bool invert)
{
    const int n = static_cast<int>(a.size());

    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        int wlen = mod_pow_ll<MOD, G>(G, (MOD - 1) / len);
        if (invert) wlen = mod_pow_ll<MOD, G>(wlen, MOD - 2);

        for (int i = 0; i < n; i += len) {
            long long w = 1;
            for (int j = 0; j < len / 2; ++j) {
                int u = a[i + j];
                int v = static_cast<int>(a[i + j + len / 2] * w % MOD);
                int x = u + v;
                if (x >= MOD) x -= MOD;
                int y = u - v;
                if (y < 0) y += MOD;
                a[i + j] = x;
                a[i + j + len / 2] = y;
                w = w * wlen % MOD;
            }
        }
    }

    if (invert) {
        int inv_n = mod_pow_ll<MOD, G>(n, MOD - 2);
        for (int& x : a) x = static_cast<int>(static_cast<long long>(x) * inv_n % MOD);
    }
}

template <int MOD, int G>
static std::vector<int> convolution_mod(const std::vector<std::pair<power, coeff>>& a_terms,
                                        const std::vector<std::pair<power, coeff>>& b_terms,
                                        size_t result_size)
{
    size_t n = 1;
    while (n < result_size) n <<= 1;

    std::vector<int> a(n, 0), b(n, 0);
    for (const auto& t : a_terms) {
        long long x = t.second % MOD;
        if (x < 0) x += MOD;
        a[t.first] = static_cast<int>(x);
    }
    for (const auto& t : b_terms) {
        long long x = t.second % MOD;
        if (x < 0) x += MOD;
        b[t.first] = static_cast<int>(x);
    }

    ntt<MOD, G>(a, false);
    ntt<MOD, G>(b, false);
    for (size_t i = 0; i < n; ++i) {
        a[i] = static_cast<int>(static_cast<long long>(a[i]) * b[i] % MOD);
    }
    ntt<MOD, G>(a, true);
    a.resize(result_size);
    return a;
}

static std::vector<std::pair<power, coeff>> ntt_multiply_dense(
    const std::vector<std::pair<power, coeff>>& a_terms,
    const std::vector<std::pair<power, coeff>>& b_terms,
    size_t result_size)
{
    static constexpr long long M1 = 998244353LL;   // primitive root 3
    static constexpr long long M2 = 1004535809LL;  // primitive root 3
    static constexpr long long INV_M1_MOD_M2 = 669690699LL; // inverse of M1 modulo M2

    std::vector<int> r1 = convolution_mod<998244353, 3>(a_terms, b_terms, result_size);
    std::vector<int> r2 = convolution_mod<1004535809, 3>(a_terms, b_terms, result_size);

    const long long MOD_PRODUCT = M1 * M2;
    std::vector<std::pair<power, coeff>> out;
    out.reserve(result_size / 2 + 1);

    for (size_t k = result_size; k > 0; --k) {
        size_t idx = k - 1;
        long long a = r1[idx];
        long long b = r2[idx];
        long long t = (b - a) % M2;
        if (t < 0) t += M2;
        t = (t * INV_M1_MOD_M2) % M2;
        long long x = a + M1 * t; // 0 <= x < M1*M2, fits in int64_t
        if (x > MOD_PRODUCT / 2) x -= MOD_PRODUCT;
        if (x != 0) out.push_back({idx, static_cast<coeff>(x)});
    }
    return out;
}


struct SparseMultData {
    const std::vector<std::pair<power, coeff>>* a;
    const std::vector<std::pair<power, coeff>>* b;
    size_t begin;
    size_t end;
    std::unordered_map<power, long long>* local;
};

static void* sparse_mult_worker(void* arg)
{
    SparseMultData* d = static_cast<SparseMultData*>(arg);
    const auto& a = *d->a;
    const auto& b = *d->b;
    auto& local = *d->local;

    size_t expected = (d->end - d->begin) * b.size();
    if (expected > 0 && expected < 2000000) local.reserve(expected * 2);

    for (size_t i = d->begin; i < d->end; ++i) {
        const power ap = a[i].first;
        const long long ac = a[i].second;
        for (const auto& bt : b) {
            local[ap + bt.first] += ac * static_cast<long long>(bt.second);
        }
    }
    return nullptr;
}


void polynomial::canonicalize()
{
    if (terms.empty()) {
        terms.push_back({0, 0});
        return;
    }

    std::sort(terms.begin(), terms.end(),
              [](const std::pair<power, coeff>& a, const std::pair<power, coeff>& b) {
                  return a.first > b.first;
              });

    std::vector<std::pair<power, coeff>> out;
    out.reserve(terms.size());

    for (const auto& t : terms) {
        if (!out.empty() && out.back().first == t.first) {
            long long sum = static_cast<long long>(out.back().second) + t.second;
            out.back().second = static_cast<coeff>(sum);
            if (out.back().second == 0) out.pop_back();
        } else if (t.second != 0) {
            out.push_back(t);
        }
    }

    if (out.empty()) out.push_back({0, 0});
    terms.swap(out);
}

polynomial::polynomial() : terms{{0, 0}} {}

polynomial::polynomial(const polynomial& other) : terms(other.terms) {}

polynomial& polynomial::operator=(const polynomial& other)
{
    if (this != &other) terms = other.terms;
    return *this;
}

void polynomial::print() const
{
    for (const auto& t : terms) std::cout << t.second << "x^" << t.first << " ";
    std::cout << '\n';
}

size_t polynomial::find_degree_of() const
{
    return terms.empty() ? 0 : terms[0].first;
}

std::vector<std::pair<power, coeff>> polynomial::canonical_form() const
{
    return terms;
}

polynomial polynomial::operator+(const polynomial& other) const
{
    polynomial result;
    result.terms.clear();
    result.terms.reserve(terms.size() + other.terms.size());

    size_t i = 0, j = 0;
    while (i < terms.size() && j < other.terms.size()) {
        if (terms[i].first > other.terms[j].first) {
            result.terms.push_back(terms[i++]);
        } else if (terms[i].first < other.terms[j].first) {
            result.terms.push_back(other.terms[j++]);
        } else {
            long long s = static_cast<long long>(terms[i].second) + other.terms[j].second;
            if (s != 0) result.terms.push_back({terms[i].first, static_cast<coeff>(s)});
            ++i;
            ++j;
        }
    }
    while (i < terms.size()) result.terms.push_back(terms[i++]);
    while (j < other.terms.size()) result.terms.push_back(other.terms[j++]);

    if (result.terms.empty()) result.terms.push_back({0, 0});
    return result;
}

polynomial polynomial::operator+(int value) const
{
    if (value == 0) return *this;
    polynomial c;
    c.terms[0] = {0, value};
    return *this + c;
}

polynomial operator+(int value, const polynomial& poly)
{
    return poly + value;
}

polynomial polynomial::operator*(int value) const
{
    if (value == 0 || is_zero_terms(terms)) return polynomial();

    polynomial result;
    result.terms.clear();
    result.terms.reserve(terms.size());
    for (const auto& t : terms) {
        long long v = static_cast<long long>(t.second) * value;
        if (v != 0) result.terms.push_back({t.first, static_cast<coeff>(v)});
    }
    if (result.terms.empty()) result.terms.push_back({0, 0});
    return result;
}

polynomial operator*(int value, const polynomial& poly)
{
    return poly * value;
}

polynomial polynomial::operator*(const polynomial& other) const
{
    if (is_zero_terms(terms) || is_zero_terms(other.terms)) return polynomial();

    const size_t a_terms = terms.size();
    const size_t b_terms = other.terms.size();
    const size_t a_size = terms[0].first + 1;
    const size_t b_size = other.terms[0].first + 1;
    const size_t result_size = a_size + b_size - 1;
    const unsigned long long total_pairs = static_cast<unsigned long long>(a_terms) * b_terms;

    const bool dense_a = (a_terms * 4 >= a_size);
    const bool dense_b = (b_terms * 4 >= b_size);

    polynomial product;
    product.terms.clear();

    if (dense_a && dense_b && result_size >= NTT_THRESHOLD) {
        product.terms = ntt_multiply_dense(terms, other.terms, result_size);
        if (product.terms.empty()) product.terms.push_back({0, 0});
        return product;
    }

    if (total_pairs < SPARSE_THREAD_THRESHOLD) {
        std::map<power, long long> acc;
        for (const auto& a : terms) {
            for (const auto& b : other.terms) {
                acc[a.first + b.first] += static_cast<long long>(a.second) * b.second;
            }
        }
        for (auto it = acc.rbegin(); it != acc.rend(); ++it) {
            if (it->second != 0) product.terms.push_back({it->first, static_cast<coeff>(it->second)});
        }
        if (product.terms.empty()) product.terms.push_back({0, 0});
        return product;
    }

    const auto* outer = &terms;
    const auto* inner = &other.terms;
    if (outer->size() < inner->size()) std::swap(outer, inner);

    int nt = detect_thread_count();
    if (static_cast<size_t>(nt) > outer->size()) nt = static_cast<int>(outer->size());
    if (nt < 1) nt = 1;

    std::vector<pthread_t> tids(nt);
    std::vector<SparseMultData> jobs(nt);
    std::vector<std::unordered_map<power, long long>> locals(nt);

    size_t base = outer->size() / nt;
    size_t rem = outer->size() % nt;
    size_t start = 0;
    for (int t = 0; t < nt; ++t) {
        size_t chunk = base + (static_cast<size_t>(t) < rem ? 1 : 0);
        jobs[t] = {outer, inner, start, start + chunk, &locals[t]};
        pthread_create(&tids[t], nullptr, sparse_mult_worker, &jobs[t]);
        start += chunk;
    }
    for (int t = 0; t < nt; ++t) pthread_join(tids[t], nullptr);

    std::map<power, long long> final_acc;
    for (auto& local : locals) {
        for (const auto& kv : local) final_acc[kv.first] += kv.second;
    }

    for (auto it = final_acc.rbegin(); it != final_acc.rend(); ++it) {
        if (it->second != 0) product.terms.push_back({it->first, static_cast<coeff>(it->second)});
    }
    if (product.terms.empty()) product.terms.push_back({0, 0});
    return product;
}

polynomial polynomial::operator%(const polynomial& divisor) const
{
    if (is_zero_terms(divisor.terms)) {
        throw std::invalid_argument("Modulo by zero polynomial");
    }
    if (is_zero_terms(terms)) return polynomial();

    std::vector<std::pair<power, coeff>> r = terms;
    const std::vector<std::pair<power, coeff>>& d = divisor.terms;
    const power d_deg = d[0].first;
    const coeff d_lead = d[0].second;

    while (!r.empty() && !(r.size() == 1 && r[0].second == 0) && r[0].first >= d_deg) {
        if (r[0].second % d_lead != 0) break;
        coeff q_coeff = r[0].second / d_lead;
        if (q_coeff == 0) break;
        power q_pow = r[0].first - d_deg;

        std::vector<std::pair<power, coeff>> next;
        next.reserve(r.size() + d.size());
        size_t i = 0, j = 0;

        while (i < r.size() || j < d.size()) {
            power rp = (i < r.size()) ? r[i].first : static_cast<power>(0);
            power sp = (j < d.size()) ? d[j].first + q_pow : static_cast<power>(0);

            if (j >= d.size() || (i < r.size() && rp > sp)) {
                if (r[i].second != 0) next.push_back(r[i]);
                ++i;
            } else if (i >= r.size() || sp > rp) {
                long long val = -static_cast<long long>(d[j].second) * q_coeff;
                if (val != 0) next.push_back({sp, static_cast<coeff>(val)});
                ++j;
            } else {
                long long val = static_cast<long long>(r[i].second)
                              - static_cast<long long>(d[j].second) * q_coeff;
                if (val != 0) next.push_back({rp, static_cast<coeff>(val)});
                ++i;
                ++j;
            }
        }

        if (next.empty()) next.push_back({0, 0});
        r.swap(next);
    }

    polynomial result;
    result.terms = std::move(r);
    result.canonicalize();
    return result;
}
