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
#include <unistd.h>
#include <complex>

using power = size_t;
using coeff = int;

static constexpr double PI = 3.14159265358979323846;

static constexpr int MAX_THREADS = 8;
static constexpr size_t J_BLOCK = 4096;

static int detect_thread_count()
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n < 1) n = 1;
    if (n > MAX_THREADS) n = MAX_THREADS;
    return static_cast<int>(n);
}

// struct ThreadData {
//     const polynomial* subset_data;
//     const polynomial* full_data;
//     int start;
//     int end;
//     std::vector<std::pair<power, coeff>>* result;
// };
using cd = std::complex<double>;

static void fft(std::vector<cd>& a, bool invert)
{
    int n = (int)a.size();
    // bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2.0 * PI / len * (invert ? -1 : 1);
        cd wlen(std::cos(ang), std::sin(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i + j];
                cd v = a[i + j + len / 2] * w;
                a[i + j]         = u + v;
                a[i + j + len/2] = u - v;
                w *= wlen;
            }
        }
    }
    if (invert) {
        for (cd& x : a) x /= n;
    }
}

// Parallel FFT worker: computes butterfly passes startLen..endLen (exclusive)
struct FFTWork {
    cd*    a;
    int    n;
    int    start_len;   // first len (power of 2) this worker handles
    int    end_len;     // one past last len
    bool   invert;
};

static void* fft_worker(void* arg)
{
    FFTWork* w = static_cast<FFTWork*>(arg);
    cd* a = w->a;
    int n = w->n;
    for (int len = w->start_len; len <= w->end_len; len <<= 1) {
        double ang = 2.0 * PI / len * (w->invert ? -1 : 1);
        cd wlen(std::cos(ang), std::sin(ang));
        for (int i = 0; i < n; i += len) {
            cd pw(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i + j];
                cd v = a[i + j + len/2] * pw;
                a[i + j]         = u + v;
                a[i + j + len/2] = u - v;
                pw *= wlen;
            }
        }
    }
    return nullptr;
}

// ── FFT-based polynomial multiplication ─────────────────────────────────────
// Uses split-radix trick to reduce floating-point error:
// split each coefficient c into c_lo = c & mask, c_hi = c >> SPLIT
// so that intermediate values stay small enough for double precision.

static constexpr int SPLIT = 15;                // 2^15 = 32768
static constexpr int SPLIT_BASE = 1 << SPLIT;

// Returns canonical-form terms (descending power, no zeros); caller sets product.terms
static std::vector<std::pair<power,coeff>> fft_multiply_impl(
    const std::vector<std::pair<power,coeff>>& at,
    const std::vector<std::pair<power,coeff>>& bt,
    size_t result_size)
{
    // Round up to next power of 2
    size_t n = 1;
    while (n < result_size) n <<= 1;

    std::vector<cd> fa_lo(n), fa_hi(n), fb_lo(n), fb_hi(n);

    for (const auto& t : at) {
        int c = t.second;
        int lo = c & (SPLIT_BASE - 1);
        int hi = (c - lo) >> SPLIT; // arithmetic shift
        fa_lo[t.first] = lo;
        fa_hi[t.first] = hi;
    }
    for (const auto& t : bt) {
        int c = t.second;
        int lo = c & (SPLIT_BASE - 1);
        int hi = (c - lo) >> SPLIT;
        fb_lo[t.first] = lo;
        fb_hi[t.first] = hi;
    }

    // FFT the four arrays
    fft(fa_lo, false);
    fft(fa_hi, false);
    fft(fb_lo, false);
    fft(fb_hi, false);

    // result = a_lo*b_lo + (a_lo*b_hi + a_hi*b_lo)*BASE + a_hi*b_hi*BASE^2
    std::vector<cd> fres(n);
    for (size_t i = 0; i < n; i++) {
        cd ll = fa_lo[i] * fb_lo[i];
        cd lh = fa_lo[i] * fb_hi[i] + fa_hi[i] * fb_lo[i];
        cd hh = fa_hi[i] * fb_hi[i];
        fres[i] = ll + lh * (double)SPLIT_BASE + hh * (double)((long long)SPLIT_BASE * SPLIT_BASE);
    }
    fft(fres, true);

    std::vector<std::pair<power,coeff>> res;
    res.reserve(result_size);
    for (size_t k = result_size; k > 0; --k) {
        size_t idx = k - 1;
        long long v = llround(fres[idx].real());
        if (v != 0)
            res.push_back({idx, static_cast<coeff>(v)});
    }
    return res; // may be empty (caller interprets as zero poly)
}

void polynomial::canonicalize() {
    if (terms.empty()) {
        terms.push_back({0, 0});
        return;
    }

    std::sort(terms.begin(), terms.end(),
              [](const std::pair<power,coeff>& a, const std::pair<power,coeff>& b)
              { return a.first > b.first; });

    std::vector<std::pair<power,coeff>> combined;
    combined.reserve(terms.size());
    for (const auto& t : terms) {
        if (!combined.empty() && combined.back().first == t.first)
            combined.back().second += t.second;
        else
            combined.push_back(t);
    }

    std::vector<std::pair<power,coeff>> result;
    result.reserve(combined.size());
    for (const auto& t : combined)
        if (t.second != 0)
            result.push_back(t);

    terms = result.empty() ? std::vector<std::pair<power,coeff>>{{0,0}} : std::move(result);
}

polynomial::polynomial()
{
    terms.push_back({0,0});
}

polynomial::polynomial(const polynomial& other)
    : terms(other.terms) {}

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

// struct DenseMultData {
//     const int64_t *a_data;
//     const int64_t *b_data;
//     int64_t       *local_data;
//     size_t i_start;
//     size_t i_end;
//     size_t b_size;
// };

// static void *dense_mult_worker(void *arg) {
//     DenseMultData *d = static_cast<DenseMultData *>(arg);
//     const int64_t * __restrict__ a   = d->a_data;
//     const int64_t * __restrict__ b   = d->b_data;
//     int64_t       * __restrict__ out = d->local_data;
//     const size_t bs    = d->b_size;
//     const size_t i0    = d->i_start;
//     const size_t i1    = d->i_end;
//     for (size_t jb = 0; jb < bs; jb += J_BLOCK) {
//         const size_t je = (jb + J_BLOCK < bs) ? (jb + J_BLOCK) : bs;
//         for (size_t i = i0; i < i1; ++i) {
//             const int64_t av = a[i];
//             if (av == 0) continue;
//             int64_t * __restrict__ row = out + (i - i0);
//             for (size_t j = jb; j < je; ++j) {
//                 row[j] += av * b[j];
//             }
//         }
//     }
//     return nullptr;
// }

struct SparseMultData {
    const std::vector<std::pair<power,coeff>>* a;
    const std::vector<std::pair<power,coeff>>* b;
    size_t a_start;
    size_t a_end;
    std::unordered_map<power, int64_t>* local;
};

static void* sparse_mult_worker(void* arg)
{
    SparseMultData* d = static_cast<SparseMultData*>(arg);
    const auto& a = *d->a;
    const auto& b = *d->b;
    auto& local = *d->local;
    local.reserve((d->a_end - d->a_start) * b.size() + 16);
    for (size_t i = d->a_start; i < d->a_end; ++i) {
        int64_t ac = a[i].second;
        if (ac == 0) continue;
        power ap = a[i].first;
        for (const auto& bj : b) {
            if (bj.second == 0) continue;
            local[ap + bj.first] += ac * (int64_t)bj.second;
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
// polynomial polynomial::operator*(const polynomial &other) const
// {
//     bool this_zero  = (terms.size() == 1 && terms[0].second == 0);
//     bool other_zero = (other.terms.size() == 1 && other.terms[0].second == 0);
//     if (this_zero || other_zero)
//         return polynomial();
 
//     size_t a_terms = terms.size();
//     size_t b_terms = other.terms.size();
 
//     size_t a_deg = terms[0].first;
//     size_t b_deg = other.terms[0].first;
//     size_t a_size = a_deg + 1;
//     size_t b_size = b_deg + 1;
//     size_t result_size = a_size + b_size - 1;
 
//     bool dense_a   = (a_terms * 4 >= a_size);
//     bool dense_b   = (b_terms * 4 >= b_size);
//     bool result_ok = (result_size <= 2'000'000);
//     bool use_dense = dense_a && dense_b && result_ok;
 
//     polynomial product;
//     product.terms.clear();
 
//     if (use_dense)
//     {
//         // Build dense vectors.  Always make `a` the larger one so the
//         // outer (parallelizable) dimension has more work to spread.
//         const std::vector<std::pair<power, coeff>> *outer_terms = &terms;
//         const std::vector<std::pair<power, coeff>> *inner_terms = &other.terms;
//         size_t outer_size = a_size;
//         size_t inner_size = b_size;
//         if (a_size < b_size) {
//             std::swap(outer_terms, inner_terms);
//             std::swap(outer_size, inner_size);
//         }
 
//         std::vector<int64_t> A(outer_size, 0);
//         std::vector<int64_t> B(inner_size, 0);
//         for (const auto &t : *outer_terms) A[t.first] = t.second;
//         for (const auto &t : *inner_terms) B[t.first] = t.second;
//         std::vector<int64_t> result_dense(result_size, 0);
 
//         size_t total_work = outer_size * inner_size;
 
//         if (total_work < 50'000)
//         {
//             // Sequential -- thread overhead would dominate.
//             const int64_t * __restrict__ bp = B.data();
//             int64_t       * __restrict__ rp = result_dense.data();
//             for (size_t jb = 0; jb < inner_size; jb += J_BLOCK) {
//                 const size_t je = (jb + J_BLOCK < inner_size) ? (jb + J_BLOCK) : inner_size;
//                 for (size_t i = 0; i < outer_size; ++i) {
//                     const int64_t ai = A[i];
//                     if (ai == 0) continue;
//                     int64_t * __restrict__ row = rp + i;
//                     for (size_t j = jb; j < je; ++j)
//                         row[j] += ai * bp[j];
//                 }
//             }
//         }
//         else
//         {
//             int num_threads = detect_thread_count();
//             if (static_cast<size_t>(num_threads) > outer_size)
//                 num_threads = static_cast<int>(outer_size);
 
//             std::vector<pthread_t>      threads(num_threads);
//             std::vector<DenseMultData>  tdata(num_threads);
//             std::vector<std::vector<int64_t>> locals(num_threads);
 
//             size_t base  = outer_size / num_threads;
//             size_t rem   = outer_size % num_threads;
//             size_t start = 0;
//             for (int i = 0; i < num_threads; ++i)
//             {
//                 size_t chunk = base + (static_cast<size_t>(i) < rem ? 1 : 0);
//                 size_t local_len = chunk + inner_size - 1;
//                 locals[i].assign(local_len, 0);
//                 tdata[i] = {A.data(), B.data(), locals[i].data(),
//                             start, start + chunk, inner_size};
//                 pthread_create(&threads[i], nullptr,
//                                dense_mult_worker, &tdata[i]);
//                 start += chunk;
//             }
//             for (int i = 0; i < num_threads; ++i)
//                 pthread_join(threads[i], nullptr);
 
//             // Merge: each thread's local buffer covers global indices
//             // [i_start, i_start + local.size()).  Adjacent threads overlap
//             // by (inner_size - 1) entries -- the addition handles that.
//             for (int t = 0; t < num_threads; ++t)
//             {
//                 const int64_t *lp = locals[t].data();
//                 size_t off = tdata[t].i_start;
//                 size_t len = locals[t].size();
//                 int64_t *rp = result_dense.data() + off;
//                 for (size_t k = 0; k < len; ++k)
//                     rp[k] += lp[k];
//             }
//         }
 
//         product.terms.reserve(result_size);
//         for (size_t k = result_size; k > 0; --k)
//         {
//             size_t idx = k - 1;
//             int64_t v = result_dense[idx];
//             if (v != 0)
//                 product.terms.push_back({idx, static_cast<coeff>(v)});
//         }
//     }
//     else
//     {
//         size_t total_work = a_terms * b_terms;
 
//         if (total_work < 5'000)
//         {
//             std::map<power, int64_t> acc;
//             for (const auto &ai : terms)
//             {
//                 if (ai.second == 0) continue;
//                 int64_t ac = ai.second;
//                 for (const auto &bj : other.terms)
//                 {
//                     if (bj.second == 0) continue;
//                     acc[ai.first + bj.first] += ac * static_cast<int64_t>(bj.second);
//                 }
//             }
//             for (auto it = acc.rbegin(); it != acc.rend(); ++it)
//                 if (it->second != 0)
//                     product.terms.push_back({it->first, static_cast<coeff>(it->second)});
//         }
//         else
//         {
//             int num_threads = detect_thread_count();
//             if (static_cast<size_t>(num_threads) > a_terms)
//                 num_threads = static_cast<int>(a_terms);
 
//             std::vector<pthread_t>     threads(num_threads);
//             std::vector<SparseMultData> tdata(num_threads);
//             std::vector<std::unordered_map<power, int64_t>> locals(num_threads);
 
//             size_t base  = a_terms / num_threads;
//             size_t rem   = a_terms % num_threads;
//             size_t start = 0;
//             for (int i = 0; i < num_threads; ++i)
//             {
//                 size_t chunk = base + (static_cast<size_t>(i) < rem ? 1 : 0);
//                 tdata[i] = {&terms, &other.terms,
//                             start, start + chunk,
//                             &locals[i]};
//                 pthread_create(&threads[i], nullptr,
//                                sparse_mult_worker, &tdata[i]);
//                 start += chunk;
//             }
//             for (int i = 0; i < num_threads; ++i)
//                 pthread_join(threads[i], nullptr);
 
//             std::map<power, int64_t> final_map;
//             for (auto &local : locals)
//                 for (const auto &kv : local)
//                     final_map[kv.first] += kv.second;
 
//             for (auto it = final_map.rbegin(); it != final_map.rend(); ++it)
//                 if (it->second != 0)
//                     product.terms.push_back({it->first, static_cast<coeff>(it->second)});
//         }
//     }
 
//     if (product.terms.empty())
//         product.terms.push_back({0, 0});
//     return product;
// }

// polynomial polynomial::operator*(int scalar) const
// {

//     //polynomial product;
//     if (scalar == 0) {
//         return polynomial();
//     }

//     polynomial result;
//     result.terms.clear();
//     result.terms.reserve(terms.size());
//     for (const auto &i : this->terms)
//     {
//         int64_t c = static_cast<int64_t>(i.second) * scalar;
//         if (c != 0) {
//             result.terms.push_back({i.first, static_cast<coeff>(c)});
//         }
//         //product.terms.push_back({i.first, scalar * i.second});
//     }
//     if (result.terms.empty()) {
//         result.terms.push_back({0, 0});
//     }
//     return result;
// }

// polynomial operator*(int scalar, const polynomial &other)
// {
//     return other * scalar;
// }


// /*polynomial polynomial::operator+(const polynomial &other) const {
//     polynomial result;
//     result.terms.clear();
//     result.terms.insert(result.terms.end(), terms.begin(), terms.end());
//     result.terms.insert(result.terms.end(), other.terms.begin(), other.terms.end());

//     return result;
// }*/
// polynomial polynomial::operator+(const polynomial &other) const
// {
//     polynomial result;
//     result.terms.clear();
//     result.terms.reserve(terms.size() + other.terms.size());

//     auto it1 = terms.begin();
//     auto it2 = other.terms.begin();

//     while (it1 != terms.end() && it2 != other.terms.end())
//     {
//         if (it1->second == 0) { ++it1; continue; }
//         if (it2->second == 0) { ++it2; continue; }

//         if (it1->first > it2->first)
//         {
//             result.terms.push_back(*it1);
//             ++it1;
//         }
//         else if (it1->first < it2->first)
//         {
//             result.terms.push_back(*it2);
//             ++it2;
//         }
//         else
//         {
//             int sum = it1->second + it2->second;
//             if (sum != 0)
//                 result.terms.push_back({it1->first, sum});
//             ++it1; ++it2;
//         }
//     }
//     while (it1 != terms.end())       { if (it1->second != 0) result.terms.push_back(*it1); ++it1; }
//     while (it2 != other.terms.end()) { if (it2->second != 0) result.terms.push_back(*it2); ++it2; }

//     if (result.terms.empty())
//         result.terms.push_back({0, 0});
//     return result;
// }

// polynomial polynomial::operator+(int value) const {
// //     polynomial result = *this;
// //     result.terms.push_back({0, value});
// //     return result;
//     if (value == 0)
//         return *this;
//     polynomial constant;
//     constant.terms.clear();
//     constant.terms.push_back({0, value});
//     return *this + constant;
// }

// polynomial operator+(int value, const polynomial &poly) {
//     return poly + value;
// }

// /*polynomial polynomial::operator%(const polynomial &divisor) const {
//     auto r = this->canonical_form();
//     auto d = divisor.canonical_form();

//     if(d.size() == 1 && d[0].second == 0) {
//         throw std::invalid_argument("Modulo by a zero polynomial");
//     }

//     while (!r.empty() && r[0].first >= d[0].first) {
//         power power_diff = r[0].first - d[0].first;
//         coeff coeff_ratio = r[0].second / d[0].second;

//         std::vector<std::pair<power, coeff>> temp;

//         for (const auto &term : d) {
//             temp.push_back({term.first + power_diff, term.second * coeff_ratio});
//         }
//         std::map<power, coeff> acc;

//         for (const auto &term : r) {
//             acc[term.first] += term.second;
//         }

//         for (const auto &term : d)
//         {
//             acc[term.first + power_diff] -= term.second * coeff_ratio;
//         }

//         r.clear();

//         for (const auto &p : acc) {
//             if (p.second != 0) {
//                 r.push_back({p.first, p.second});
//             }
//         }

//         std::sort(r.begin(), r.end(), [](auto &a, auto &b)
//         {
//             return a.first > b.first;
//         });
//     }
//     if (r.empty()) {
//         return polynomial();
//     }
//     polynomial result;
//     result.terms = r;
//     return result;
// }*/

// polynomial polynomial::operator%(const polynomial &divisor) const
// {
//     if (divisor.terms.size() == 1 && divisor.terms[0].second == 0)
//         throw std::invalid_argument("Modulo by zero polynomial");

//     std::vector<std::pair<power, coeff>> r = terms;
//     const auto &d = divisor.terms;

//     if (r.size() == 1 && r[0].second == 0)
//         return polynomial();

//     power d_deg  = d[0].first;
//     coeff d_lead = d[0].second;

//     while (!r.empty() && r[0].first >= d_deg)
//     {
//         if (r[0].second == 0) { r.erase(r.begin()); continue; }
//         if (r[0].second % d_lead != 0) break;

//         coeff coeff_ratio = r[0].second / d_lead;
//         if (coeff_ratio == 0) break;
//         power power_diff = r[0].first - d_deg;

//         std::vector<std::pair<power, coeff>> new_r;
//         new_r.reserve(r.size() + d.size());

//         size_t i = 0, j = 0;
//         while (i < r.size() && j < d.size())
//         {
//             power r_p = r[i].first;
//             power d_p = d[j].first + power_diff;
//             if (r_p > d_p)
//             {
//                 if (r[i].second != 0) new_r.push_back(r[i]);
//                 ++i;
//             }
//             else if (r_p < d_p)
//             {
//                 coeff val = -d[j].second * coeff_ratio;
//                 if (val != 0) new_r.push_back({d_p, val});
//                 ++j;
//             }
//             else
//             {
//                 coeff val = r[i].second - d[j].second * coeff_ratio;
//                 if (val != 0) new_r.push_back({r_p, val});
//                 ++i; ++j;
//             }
//         }
//         while (i < r.size())
//         {
//             if (r[i].second != 0) new_r.push_back(r[i]);
//             ++i;
//         }
//         while (j < d.size())
//         {
//             coeff val = -d[j].second * coeff_ratio;
//             if (val != 0) new_r.push_back({d[j].first + power_diff, val});
//             ++j;
//         }
//         r = std::move(new_r);
//     }

//     if (r.empty())
//         return polynomial();
//     polynomial result;
//     result.terms = std::move(r);
//     return result;
// }
polynomial polynomial::operator*(const polynomial& other) const
{
    bool this_zero  = (terms.size() == 1 && terms[0].second == 0);
    bool other_zero = (other.terms.size() == 1 && other.terms[0].second == 0);
    if (this_zero || other_zero) return polynomial();

    size_t a_sz = terms.size();
    size_t b_sz = other.terms.size();
    size_t a_deg = terms[0].first;
    size_t b_deg = other.terms[0].first;
    size_t a_size = a_deg + 1;    // dense array size needed for A
    size_t b_size = b_deg + 1;
    size_t result_size = a_size + b_size - 1;

    bool dense_a = (a_sz * 4 >= a_size);
    bool dense_b = (b_sz * 4 >= b_size);

    polynomial product;
    product.terms.clear();

    // ── FFT path: large dense polynomials ────────────────────────────────────
    // FFT is O(n log n), much faster than O(n^2) for large n.
    // if (dense_a && dense_b && result_size > 4096) {
    //     auto res = fft_multiply_impl(terms, other.terms, result_size);
    //     polynomial prod;
    //     if (res.empty()) return prod; // zero
    //     prod.terms = std::move(res);
    //     return prod;
    // }

    // ── Dense O(n^2) path: small/medium dense polynomials ────────────────────
    if (dense_a && dense_b && result_size <= 2'000'000) {
        // Ensure outer dimension is larger for better parallelism
        const std::vector<std::pair<power,coeff>>* outer_terms = &terms;
        const std::vector<std::pair<power,coeff>>* inner_terms = &other.terms;
        size_t outer_size = a_size;
        size_t inner_size = b_size;
        if (a_size < b_size) {
            std::swap(outer_terms, inner_terms);
            std::swap(outer_size, inner_size);
        }

        std::vector<int64_t> A(outer_size, 0), B(inner_size, 0);
        for (const auto& t : *outer_terms) A[t.first] = t.second;
        for (const auto& t : *inner_terms) B[t.first] = t.second;
        std::vector<int64_t> result_dense(result_size, 0);

        // Sequential for small work, threaded for large
        size_t total_work = outer_size * inner_size;
        if (total_work < 50'000) {
            for (size_t i = 0; i < outer_size; ++i) {
                if (A[i] == 0) continue;
                for (size_t j = 0; j < inner_size; ++j)
                    result_dense[i + j] += A[i] * B[j];
            }
        } else {
            // Simple threaded approach: each thread handles a slice of i
            int nt = detect_thread_count();
            if ((size_t)nt > outer_size) nt = (int)outer_size;

            struct DenseSlice { const int64_t* A; const int64_t* B; int64_t* R; size_t i0, i1, bsz; };
            auto worker = [](void* arg) -> void* {
                DenseSlice* s = (DenseSlice*)arg;
                for (size_t i = s->i0; i < s->i1; ++i) {
                    if (s->A[i] == 0) continue;
                    for (size_t j = 0; j < s->bsz; ++j)
                        s->R[i + j] += s->A[i] * s->B[j];
                }
                return nullptr;
            };

            std::vector<pthread_t> tids(nt);
            std::vector<DenseSlice> slices(nt);
            // Each thread writes to a separate region of R (R[i0..i1+bsz-2]),
            // but adjacent regions overlap by (bsz-1). We handle this by
            // giving each thread its own scratch buffer and merging.
            std::vector<std::vector<int64_t>> scratch(nt);

            size_t base = outer_size / nt, rem = outer_size % nt, start = 0;
            for (int t = 0; t < nt; ++t) {
                size_t chunk = base + (t < (int)rem ? 1 : 0);
                size_t slen  = chunk + inner_size - 1;
                scratch[t].assign(slen, 0);
                slices[t] = {A.data(), B.data(), scratch[t].data(), 0, chunk, inner_size};
                // Worker sees local indices: A[i] means A[start+i], R[i+j]
                // Capture start via a lambda wrapper
                size_t gi = start;
                auto& sl = slices[t];
                sl = {A.data() + gi, B.data(), scratch[t].data(), 0, chunk, inner_size};
                pthread_create(&tids[t], nullptr, worker, &slices[t]);
                start += chunk;
            }
            for (int t = 0; t < nt; ++t) pthread_join(tids[t], nullptr);

            // Merge scratch buffers into result_dense
            start = 0;
            size_t base2 = outer_size / nt, rem2 = outer_size % nt;
            for (int t = 0; t < nt; ++t) {
                size_t chunk = base2 + (t < (int)rem2 ? 1 : 0);
                size_t slen  = scratch[t].size();
                for (size_t k = 0; k < slen; ++k)
                    result_dense[start + k] += scratch[t][k];
                start += chunk;
            }
        }

        product.terms.reserve(result_size);
        for (size_t k = result_size; k > 0; --k) {
            int64_t v = result_dense[k - 1];
            if (v != 0) product.terms.push_back({k - 1, (coeff)v});
        }
        if (product.terms.empty()) product.terms.push_back({0, 0});
        return product;
    }

    // ── Sparse path ──────────────────────────────────────────────────────────
    {
        size_t total_work = a_sz * b_sz;
        if (total_work < 5'000) {
            std::map<power, int64_t> acc;
            for (const auto& ai : terms) {
                if (ai.second == 0) continue;
                for (const auto& bj : other.terms) {
                    if (bj.second == 0) continue;
                    acc[ai.first + bj.first] += (int64_t)ai.second * bj.second;
                }
            }
            for (auto it = acc.rbegin(); it != acc.rend(); ++it)
                if (it->second != 0)
                    product.terms.push_back({it->first, (coeff)it->second});
        } else {
            int nt = detect_thread_count();
            if ((size_t)nt > a_sz) nt = (int)a_sz;

            std::vector<pthread_t> tids(nt);
            std::vector<SparseMultData> tdata(nt);
            std::vector<std::unordered_map<power, int64_t>> locals(nt);

            size_t base = a_sz / nt, rem = a_sz % nt, start = 0;
            for (int t = 0; t < nt; ++t) {
                size_t chunk = base + (t < (int)rem ? 1 : 0);
                tdata[t] = {&terms, &other.terms, start, start + chunk, &locals[t]};
                pthread_create(&tids[t], nullptr, sparse_mult_worker, &tdata[t]);
                start += chunk;
            }
            for (int t = 0; t < nt; ++t) pthread_join(tids[t], nullptr);

            std::map<power, int64_t> final_map;
            for (auto& lm : locals)
                for (const auto& kv : lm)
                    final_map[kv.first] += kv.second;

            for (auto it = final_map.rbegin(); it != final_map.rend(); ++it)
                if (it->second != 0)
                    product.terms.push_back({it->first, (coeff)it->second});
        }
    }

    if (product.terms.empty()) product.terms.push_back({0, 0});
    return product;
}

// ── scalar and free-function overloads ───────────────────────────────────────

polynomial polynomial::operator*(int scalar) const
{
    if (scalar == 0) return polynomial();
    polynomial result;
    result.terms.clear();
    result.terms.reserve(terms.size());
    for (const auto& t : terms) {
        int64_t c = (int64_t)t.second * scalar;
        if (c != 0) result.terms.push_back({t.first, (coeff)c});
    }
    if (result.terms.empty()) result.terms.push_back({0, 0});
    return result;
}

polynomial operator*(int scalar, const polynomial& other)
{
    return other * scalar;
}

polynomial polynomial::operator+(const polynomial& other) const
{
    polynomial result;
    result.terms.clear();
    result.terms.reserve(terms.size() + other.terms.size());

    auto it1 = terms.begin();
    auto it2 = other.terms.begin();
    while (it1 != terms.end() && it2 != other.terms.end()) {
        if (it1->second == 0) { ++it1; continue; }
        if (it2->second == 0) { ++it2; continue; }
        if (it1->first > it2->first) { result.terms.push_back(*it1++); }
        else if (it1->first < it2->first) { result.terms.push_back(*it2++); }
        else {
            int s = it1->second + it2->second;
            if (s != 0) result.terms.push_back({it1->first, s});
            ++it1; ++it2;
        }
    }
    while (it1 != terms.end())       { if (it1->second != 0) result.terms.push_back(*it1); ++it1; }
    while (it2 != other.terms.end()) { if (it2->second != 0) result.terms.push_back(*it2); ++it2; }

    if (result.terms.empty()) result.terms.push_back({0, 0});
    return result;
}

polynomial polynomial::operator+(int value) const
{
    if (value == 0) return *this;
    polynomial constant;
    constant.terms.clear();
    constant.terms.push_back({0, value});
    return *this + constant;
}

polynomial operator+(int value, const polynomial& poly)
{
    return poly + value;
}

// ── Modulo ───────────────────────────────────────────────────────────────────

polynomial polynomial::operator%(const polynomial& divisor) const
{
    if (divisor.terms.size() == 1 && divisor.terms[0].second == 0)
        throw std::invalid_argument("Modulo by zero polynomial");

    std::vector<std::pair<power,coeff>> r = terms;
    const auto& d = divisor.terms;

    if (r.size() == 1 && r[0].second == 0) return polynomial();

    power d_deg  = d[0].first;
    coeff d_lead = d[0].second;

    while (!r.empty() && r[0].first >= d_deg) {
        if (r[0].second == 0) { r.erase(r.begin()); continue; }
        if (r[0].second % d_lead != 0) break;

        coeff ratio = r[0].second / d_lead;
        power pdiff = r[0].first - d_deg;

        std::vector<std::pair<power,coeff>> new_r;
        new_r.reserve(r.size() + d.size());

        size_t i = 0, j = 0;
        while (i < r.size() && j < d.size()) {
            power rp = r[i].first;
            power dp = d[j].first + pdiff;
            if (rp > dp) {
                if (r[i].second != 0) new_r.push_back(r[i]);
                ++i;
            } else if (rp < dp) {
                coeff val = -d[j].second * ratio;
                if (val != 0) new_r.push_back({dp, val});
                ++j;
            } else {
                coeff val = r[i].second - d[j].second * ratio;
                if (val != 0) new_r.push_back({rp, val});
                ++i; ++j;
            }
        }
        while (i < r.size()) { if (r[i].second != 0) new_r.push_back(r[i]); ++i; }
        while (j < d.size()) {
            coeff val = -d[j].second * ratio;
            if (val != 0) new_r.push_back({d[j].first + pdiff, val});
            ++j;
        }
        r = std::move(new_r);
    }

    if (r.empty()) return polynomial();
    polynomial result;
    result.terms = std::move(r);
    return result;
}