#include <iostream>
#include <chrono>
#include <optional>
#include <vector>

#include "poly.h"

std::optional<double> poly_test(polynomial& p1,
                                polynomial& p2,
                                std::vector<std::pair<power, coeff>> solution)

{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    polynomial p3 = p1 * p2;

    auto p3_can_form = p3.canonical_form();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    p3.print();

    if (p3_can_form != solution)
    {
        return std::nullopt;
    }

    return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
}

std::vector<std::pair<power, coeff>> make_poly(int n)
{
    std::vector<std::pair<power, coeff>> poly;

    for (int i = 0; i < n; i++)
    {
        poly.push_back({i, 1});  // x^i with coefficient 1
    }

    return poly;
}

int main()
{
    // /** We're doing (x+1)^2, so solution is x^2 + 2x + 1*/
    // std::vector<std::pair<power, coeff>> solution = {{2,1}, {1,2}, {0,1}};

    // /** This holds (x+1), which we'll pass to each polynomial */
    // std::vector<std::pair<power, coeff>> poly_input = {{1,1}, {0,1}};

    // polynomial p1(poly_input.begin(), poly_input.end());
    // polynomial p2(poly_input.begin(), poly_input.end());

    std::vector<std::pair<power, coeff>> poly_input = make_poly(100);
    std::vector<std::pair<power, coeff>> poly_input2 = make_poly(100);

    polynomial p1(poly_input.begin(), poly_input.end());
    polynomial p2(poly_input2.begin(), poly_input2.end());

    polynomial solution_poly = p1 * p2;
    auto solution = solution_poly.canonical_form();

    std::optional<double> result = poly_test(p1, p2, solution);

    if (result.has_value())
    {
        std::cout << "Passed test, took " << result.value()/1000 << " seconds" << std::endl;
    } 
    else 
    {
        std::cout << "Failed test" << std::endl;
    }
}