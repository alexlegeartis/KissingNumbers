// Enumerate a prescribed residue class in an exact enclosing LDL ellipsoid.
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

using Integer = __int128_t;
using Unsigned = __uint128_t;
constexpr int dimension = 40;
using Row = std::array<int64_t, dimension>;
using Matrix = std::array<Row, dimension>;

Integer parse_radius(const std::string& text) {
    Integer radius = 0;
    for (char ch : text) {
        assert(ch >= '0' && ch <= '9' && radius < (Integer{1} << 116));
        radius = 10 * radius + (ch - '0');
    }
    assert(radius > 0 && radius < (Integer{1} << 120));
    return radius;
}

int64_t floor_divide(int64_t x, int64_t y) {
    return x / y - (x % y < 0);
}

uint64_t integer_sqrt(Unsigned n) {
    uint64_t root = static_cast<uint64_t>(std::sqrt(static_cast<long double>(n)));
    while (Unsigned{root} * root > n) --root;
    while (Unsigned{root + 1} * (root + 1) <= n) ++root;
    return root;
}

class Enumeration {
    Matrix gram{}, lower{};
    Row diagonal{}, residue{}, point{};
    int64_t scale = 0;
    Integer radius = 0;
    std::ostream& output;
    uint64_t count = 0;

    void visit(int k, Integer remaining) {
        if (k < 0) {
            Integer norm = 0;
            for (int i = 0; i < dimension; ++i)
                for (int j = 0; j < dimension; ++j)
                    norm += Integer{point[i]} * gram[i][j] * point[j];
            if (norm == 54) {
                std::array<int8_t, dimension> bytes{};
                for (int i = 0; i < dimension; ++i) {
                    assert(point[i] >= -127 && point[i] <= 127);
                    bytes[i] = static_cast<int8_t>(point[i]);
                }
                output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
                ++count;
            }
            return;
        }
        Integer translation = 0;
        for (int j = k + 1; j < dimension; ++j)
            translation += Integer{lower[j][k]} * point[j];
        assert(translation > -(Integer{1} << 60) && translation < (Integer{1} << 60));
        assert(remaining >= 0 && remaining < (Integer{1} << 120));
        const int64_t shift = static_cast<int64_t>(translation);
        const int64_t limit = integer_sqrt(static_cast<Unsigned>(remaining / diagonal[k]));
        assert(limit < (int64_t{1} << 60));
        int64_t lo = -floor_divide(limit + shift, scale);
        const int64_t hi = floor_divide(limit - shift, scale);
        lo += ((residue[k] - lo) % 3 + 3) % 3;
        assert(lo > -(int64_t{1} << 20) && hi < (int64_t{1} << 20));
        for (int64_t value = lo; value <= hi; value += 3) {
            point[k] = value;
            const Integer offset = Integer{scale} * value + shift;
            assert(offset > -(Integer{1} << 44) && offset < (Integer{1} << 44));
            const Integer cost = Integer{diagonal[k]} * offset * offset;
            assert(cost >= 0 && cost <= remaining);
            visit(k - 1, remaining - cost);
        }
    }

public:
    Enumeration(std::istream& input, std::ostream& output) : output(output) {
        std::string text;
        input >> scale >> text;
        assert(scale > 0 && scale < (int64_t{1} << 40));
        radius = parse_radius(text);
        for (int64_t& value : diagonal) {
            input >> value;
            assert(value > 0 && value < (int64_t{1} << 36));
        }
        for (Row& row : lower)
            for (int64_t& value : row) {
                input >> value;
                assert(value > -(int64_t{1} << 40) && value < (int64_t{1} << 40));
            }
        for (Row& row : gram)
            for (int64_t& value : row) {
                input >> value;
                assert(value > -(int64_t{1} << 20) && value < (int64_t{1} << 20));
            }
        for (int64_t& value : residue) input >> value;
        assert(input);
    }

    uint64_t run() {
        visit(dimension - 1, radius);
        assert(output);
        return count;
    }
};

int main(int argc, char** argv) {
    if (argc != 3) return 1;
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2], std::ios::binary);
    if (!input || !output) return 2;
    std::cout << Enumeration(input, output).run() << '\n';
}
