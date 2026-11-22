// Bounded-width reverse meander transfer. Each parity writes to disjoint cells.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

using Word = std::uint32_t;
using Index = std::uint32_t;
constexpr Index missing = (1u << 24) - 1;

struct ThreeBytes {
    std::uint8_t a, b, c;
    void set(Index rank) { a = rank; b = rank >> 8; c = rank >> 16; }
    Index get() const { return Index(a) | (Index(b) << 8) | (Index(c) << 16); }
};
static_assert(sizeof(ThreeBytes) == 3);

struct Count {
    std::uint64_t low{}, middle{}, high{};
    void add(const Count& other) {
        std::uint64_t l, m, h;
        const bool c0 = __builtin_add_overflow(low, other.low, &l);
        const bool c1 = __builtin_add_overflow(middle, other.middle, &m);
        const bool c2 = __builtin_add_overflow(m, std::uint64_t(c0), &m);
        const bool c3 = __builtin_add_overflow(high, other.high, &h);
        const bool c4 = __builtin_add_overflow(h, std::uint64_t(c1 || c2), &h);
        if (c3 || c4) low = middle = high = ~std::uint64_t(0);
        else { low = l; middle = m; high = h; }
    }
    std::string decimal() const {
        if ((low & middle & high) == ~std::uint64_t(0)) return "SATURATED";
        auto v = *this;
        std::string result;
        while (v.low || v.middle || v.high) {
            unsigned __int128 p = v.high;
            v.high = p / 10; p = ((p % 10) << 64) | v.middle;
            v.middle = p / 10; p = ((p % 10) << 64) | v.low;
            v.low = p / 10; result.push_back(char('0' + p % 10));
        }
        std::reverse(result.begin(), result.end());
        return result.empty() ? "0" : result;
    }
};
static_assert(sizeof(Count) == 24);

// Rank a Dyck word in O(1). Its first half determines the initial offset;
// the last half has a unique possible start height if it ends on the line.
struct Layer {
    int half{};
    Word mask{};
    Index base{};
    std::vector<Index> first, last;
    std::vector<Word> words;
    std::vector<Index> reflection;
    std::vector<ThreeBytes> down;
    Index rank(Word p) const { return first[p & mask] + last[p >> half]; }
    Index cell(Index r, int u) const { return base + r * (half + 1) + u; }
    void generate(int k, Index offset) {
        half = k; base = offset; mask = (1u << half) - 1;
        first.resize(mask + 1); last.resize(mask + 1);
        std::vector<std::vector<Word>> suffixes(k + 1);
        for (Word bits = 0; bits <= mask; ++bits) {
            int net = 0, minimum = 0;
            for (int i = 0; i < k; ++i) {
                net += (bits & (1u << i)) ? 1 : -1;
                minimum = std::min(minimum, net);
            }
            if (-net >= 0 && minimum >= net) {
                auto& valid = suffixes[-net];
                last[bits] = valid.size();
                valid.push_back(bits);
            }
        }
        Index count = 0;
        for (Word bits = 0; bits <= mask; ++bits) {
            int level = 0;
            for (int i = 0; i < k && level >= 0; ++i)
                level += (bits & (1u << i)) ? 1 : -1;
            if (level >= 0) { first[bits] = count; count += suffixes[level].size(); }
        }
        words.reserve(count);
        for (Word bits = 0; bits <= mask; ++bits) {
            int level = 0;
            for (int i = 0; i < k && level >= 0; ++i)
                level += (bits & (1u << i)) ? 1 : -1;
            if (level >= 0) for (auto suffix : suffixes[level])
                words.push_back(bits | (suffix << k));
        }
        reflection.resize(count);
        down.resize(std::size_t(count) * half);
    }
};

void make_edges(std::vector<Layer>& layers, int k) {
    auto& layer = layers[k];
    if (k > 1 && layers[k-1].words.size() >= missing)
        throw std::overflow_error("3-byte matching rank");
    #pragma omp parallel for schedule(static)
    for (Index r = 0; r < layer.words.size(); ++r) {
        const Word p = layer.words[r];
        int mate[32], stack[16], top = 0;
        Word reflected = 0;
        for (int i = 0; i < 2*k; ++i) {
            if (p & (1u << i)) stack[top++] = i;
            else {
                const int j = stack[--top];
                mate[i] = j; mate[j] = i;
                reflected |= 1u << (2*k-1-i);
            }
        }
        layer.reflection[r] = layer.rank(reflected);
        for (int u = 1; u <= k; ++u) {
            const int a = mate[u-1], b = mate[u];
            Index target = missing;
            if (a != u) {
                const Word edited = (p | (1u << std::min(a,b))) & ~(1u << std::max(a,b));
                const Word mask = (1u << (u-1)) - 1;
                target = layers[k-1].rank((edited & mask) | ((edited >> (u+1)) << (u-1)));
            }
            layer.down[std::size_t(r)*k + u-1].set(target);
        }
    }
}

void sweep_layer(const std::vector<Layer>& layers, int k, int cap,
                 int parity, std::vector<Count>& weights) {
    const auto& layer = layers[k];
    const auto& smaller = layers[k-1];
    const auto& larger = layers[std::min(k+1,cap)];
    #pragma omp parallel for schedule(static)
    for (Index r = 0; r < layer.words.size(); ++r) {
        const Word p = layer.words[r];
        const Index start = layer.cell(r,0);
        for (int u = parity; u <= k; u += 2) {
            Count z{};
            if (u) z = weights[start+u-1];
            if (u < k) z.add(weights[start+u+1]);
            else z.add(weights[layer.cell(layer.reflection[r],u-1)]);
            if (k < cap) {
                const Word mask = (1u << u)-1;
                const Word opened = (p & mask) | (1u << u) | ((p >> u) << (u+2));
                z.add(weights[larger.cell(larger.rank(opened),u+1)]);
            }
            if (u) {
                const Index d = layer.down[std::size_t(r)*k+u-1].get();
                if (d != missing) z.add(weights[smaller.cell(d,u-1)]);
            }
            weights[start+u] = z;
        }
    }
}

int main(int argc, char** argv) {
    if (argc != 3) return 1;
    const int width = std::stoi(argv[1]), orders = std::stoi(argv[2]);
    if (width < 2 || width % 2 || width > 32 || orders < 2) return 2;
    const int cap = width / 2;
    const auto started = std::chrono::steady_clock::now();
    auto seconds = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count(); };
    std::vector<Layer> layers(cap+1);
    std::uint64_t total = 0;
    for (int k=1; k<=cap; ++k) {
        layers[k].generate(k, total);
        total += layers[k].words.size() * (k+1);
        if (total > std::numeric_limits<Index>::max()) throw std::overflow_error("state index");
    }
    for (int k=1; k<=cap; ++k) make_edges(layers,k);
    std::vector<Count> weights(total);
    const Index root = layers[1].cell(0,1), shifted = layers[1].cell(0,0);
    weights[root].low = 1;
    std::cerr << "initialized " << seconds() << " seconds\n";
    std::cout << "# width " << width << " orders " << orders << " states " << total << std::endl;
    for (int dist = 0; dist <= 2*orders-2; ++dist) {
        if (dist % 2) std::cout << "B " << (dist+3)/2 << " " << weights[shifted].decimal() << std::endl;
        else std::cout << "C " << (dist+2)/2 << " " << weights[root].decimal() << std::endl;
        if (dist % 10 == 0) std::cerr << "distance " << dist << " elapsed " << seconds() << "\n";
        if (dist == 2*orders-2) break;
        // Opposite parities occupy different cells. Higher layers are still zero.
        for (int k = 1; k <= std::min(cap,dist+2); ++k)
            sweep_layer(layers,k,cap,dist % 2,weights);
    }
}
