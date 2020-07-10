#ifndef __CONTAINER__UTILS__
#define __CONTAINER__UTILS__

#include "lwe/randomness/prg.hpp"
#include <cassert>
#include <cstdint>
#include <stack>

namespace LWE {
    inline __int128_t mod_inv_proto(__int128_t A, __int128_t B, __int128_t S0,
                                    __int128_t S1) {
        if (B == 0 && A == 1)
            return S0;
        else if (B == 0 && A == 0)
            assert(0);
        return mod_inv_proto(B, A % B, S1, S0 - S1 * (A / B));
    }

    inline __int128_t modular_inverse(__int128_t A, __int128_t B) {
        __int128_t RES = mod_inv_proto(A, B, 1, 0);
        return RES + (RES < 0 ? B : 0);
    }

    constexpr size_t log2_ceil(__uint128_t n) {
        return ((n < 2) ? 1 : 1 + log2_ceil(n >> 1u));
    }

    constexpr uint64_t pow2(size_t n) {
        return ((n < 1) ? 1 : 2 * pow2(n - 1));
    }
}

namespace libsnark {
    std::ostream &operator<<(std::ostream &o, const __uint128_t &val) {
        __uint128_t my_var = val;
        std::stack<char> storing;
        auto stack_it = [&storing, &my_var](int flag) -> void {
            do {
                storing.push("0123456789abcdef"[my_var % flag]);
                my_var /= flag;
            } while (my_var != 0);
        };
        if (o.flags() & std::ostream::oct) {
            stack_it(8);
            o << "0";
        } else if (o.flags() & std::ostream::hex) {
            stack_it(16);
            o << "0x";
        } else if (o.flags() & std::ostream::dec)
            stack_it(10);
        else if (o.flags() & std::ostream::binary) {
            stack_it(2);
            o << "0b";
        }
        while (!storing.empty()) {
            o << storing.top();
            storing.pop();
        }
        return o;
    }

    std::istream &operator>>(std::istream &in, __uint128_t &val) {
        uint64_t top, bot;
        in >> top >> bot;
        __uint128_t _top = top, _bot = bot;
        val = (_top << 64) + _bot;
        return in;
    }
}

namespace Parse128T {
    constexpr uint8_t hex_to_digit(char c) noexcept {
        return c >= 'a' ? (10 + c - 'a') : c >= 'A' ? (10 + c - 'A') : c - '0';
    }

    template <int BASE, __uint128_t V>
    constexpr __uint128_t parse_128T() noexcept {
        return V;
    }

    template <int BASE, __uint128_t V, char C, char... Cs>
    constexpr __uint128_t parse_128T() noexcept {
        static_assert(BASE != 16 || sizeof...(Cs) <= 32 - 1,
                      "Literal too large for BASE=16");
        static_assert(BASE != 10 || sizeof...(Cs) <= 39 - 1,
                      "Literal too large for BASE=10");
        static_assert(BASE != 8 || sizeof...(Cs) <= 44 - 1,
                      "Literal too large for BASE=8");
        static_assert(BASE != 2 || sizeof...(Cs) <= 128 - 1,
                      "Literal too large for BASE=2");
        return parse_128T<BASE, BASE * V + hex_to_digit(C), Cs...>();
    }

    template <char... Cs> struct Parse128T {
        static constexpr __uint128_t eval() noexcept {
            return parse_128T<10, 0, Cs...>();
        }
    };

    template <char... Cs> struct Parse128T<'0', 'x', Cs...> {
        static constexpr __uint128_t eval() noexcept {
            return parse_128T<16, 0, Cs...>();
        }
    };

    template <char... Cs> struct Parse128T<'0', 'b', Cs...> {
        static constexpr __uint128_t eval() noexcept {
            return parse_128T<2, 0, Cs...>();
        }
    };

    template <char... Cs> struct Parse128T<'0', Cs...> {
        static constexpr __uint128_t eval() noexcept {
            return parse_128T<8, 0, Cs...>();
        }
    };

    template <char... Cs> constexpr __uint128_t operator"" _U128T() noexcept {
        return Parse128T<Cs...>::eval();
    }
}

template <char... Cs> constexpr __uint128_t operator"" _U128T() noexcept {
    return Parse128T::operator""_U128T<Cs...>();
}

#endif
