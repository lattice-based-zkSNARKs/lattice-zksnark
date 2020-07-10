#ifndef __PSEUDO_RANDOM_GENERATOR_HPP__
#define __PSEUDO_RANDOM_GENERATOR_HPP__

#include "aes.hpp"
#include "lwe/lwe_params.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

namespace LWERandomness {

    /* Randomness is derived from internal PRG.
       Implemented with AES in CTR mode.
       BEWARE! THIS IS NOT THREAD SAFE! */

    using byte = unsigned char;

    const uint32_t AES_KEY_BYTES = 16;

    constexpr const double Pi = 3.141592653589793238462L;

    class PseudoRandomGenerator {
    private:
        __uint128_t _counter;
        AES_KEY PRG_key;

    public:
        PseudoRandomGenerator() : _counter(0), PRG_key{} {
            static std::ifstream urandom("/dev/urandom", std::ios::binary);
            byte buffer[AES_KEY_BYTES];
            urandom.read(reinterpret_cast<char *>(buffer), AES_KEY_BYTES);
            urandom.close();
            AES_128_Key_Expansion(buffer, &this->PRG_key);
        }

        PseudoRandomGenerator(PseudoRandomGenerator &&prg) = default;

        explicit PseudoRandomGenerator(const AES_KEY &_key)
            : _counter(0), PRG_key{} {
            for (int i = 0; i < 15; i++)
                this->PRG_key.rd_key[i] = _key.rd_key[i];
        }

        __uint128_t next_prg_block() {
            auto blk = block(this->_counter++);
            AES_ecb_encrypt_blk(&blk, &this->PRG_key);
            return __uint128_t(blk);
        }

        template <uint64_t LENGTH, typename T>
        inline void prg_mem_randomize(std::array<T, LENGTH> &_dest) {
            size_t T_size = sizeof(T);
            uint64_t full_block_num = T_size * LENGTH / AES_KEY_BYTES,
                     remain_size = T_size * LENGTH % AES_KEY_BYTES;
            auto _block_dest = reinterpret_cast<block *>(&_dest);
            for (uint64_t i = 0; i < full_block_num; i++) {
                memcpy(_block_dest + i, &this->_counter, sizeof(block));
                this->_counter++;
            }
            AES_ecb_encrypt_blks(_block_dest, full_block_num, &this->PRG_key);
            if (remain_size != 0) {
                auto remainder = this->next_prg_block();
                memcpy(_block_dest + full_block_num, &remainder, remain_size);
            }
        }

        __uint128_t bounded(__uint128_t bound) {
            if (!(bound & (bound - 1)))
                return this->next_prg_block() & (bound - 1);
            else {
                // REJECT SAMPLING
                __uint128_t w_param = (__uint128_t(0) - 1) / bound,
                            w_bound = w_param * bound,
                            buffer = this->next_prg_block();
                while (buffer >= w_bound)
                    buffer = this->next_prg_block();
                return buffer % bound;
            }
        }

        __int128_t pm_bounded(__uint128_t bound) {
            // bound should be less than 2^127
            // REJECT SAMPLING
            __uint128_t w_param = (~(__uint128_t(1) << 127)) / bound;
            __int128_t w_bound = w_param * bound,
                       buffer = __int128_t(this->next_prg_block());
            while (buffer >= w_bound || buffer <= -w_bound)
                buffer = __int128_t(this->next_prg_block());
            return buffer % bound;
        }
    };

    class DiscreteGaussian {
    public:
        double width;
        double expand;
        std::map<uint64_t, int> probability_interval;
        PseudoRandomGenerator &prg;

        DiscreteGaussian(const double _width, const double _expand,
                         PseudoRandomGenerator &_prg)
            : width(_width), expand(_expand), prg(_prg) {
            int boundary = (int) ceil(this->width * this->expand);
            std::map<int, double> temp;
            double normalize = 0.0;
            uint64_t accumulate = 0;

            for (int i = 0; i <= boundary; i++) {
                double prob = std::exp(-Pi * i * i / std::pow(this->width, 2));
                normalize += prob * (i == 0 ? 1 : 2);
                temp[i] = prob;
                temp[-i] = prob;
            }
            for (const auto &i : temp) {
                this->probability_interval[accumulate] = i.first;
                uint64_t bucket = (i.second / normalize) * double(UINT64_MAX);
                accumulate += bucket;
            }
        }

        DiscreteGaussian(DiscreteGaussian &&dg) = default;

        int sample() {
            __uint128_t sample_res = prg.next_prg_block();
            auto rnd = sample_res >> 64u;
            auto which_bucket = --this->probability_interval.lower_bound(rnd);
            return which_bucket->second;
        }
    };
}

#endif
