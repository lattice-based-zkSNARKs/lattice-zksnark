#ifndef __AES_H__
#define __AES_H__

#include <emmintrin.h>
#include <wmmintrin.h>
#include <xmmintrin.h>

#pragma GCC optimize("O3")

#define UNROLL_LOOP4

namespace LWERandomness {

    typedef __m128i block;

    typedef struct {
        block rd_key[15];
    } AES_KEY;

#define EXPAND_ASSIST(v1, v2, v3, v4, shuff_const, aes_const)                  \
    v2 = _mm_aeskeygenassist_si128(v4, aes_const);                             \
    v3 = _mm_castps_si128(                                                     \
        _mm_shuffle_ps(_mm_castsi128_ps(v3), _mm_castsi128_ps(v1), 16));       \
    v1 = _mm_xor_si128(v1, v3);                                                \
    v3 = _mm_castps_si128(                                                     \
        _mm_shuffle_ps(_mm_castsi128_ps(v3), _mm_castsi128_ps(v1), 140));      \
    v1 = _mm_xor_si128(v1, v3);                                                \
    v2 = _mm_shuffle_epi32(v2, shuff_const);                                   \
    v1 = _mm_xor_si128(v1, v2)

    static inline void AES_128_Key_Expansion(const unsigned char *userkey,
                                             AES_KEY *key) {
        __m128i x0, x1, x2;
        auto *kp = (__m128i *) key;
        kp[0] = x0 = _mm_loadu_si128((__m128i *) userkey);
        x2 = _mm_setzero_si128();
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 1);
        kp[1] = x0;
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 2);
        kp[2] = x0;
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 4);
        kp[3] = x0;
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 8);
        kp[4] = x0;
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 16);
        kp[5] = x0;
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 32);
        kp[6] = x0;
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 64);
        kp[7] = x0;
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 128);
        kp[8] = x0;
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 27);
        kp[9] = x0;
        EXPAND_ASSIST(x0, x1, x2, x0, 255, 54);
        kp[10] = x0;
    }

    static inline void AES_ecb_encrypt_blk(block *blk, const AES_KEY *key) {
        const __m128i *sched = ((__m128i *) (key->rd_key));
        *blk = _mm_xor_si128(*blk, sched[0]);
        *blk = _mm_aesenc_si128(*blk, sched[1]);
        *blk = _mm_aesenc_si128(*blk, sched[2]);
        *blk = _mm_aesenc_si128(*blk, sched[3]);
        *blk = _mm_aesenc_si128(*blk, sched[4]);
        *blk = _mm_aesenc_si128(*blk, sched[5]);
        *blk = _mm_aesenc_si128(*blk, sched[6]);
        *blk = _mm_aesenc_si128(*blk, sched[7]);
        *blk = _mm_aesenc_si128(*blk, sched[8]);
        *blk = _mm_aesenc_si128(*blk, sched[9]);
        *blk = _mm_aesenclast_si128(*blk, sched[10]);
    }

#ifndef UNROLL_LOOP4

    static inline void AES_ecb_encrypt_blks(block *blks, int nblks,
                                            const AES_KEY *key) {
        const __m128i *sched = ((__m128i *) (key->rd_key));
        for (int i = 0; i < nblks; ++i) {
            blks[i] = _mm_xor_si128(blks[i], sched[0]);
            blks[i] = _mm_aesenc_si128(blks[i], sched[1]);
            blks[i] = _mm_aesenc_si128(blks[i], sched[2]);
            blks[i] = _mm_aesenc_si128(blks[i], sched[3]);
            blks[i] = _mm_aesenc_si128(blks[i], sched[4]);
            blks[i] = _mm_aesenc_si128(blks[i], sched[5]);
            blks[i] = _mm_aesenc_si128(blks[i], sched[6]);
            blks[i] = _mm_aesenc_si128(blks[i], sched[7]);
            blks[i] = _mm_aesenc_si128(blks[i], sched[8]);
            blks[i] = _mm_aesenc_si128(blks[i], sched[9]);
            blks[i] = _mm_aesenclast_si128(blks[i], sched[10]);
        }
    }

#else

#define UNROLL_XOR(srt, index)                                                 \
    blks[srt] = _mm_xor_si128(blks[srt], sched[index]);                        \
    blks[srt + 1] = _mm_xor_si128(blks[srt + 1], sched[index]);                \
    blks[srt + 2] = _mm_xor_si128(blks[srt + 2], sched[index]);                \
    blks[srt + 3] = _mm_xor_si128(blks[srt + 3], sched[index])

#define UNROLL_CHUNK(srt, index)                                               \
    blks[srt] = _mm_aesenc_si128(blks[srt], sched[index]);                     \
    blks[srt + 1] = _mm_aesenc_si128(blks[srt + 1], sched[index]);             \
    blks[srt + 2] = _mm_aesenc_si128(blks[srt + 2], sched[index]);             \
    blks[srt + 3] = _mm_aesenc_si128(blks[srt + 3], sched[index])

#define UNROLL_LAST(srt, index)                                                \
    blks[srt] = _mm_aesenclast_si128(blks[srt], sched[index]);                 \
    blks[srt + 1] = _mm_aesenclast_si128(blks[srt + 1], sched[index]);         \
    blks[srt + 2] = _mm_aesenclast_si128(blks[srt + 2], sched[index]);         \
    blks[srt + 3] = _mm_aesenclast_si128(blks[srt + 3], sched[index])

    static inline void AES_ecb_encrypt_blks(block *blks, int nblks,
                                            const AES_KEY *key) {
        const __m128i *sched = ((__m128i *) (key->rd_key));
        int nblk_pipeline = nblks >> 2, nblk_remain = nblks & 3;
        for (int i = 0; i < nblk_pipeline; ++i) {
            int srt = i << 2;
            UNROLL_XOR(srt, 0);
            UNROLL_CHUNK(srt, 1);
            UNROLL_CHUNK(srt, 2);
            UNROLL_CHUNK(srt, 3);
            UNROLL_CHUNK(srt, 4);
            UNROLL_CHUNK(srt, 5);
            UNROLL_CHUNK(srt, 6);
            UNROLL_CHUNK(srt, 7);
            UNROLL_CHUNK(srt, 8);
            UNROLL_CHUNK(srt, 9);
            UNROLL_LAST(srt, 10);
        }
        for (int i = 0; i < nblk_remain; i++)
            AES_ecb_encrypt_blk(&blks[i + (nblk_pipeline << 2)], key);
    }

#endif

}

#endif
