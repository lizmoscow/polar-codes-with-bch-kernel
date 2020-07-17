#include <math.h>
#include "SoftProcessing.h"
#include <crtdbg.h>
#include <cassert>
#if defined AVX
#include <immintrin.h>
#endif

#include <iostream>

// Saturate addition for integers.
#if defined(USE_FLOAT) || defined(USE_INT32)
MType adds(const MType &a, const MType &b) {
    return a + b;
}

MType subs(const MType &a, const MType &b) {
    return a - b;
}
#elif defined(USE_INT8)
// -(-128) == -128 problem
int8_t adds(const int8_t &a, const int8_t &b) {
    int16_t result = a + b;
    if (result >  127) result =  127;
    if (result < -127) result = -127;
    return result;
}

int8_t subs(const int8_t &a, const int8_t &b) {
    int16_t result = a - b;
    if (result >  127) result =  127;
    if (result < -127) result = -127;
    return result;
}
#endif

/// Calculation of G function: 
/// g(lambda, beta) = (-1)^C_l[beta][0]*M_{l-1}[2*beta] + M_{l-1}[2*beta+1]
inline void SoftCombine2Elements(MType *M_l, const MType &M_l_1, const MType &M_l_2, const tBit &C_l) 
{
	assert((C_l & 0x7f) == 0);
    bool u1 = (C_l != 0);
    if (u1) {
        *M_l = subs(M_l_2, M_l_1); SUM_COUNT;
    }
    else {
        *M_l = adds(M_l_2, M_l_1); SUM_COUNT;
    }
    assert((*M_l >= 0) || (*M_l < 0));
}

/// Calculation of Q function:
/// Q(a, b) = sign(a)*sign(b)*min(a,b)
inline void SoftXOR2Elements(MType *M_l, const MType &M_l_1, const MType &M_l_2) {
    MType mM_l_1 = umin(M_l_1), mM_l_2 = umin(M_l_2);
    CMP_COUNT;
    if (M_l_1 < M_l_2)
    {
        if (M_l_1 > mM_l_2)
        {
            *M_l = M_l_1;
        }
        else
        {
            *M_l = mM_l_2;
        };
    }
    else
    {
        if (M_l_1 > mM_l_2)
        {
            *M_l = M_l_2;
        }
        else
        {
            *M_l = mM_l_1;
        }
    }
    assert((*M_l >= 0) || (*M_l < 0));
}

#if defined(AVX) 
#if defined(USE_FLOAT)
const __m128i UnpackMask = _mm_set_epi8(3, -1, -1, -1, 2, -1, -1, -1, 1, -1, -1, -1, 0, -1, -1, -1);
void SoftCombine(float *M_l, const float *M_l_1, const float *M_l_2, const tBit *C_l, const uint32_t &N) 
{
	if (N >= 8) {
		SUM_BLOCK_COUNT(N);
#pragma unroll(4)
		for (unsigned beta = 0; beta < N; beta += 8) {
			__m128i C_1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(C_l + beta));
			__m128i C_2 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(C_l + beta + 4));
			//assert((C_l[beta] | C_l[beta + 1] | C_l[beta + 2] | C_l[beta + 3] | C_l[beta + 4] | C_l[beta + 5] | C_l[beta + 6] | C_l[beta + 7]) < 2);
			C_1 = _mm_shuffle_epi8(C_1, UnpackMask);
			C_2 = _mm_shuffle_epi8(C_2, UnpackMask);
			__m256 C = _mm256_set_m128(_mm_castsi128_ps(C_2), _mm_castsi128_ps(C_1));
			__m256 A = _mm256_load_ps(M_l_1 + beta); __m256 B = _mm256_load_ps(M_l_2 + beta);
			__m256 Sum = _mm256_add_ps(A, B); __m256 Diff = _mm256_sub_ps(B, A);
			__m256 M = _mm256_blendv_ps(Sum, Diff, C);
			_mm256_store_ps(M_l + beta, M);
            for (int i = 0; i < 8; ++i)
			    assert((M_l[beta + i] >= 0) || (M_l[beta + i] < 0));
		}
	}
	else if (N == 4)
	{
		SUM_BLOCK_COUNT(N);
		__m128i C = _mm_set1_epi32(*(reinterpret_cast<const int*>(C_l)));
		C = _mm_shuffle_epi8(C, UnpackMask);
		__m128 M = _mm_load_ps(M_l_1);
		M = _mm_xor_ps(_mm_castsi128_ps(C), M);
		M = _mm_add_ps(M, _mm_load_ps(M_l_2));
		_mm_store_ps(M_l, M);
        for (int i = 0; i < 4; ++i)
    		assert((M_l[i] >= 0) || (M_l[i] < 0));
	}
	else if (N == 2)
	{
		SUM_BLOCK_COUNT(N);
		__m128i C = _mm_set1_epi32(*(reinterpret_cast<const int*>(C_l)));
		C = _mm_shuffle_epi8(C, UnpackMask);
		__m128 M = _mm_load_ps(M_l_1);
		M = _mm_xor_ps(_mm_castsi128_ps(C), M);
		M = _mm_add_ps(M, _mm_loadu_ps(M_l_2));
		_mm_storel_pd(reinterpret_cast<double*>(M_l), _mm_castps_pd(M));
        for (int i = 0; i < 2; ++i)
    		assert((M_l[i] >= 0) || (M_l[i] < 0));
	}
	else 
	{
		for(unsigned i=0;i<N;i++)
			SoftCombine2Elements(&M_l[i], M_l_1[i], M_l_2[i], C_l[i]);
	}
}

extern const __m256 AbsMask = _mm256_castsi256_ps(_mm256_set_epi32(0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF));
extern const __m128 AbsMask128 = _mm_castsi128_ps(_mm_set_epi32(0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF));

float SoftCombineWithSum(float *M_l, const float *M_l_1, const float *M_l_2, const tBit *C_l, const uint32_t &N) {
	const __m128i UnpackMask = _mm_set_epi8(3, -1, -1, -1, 2, -1, -1, -1, 1, -1, -1, -1, 0, -1, -1, -1);
	
	if (N >= 8) {
		SUM_BLOCK_COUNT(2*N);
		__m256 Sum = _mm256_set1_ps(0);
#pragma unroll(4)
		for (unsigned beta = 0; beta < N; beta += 8) 
		{
			__m128i C_1 = _mm_castpd_si128(_mm_loadu_pd(reinterpret_cast<const double*>(C_l + beta)));
			__m128i C_2 = _mm_castpd_si128(_mm_loadu_pd(reinterpret_cast<const double*>(C_l + beta + 4)));
			C_1 = _mm_shuffle_epi8(C_1, UnpackMask);
			C_2 = _mm_shuffle_epi8(C_2, UnpackMask);
			__m256 C = _mm256_set_m128(_mm_castsi128_ps(C_2), _mm_castsi128_ps(C_1));
			__m256 M = _mm256_load_ps(M_l_1 + beta);
			M = _mm256_xor_ps(C, M);
			M = _mm256_add_ps(M, _mm256_load_ps(M_l_2 + beta));
			_mm256_store_ps(M_l + beta, M);
			M = _mm256_and_ps(M, AbsMask);
			Sum = _mm256_add_ps(Sum, M);
            for (int i = 0; i < 8; ++i)
    			assert((M_l[beta + i] >= 0) || (M_l[beta + i] < 0));
		}
		Sum = _mm256_hadd_ps(Sum, Sum);
		Sum = _mm256_hadd_ps(Sum, Sum);
		return Sum.m256_f32[0] + Sum.m256_f32[4];
	}
	else if (N == 4)
	{
		SUM_BLOCK_COUNT(N);
		__m128i C = _mm_set1_epi32(*(reinterpret_cast<const int*>(C_l)));
		C = _mm_shuffle_epi8(C, UnpackMask);
		__m128 M = _mm_load_ps(M_l_1);
		M = _mm_xor_ps(_mm_castsi128_ps(C), M);
		M = _mm_add_ps(M, _mm_load_ps(M_l_2));
		_mm_store_ps(M_l, M);
		M = _mm_and_ps(M, AbsMask128);
		__m128 Sum = _mm_hadd_ps(M, M);
		Sum = _mm_hadd_ps(Sum, Sum);

        for (int i = 0; i < 4; ++i)
    		assert((M_l[i] >= 0) || (M_l[i] < 0));
		return Sum.m128_f32[0];
	}
	else if (N == 2)
	{
		__m128i C = _mm_set1_epi32(*(reinterpret_cast<const int*>(C_l)));
		C = _mm_shuffle_epi8(C, UnpackMask);
		__m128 M = _mm_load_ps(M_l_1);
		M = _mm_xor_ps(_mm_castsi128_ps(C), M);
		M = _mm_add_ps(M, _mm_loadu_ps(M_l_2));
		_mm_storel_pd(reinterpret_cast<double*>(M_l), _mm_castps_pd(M));
		M = _mm_and_ps(M, AbsMask128);
		__m128 Sum = _mm_hadd_ps(M, M);
        for (int i = 0; i < 2; ++i)
    		assert((M_l[i] >= 0) || (M_l[i] < 0));
		return Sum.m128_f32[0];
	}
	else // N == 1
	{
		tBit u1 = C_l[0];
		if (u1) {
			M_l[0] = M_l_2[0] - M_l_1[0]; SUM_COUNT;
		}
		else {
			M_l[0] = M_l_2[0] + M_l_1[0]; SUM_COUNT;
		}
		assert((M_l[0] >= 0) || (M_l[0] < 0));
		return fabs(M_l[0]);
	}
}


/// Calculation of Q function:
/// Q(a, b) = sign(a)*sign(b)*min(a,b)
void SoftXOR(float *M_l, const float *M_l_1, const float *M_l_2, const uint32_t &N) {
	
	if (N >= 8) {
		CMP_BLOCK_COUNT(N);
#pragma unroll(4)
		for (unsigned beta = 0; beta < N; beta += 8)
		{
			// sign(a)*sign(b)*min(a,b) calculation
			// extract and remove signs
			__m256 M1 = _mm256_load_ps(M_l_1 + beta);
			__m256 M2 = _mm256_load_ps(M_l_2 + beta);
			__m256 M1_abs = _mm256_and_ps(M1, AbsMask);
			__m256 M2_abs = _mm256_and_ps(M2, AbsMask);
			// xor = mul
			__m256 signs = _mm256_xor_ps(M1, M2);
			signs = _mm256_andnot_ps(AbsMask, signs);
			// get minimum values
			__m256 M = _mm256_min_ps(M1_abs, M2_abs);
			// set result signs
			M = _mm256_or_ps(M, signs);
			_mm256_store_ps((M_l + beta), M);
            for (int i = 0; i < 8; ++i)
    			assert((M_l[beta + i] >= 0) || (M_l[beta + i] < 0));
			//assert(_CrtCheckMemory());
		}
	}
	else if (N == 4)
	{
		CMP_BLOCK_COUNT(N);
		// sign(a)*sign(b)*min(a,b) calculation
		// extract and remove signs
		__m128 M1 = _mm_load_ps(M_l_1);
		__m128 M2 = _mm_load_ps(M_l_2);
		__m128 M1_abs = _mm_and_ps(M1, AbsMask128);
		__m128 M2_abs = _mm_and_ps(M2, AbsMask128);
		// xor = mul
		__m128 signs = _mm_xor_ps(M1, M2);
		signs = _mm_andnot_ps(AbsMask128, signs);
		// get minimum values
		__m128 M = _mm_min_ps(M1_abs, M2_abs);
		// set result signs
		M = _mm_or_ps(M, signs);
		_mm_store_ps(M_l, M);
        for (int i = 0; i < 4; ++i)
    		assert((M_l[i] >= 0) || (M_l[i] < 0));
	}
	else if (N == 2)
	{
		CMP_BLOCK_COUNT(N);
		// sign(a)*sign(b)*min(a,b) calculation
		// extract and remove signs
		__m128 M1 = _mm_load_ps(M_l_1);
		__m128 M2 = _mm_loadu_ps(M_l_2);
		__m128 M1_abs = _mm_and_ps(M1, AbsMask128);
		__m128 M2_abs = _mm_and_ps(M2, AbsMask128);
		// xor = mul
		__m128 signs = _mm_xor_ps(M1, M2);
		signs = _mm_andnot_ps(AbsMask128, signs);
		// get minimum values
		__m128 M = _mm_min_ps(M1_abs, M2_abs);
		// set result signs
		M = _mm_or_ps(M, signs);
		_mm_storel_pd(reinterpret_cast<double*>(M_l), _mm_castps_pd(M));
        for (int i = 0; i < 2; ++i)
    		assert((M_l[i] >= 0) || (M_l[i] < 0));
	}
	else 
	{
		for(unsigned i=0;i<N;i++)
			SoftXOR2Elements(&M_l[i], M_l_1[i], M_l_2[i]);
	}
}

/// Calculation of Q function:
/// Q(a, b, c, d) = sign(a)*sign(b)*sign(c)*sign(d)*min(|a|,|b|,|c|,|d|), as well as sign(a)*sign(c)*min(|a|,|c|), sign(b)*sign(d)*min(|b|,|d|)
/// a,b,c,d are assumed to be located in the same array in consecutive blocks of size n
void SoftXOR(float *ABCD, float* AC, const float *A, const uint32_t &N) 
{
	if (N >= 8) {
		CMP_BLOCK_COUNT(3*N);
#pragma unroll(4)
		for (unsigned beta = 0; beta < N; beta += 8)
		{
			__m256 a = _mm256_load_ps(A + beta);
			__m256 b = _mm256_load_ps(A + N + beta);
			__m256 c = _mm256_load_ps(A + 2*N + beta);
			__m256 d = _mm256_load_ps(A + 3 * N + beta);

			__m256 a_abs = _mm256_and_ps(a, AbsMask);
			__m256 b_abs = _mm256_and_ps(b, AbsMask);
			__m256 c_abs = _mm256_and_ps(c, AbsMask);
			__m256 d_abs = _mm256_and_ps(d, AbsMask);

			// xor = mul
			__m256 ac_s = _mm256_xor_ps(a, c);
			__m256 bd_s = _mm256_xor_ps(b, d);
			ac_s = _mm256_andnot_ps(AbsMask, ac_s);
			bd_s = _mm256_andnot_ps(AbsMask, bd_s);
			__m256 abcd_s = _mm256_xor_ps(ac_s, bd_s);
			__m256 ac_m = _mm256_min_ps(a_abs, c_abs);
			__m256 bd_m = _mm256_min_ps(b_abs, d_abs);
			__m256 abcd_m = _mm256_min_ps(ac_m, bd_m);

			// set result signs
			abcd_m = _mm256_or_ps(abcd_m, abcd_s);
			ac_m = _mm256_or_ps(ac_m, ac_s);
			bd_m = _mm256_or_ps(bd_m, bd_s);
			_mm256_store_ps((ABCD + beta), abcd_m);
			_mm256_store_ps((AC + beta), ac_m);
			_mm256_store_ps((AC +N+ beta), bd_m);
			//assert(_CrtCheckMemory());
		}
	}
	else if (N == 4)
	{
		CMP_BLOCK_COUNT(3*N);
		// sign(a)*sign(b)*min(a,b) calculation
		__m128 a = _mm_load_ps(A );
		__m128 b = _mm_load_ps(A + N );
		__m128 c = _mm_load_ps(A + 2 * N );
		__m128 d = _mm_load_ps(A + 3 * N );

		__m128 a_abs = _mm_and_ps(a, AbsMask128);
		__m128 b_abs = _mm_and_ps(b, AbsMask128);
		__m128 c_abs = _mm_and_ps(c, AbsMask128);
		__m128 d_abs = _mm_and_ps(d, AbsMask128);

		// xor = mul
		__m128 ac_s = _mm_xor_ps(a, c);
		__m128 bd_s = _mm_xor_ps(b, d);
		ac_s = _mm_andnot_ps(AbsMask128, ac_s);
		bd_s = _mm_andnot_ps(AbsMask128, bd_s);
		__m128 abcd_s = _mm_xor_ps(ac_s, bd_s);
		__m128 ac_m = _mm_min_ps(a_abs, c_abs);
		__m128 bd_m = _mm_min_ps(b_abs, d_abs);
		__m128 abcd_m = _mm_min_ps(ac_m, bd_m);

		// set result signs
		abcd_m = _mm_or_ps(abcd_m, abcd_s);
		ac_m = _mm_or_ps(ac_m, ac_s);
		bd_m = _mm_or_ps(bd_m, bd_s);
		_mm_store_ps((ABCD ), abcd_m);
		_mm_store_ps((AC ), ac_m);
		_mm_store_ps((AC + N ), bd_m);
	}
	else if (N == 2)
	{
		CMP_BLOCK_COUNT(N);
		// sign(a)*sign(b)*min(a,b) calculation
		// extract and remove signs
		__m128 a = _mm_load_ps(A);
		__m128 b = _mm_loadu_ps(A + N);
		__m128 c = _mm_load_ps(A + 2 * N);
		__m128 d = _mm_loadu_ps(A + 3 * N);

		__m128 a_abs = _mm_and_ps(a, AbsMask128);
		__m128 b_abs = _mm_and_ps(b, AbsMask128);
		__m128 c_abs = _mm_and_ps(c, AbsMask128);
		__m128 d_abs = _mm_and_ps(d, AbsMask128);

		// xor = mul
		__m128 ac_s = _mm_xor_ps(a, c);
		__m128 bd_s = _mm_xor_ps(b, d);
		ac_s = _mm_andnot_ps(AbsMask128, ac_s);
		bd_s = _mm_andnot_ps(AbsMask128, bd_s);
		__m128 abcd_s = _mm_xor_ps(ac_s, bd_s);
		__m128 ac_m = _mm_min_ps(a_abs, c_abs);
		__m128 bd_m = _mm_min_ps(b_abs, d_abs);
		__m128 abcd_m = _mm_min_ps(ac_m, bd_m);

		// set result signs
		abcd_m = _mm_or_ps(abcd_m, abcd_s);
		ac_m = _mm_or_ps(ac_m, ac_s);
		bd_m = _mm_or_ps(bd_m, bd_s);
		_mm_storel_pd(reinterpret_cast<double*>(ABCD), _mm_castps_pd(abcd_m));
		_mm_storel_pd(reinterpret_cast<double*>(AC), _mm_castps_pd(ac_m));
		_mm_storel_pd(reinterpret_cast<double*>(AC+N), _mm_castps_pd(bd_m));
	}
	else // N == 1
	{
		ABCD;
		SoftXOR(AC, A,2);
		SoftXOR2Elements(ABCD,AC[0], AC[1]);
	};
}

/// Calculation of Q function:
/// Q(a, b) = sign(a)*sign(b)*min(a,b)
/// return \sum |Q(a_i,b_i)|
float SoftXORWithSum(float *M_l, const float *M_l_1, const float *M_l_2, const uint32_t &N) 
{
	

    if (N >= 8) 
    {
        CMP_BLOCK_COUNT(N);
        SUM_BLOCK_COUNT(N);
        __m256 Sum = _mm256_set1_ps(0);
#pragma unroll(4)
        for (unsigned beta = 0; beta < N; beta += 8)
        {
            // sign(a)*sign(b)*min(a,b) calculation
            // extract and remove signs
            __m256 M1 = _mm256_load_ps(M_l_1 + beta);
            __m256 M2 = _mm256_load_ps(M_l_2 + beta);
            __m256 M1_abs = _mm256_and_ps(M1, AbsMask);
            __m256 M2_abs = _mm256_and_ps(M2, AbsMask);
            // xor = mul
            __m256 signs = _mm256_xor_ps(M1, M2);
            signs = _mm256_andnot_ps(AbsMask, signs);
            // get minimum values
            __m256 M = _mm256_min_ps(M1_abs, M2_abs);
            Sum = _mm256_add_ps(Sum, M);
            // set result signs
            M = _mm256_or_ps(M, signs);
            _mm256_store_ps((M_l + beta), M);
            for (int i = 0; i < 8; ++i)
                assert((M_l[beta + i] >= 0) || (M_l[beta + i] < 0));
            //assert(_CrtCheckMemory());
        }
        Sum = _mm256_hadd_ps(Sum, Sum);
        Sum = _mm256_hadd_ps(Sum, Sum);
        return Sum.m256_f32[0] + Sum.m256_f32[4];

    }
    else if (N == 4)
    {
        CMP_BLOCK_COUNT(N);
        SUM_BLOCK_COUNT(N);
        // sign(a)*sign(b)*min(a,b) calculation
        // extract and remove signs
        __m128 M1 = _mm_load_ps(M_l_1);
        __m128 M2 = _mm_load_ps(M_l_2);
        __m128 M1_abs = _mm_and_ps(M1, AbsMask128);
        __m128 M2_abs = _mm_and_ps(M2, AbsMask128);
        // xor = mul
        __m128 signs = _mm_xor_ps(M1, M2);
        signs = _mm_andnot_ps(AbsMask128, signs);
        // get minimum values
        __m128 M = _mm_min_ps(M1_abs, M2_abs);
        //compute M[0]+M[1],m[2]+M[3]
        __m128 Sum = _mm_hadd_ps(M,M);
        Sum = _mm_hadd_ps(Sum,Sum);
        // set result signs
        M = _mm_or_ps(M, signs);
        _mm_store_ps(M_l, M);
        for (int i = 0; i < 4; ++i)
            assert((M_l[i] >= 0) || (M_l[i] < 0));
        return Sum.m128_f32[0] ;
    }
    else if (N == 2)
    {
        CMP_BLOCK_COUNT(N);
        SUM_BLOCK_COUNT(1);
        // sign(a)*sign(b)*min(a,b) calculation
        // extract and remove signs
        __m128 M1 = _mm_load_ps(M_l_1);
        __m128 M2 = _mm_loadu_ps(M_l_2);
        __m128 M1_abs = _mm_and_ps(M1, AbsMask128);
        __m128 M2_abs = _mm_and_ps(M2, AbsMask128);
        // xor = mul
        __m128 signs = _mm_xor_ps(M1, M2);
        signs = _mm_andnot_ps(AbsMask128, signs);
        // get minimum values
        __m128 M = _mm_min_ps(M1_abs, M2_abs);
        //compute M[0]+M[1]
        __m128 Sum = _mm_hadd_ps(M, M);
        // set result signs
        M = _mm_or_ps(M, signs);
        _mm_storel_pd(reinterpret_cast<double*>(M_l), _mm_castps_pd(M));
        for (int i = 0; i < 2; ++i)
            assert((M_l[i] >= 0) || (M_l[i] < 0));
        return Sum.m128_f32[0];
    }
    else // N == 1
    {
        CMP_COUNT;
        if (M_l_1[0] < M_l_2[0])
        {
            CMP_COUNT;
            if (M_l_1[0] > -M_l_2[0])
            {
                M_l[0] = M_l_1[0];
            }
            else
            {
                M_l[0] = -M_l_2[0];
            };
        }
        else
        {
            CMP_COUNT;
            if (M_l_1[0] > -M_l_2[0])
            {
                M_l[0] = M_l_2[0];
            }
            else
            {
                M_l[0] = -M_l_1[0];
            }
        }
        assert((M_l[0] >= 0) || (M_l[0] < 0));
        return fabs(M_l[0]);
    }

};

///compute 'sum_i |L_i|
float  SumAbs(const float *pLLR, unsigned N)
{
	__m256 Sum = _mm256_setzero_ps();
	unsigned i = 0;
	for (i; i < (N&~7); i += 8)
	{
		__m256 llr = _mm256_load_ps(pLLR+i);
		llr = _mm256_and_ps(llr, AbsMask);
		Sum = _mm256_add_ps(Sum, llr);
	};
	Sum = _mm256_hadd_ps(Sum, Sum);
	Sum = _mm256_hadd_ps(Sum, Sum);
	float sum = Sum.m256_f32[0] + Sum.m256_f32[4];
	for (i; i < N; i++)
	{
		sum += fabs(pLLR[i]);
	};
	return sum;
}

//change LLR signs using a given mask
//if C[i]=1, L1[i]=-L[i], else L1[i]=L[i]
void ChangeLLRSigns(float *pDestLLRs, const float *pSrcLLRs, const tBit *pMask, unsigned N)
{
	if (N>= 8)
	{
		for (unsigned beta = 0; beta < N; beta += 8)
		{
			__m128i C_1 = _mm_castpd_si128(_mm_loadu_pd(reinterpret_cast<const double*>(pMask + beta)));
			__m128i C_2 = _mm_castpd_si128(_mm_loadu_pd(reinterpret_cast<const double*>(pMask + beta + 4)));
			C_1 = _mm_shuffle_epi8(C_1, UnpackMask);
			C_2 = _mm_shuffle_epi8(C_2, UnpackMask);
			__m256 C = _mm256_set_m128(_mm_castsi128_ps(C_2), _mm_castsi128_ps(C_1));
			__m256 M = _mm256_load_ps(pSrcLLRs + beta);
			M = _mm256_xor_ps(C, M);
			_mm256_store_ps(pDestLLRs + beta, M);
		}
	}else
	if (N==4)
	{
		__m128i C = _mm_set1_epi32(*(reinterpret_cast<const uint32_t*>(pMask)));
		C = _mm_shuffle_epi8(C, UnpackMask);
		__m128 M = _mm_load_ps(pSrcLLRs);
		M = _mm_xor_ps(_mm_castsi128_ps(C), M);
		_mm_store_ps(pDestLLRs , M);
	}
	else if (N == 2)
	{
		__m128i C = _mm_set1_epi32(*(reinterpret_cast<const int32_t*>(pMask)));
		C = _mm_shuffle_epi8(C, UnpackMask);
		__m128 M = _mm_load_ps(pSrcLLRs);
		M = _mm_xor_ps(_mm_castsi128_ps(C), M);
		_mm_storel_pd(reinterpret_cast<double*>(pDestLLRs), _mm_castps_pd(M));
	}
	else 
	{
		for (unsigned i = 0; i < N; i++)
		{
			if (pMask[i])
				pDestLLRs[i] = -pSrcLLRs[i];
			else
				pDestLLRs[i] = pSrcLLRs[i];
		};
	}

};

#elif defined(USE_INT32)
void SoftCombine(MType *M_l, const MType *M_l_1, const MType *M_l_2,const tBit *C_l, const uint32_t &N) {
	for (uint32_t beta = 0; beta < N; ++beta)
	{
        SoftCombine2Elements(&M_l[beta], M_l_1[beta], M_l_2[beta], C_l[beta]);
	}
	

}

void SoftXOR(MType *M_l, const MType *M_l_1, const MType *M_l_2, const uint32_t &N) {
	for (uint32_t beta = 0; beta < N; beta++)
	{
        SoftXOR2Elements(&M_l[beta], M_l_1[beta], M_l_2[beta]);
	};
	
}

int SumAbs(const int32_t *pLLR, unsigned N)
{
	__m128i Sum = _mm_setzero_si128();
	unsigned i = 0;
	for (i; i < N&~7; i += 8)
	{
		__m128i llr0 = _mm_load_si128(reinterpret_cast<const __m128i*>(pLLR + i));     // load
		__m128i llr1 = _mm_load_si128(reinterpret_cast<const __m128i*>(pLLR + i + 4));
		llr0 = _mm_abs_epi32(llr0);  // get abs values
		llr1 = _mm_abs_epi32(llr1);
		llr0 = _mm_hadd_epi32(llr0, llr1);
		Sum = _mm_add_epi32(Sum, llr0);
	};
	Sum = _mm_hadd_epi32(Sum, Sum);
	int sum = Sum.m128i_i32[0] + Sum.m128i_i32[1];
	for (i; i < N; i++)
	{
		sum += abs(pLLR[i]);
	};
	return sum;
}

void ChangeLLRSigns(MType* pDestLLRs,const MType* pSrcLLRs, const tBit* pMask, unsigned  N)
{
	for (unsigned i = 0; i < N; i++)
	{
		if (pMask[i])
		{
			pDestLLRs[i] = umin(pSrcLLRs[i]);
		}
		else
		{
			pDestLLRs[i] = pSrcLLRs[i];
		}
	}
};
#elif defined(USE_INT8)
#ifdef AVX2
const __m256i Mask = _mm256_set1_epi8(1);
const __m256i Mask1 = _mm256_set1_epi8('\x80');
const __m256i Mask2 = _mm256_set1_epi8('\xFE');
const __m256i minValue = _mm256_set1_epi8(-127);

/// Calculation of G function: 
/// g(lambda, beta) = (-1)^C_l[beta][0]*M_{l-1}[2*beta] + M_{l-1}[2*beta+1]
inline void SoftCombine_m256(const __m256i & M1, const __m256i & M2, const __m256i & C, __m256i * out) {
	__m256i sum1 = _mm256_adds_epi8(M2, M1);
	__m256i sum2 = _mm256_subs_epi8(M2, M1);
	__m256i M = _mm256_blendv_epi8(sum1, sum2, C);
	*out = _mm256_max_epi8(M, minValue);
}

void SoftCombine(int8_t *M_l, const int8_t *M_l_1, const int8_t *M_l_2, const tBit *C_l, const uint32_t &N) {
	SUM_BLOCK_COUNT(N);
	

	__m256i C = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(C_l));
	__m256i M1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_1));
	__m256i M2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_2));
	__m256i M; SoftCombine_m256(M1, M2, C, &M);
	switch (N)
	{
	case 1:
		*M_l = M.m256i_i8[0];
		break;
	case 2:
		*reinterpret_cast<int16_t *>(M_l) = M.m256i_i16[0];
		break;
	case 4:
		*reinterpret_cast<int32_t *>(M_l) = M.m256i_i32[0];
		break;
	case 8:
		*reinterpret_cast<int64_t *>(M_l) = M.m256i_i64[0];
		break;
	case 16:
		*reinterpret_cast<__m128i *>(M_l) = _mm256_extracti128_si256(M, 0x0);
		break;
	default:
		_mm256_store_si256(reinterpret_cast<__m256i *>(M_l), M);
		for (uint32_t beta = 32; beta < N; beta += 32)
		{
			__m256i C = _mm256_load_si256(reinterpret_cast<const __m256i *>(C_l + beta));
			__m256i M1 = _mm256_load_si256(reinterpret_cast<const __m256i *>(M_l_1 + beta));
			__m256i M2 = _mm256_load_si256(reinterpret_cast<const __m256i *>(M_l_2 + beta));
			__m256i M;
			SoftCombine_m256(M1, M2, C, &M);
			_mm256_store_si256(reinterpret_cast<__m256i *>(M_l + beta), M);
		}
	}

}

/// Calculation of Q function:
/// Q(a, b) = sign(a)*sign(b)*min(a,b)
inline void SoftXOR_m256(const __m256i & M1, const __m256i & M2, __m256i *absM, __m256i * out) {
	// sign(a)*sign(b)*min(a,b) calculation
	// extract and remove signs
	__m256i M1_abs = _mm256_abs_epi8(M1);
	__m256i M2_abs = _mm256_abs_epi8(M2);
	// xor sign = mul sign
	__m256i signs = _mm256_xor_si256(M1, M2);
	signs = _mm256_or_si256(signs, Mask);
	// get minimum values
	*absM = _mm256_min_epi8(M1_abs, M2_abs);
	// set result signs
	*out = _mm256_sign_epi8(*absM, signs);
}

void SoftXOR(int8_t *M_l, const int8_t *M_l_1, const int8_t *M_l_2, const uint32_t &N) {
	CMP_BLOCK_COUNT(N);
	__m256i M1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_1));
	__m256i M2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_2));
	__m256i M, absM; SoftXOR_m256(M1, M2, &absM, &M);
	switch (N)
	{
	case 1:
		*M_l = M.m256i_i8[0];
		break;
	case 2:
		*reinterpret_cast<int16_t *>(M_l) = M.m256i_i16[0];
		break;
	case 4:
		*reinterpret_cast<int32_t *>(M_l) = M.m256i_i32[0];
		break;
	case 8:
		*reinterpret_cast<int64_t *>(M_l) = M.m256i_i64[0];
		break;
	case 16:
		*reinterpret_cast<__m128i *>(M_l) = _mm256_extracti128_si256(M, 0x0);
		break;
	default:
		_mm256_store_si256(reinterpret_cast<__m256i *>(M_l), M);
		for (uint32_t beta = 32; beta < N; beta += 32)
		{
			__m256i M1 = _mm256_load_si256(reinterpret_cast<const __m256i *>(M_l_1 + beta));
			__m256i M2 = _mm256_load_si256(reinterpret_cast<const __m256i *>(M_l_2 + beta));
			__m256i M, absM;
			SoftXOR_m256(M1, M2, &absM, &M);
			_mm256_store_si256(reinterpret_cast<__m256i *>(M_l + beta), M);
		}
	}
	

}

int32_t SumAbs(const int8_t *pLLR, unsigned N)
{
	__m256i sum_256 = _mm256_setzero_si256();
	int32_t index = 0;
	__m256i llr;
	while (index < (N & (~0x1F))) {
		llr = _mm256_load_si256(reinterpret_cast<const __m256i*>(pLLR + index));
		llr = _mm256_abs_epi8(llr);
		sum_256 = _mm256_add_epi8(sum_256, llr);
		index += 32;
	}

	int sum = 0;
	for (int32_t i = 0; i < 32; ++i) {
		sum += sum_256.m256i_i8[i];
	}

	const int8_t *it = pLLR + index;
	for (int32_t i = index; i < N; ++i) {
		sum += abs(*it);
		++it;
	}
	return sum;
}

//change LLR signs using a given mask
//if C[i]=1, L1[i]=-L[i], else L1[i]=L[i]
void ChangeLLRSigns(int8_t* pDestLLRs, const int8_t* pSrcLLRs, const tBit* pMask, unsigned  N)
{
	unsigned index = 0;
	while (index < (N & (~0x1F)))
	{
		__m256i srcLLRs = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(pSrcLLRs + index));
		__m256i dstLLRs = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(pDestLLRs + index));
		__m256i mask = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(pMask + index));
		mask = _mm256_or_si256(mask, Mask);
		dstLLRs = _mm256_sign_epi8(srcLLRs, mask);
		_mm256_storeu_si256(reinterpret_cast<__m256i *>(pDestLLRs + index), dstLLRs);
		index += 32;
	}

	for (unsigned i = index; i < N; i++)
	{
		if ((bool)pMask[i])
		{
			pDestLLRs[i] = umin(pSrcLLRs[i]);
		}
		else
		{
			pDestLLRs[i] = pSrcLLRs[i];
		}
	}
};

inline void SoftXORWithSum2_m256(const __m256i & M1, const __m256i & M2, __m256i * out, int32_t * sum) {
	__m256i absM;
	SoftXOR_m256(M1, M2, &absM, out);
	absM = _mm256_cvtepi8_epi16(_mm256_extracti128_si256(absM, 0x0));
	absM = _mm256_hadd_epi16(absM, absM); // 2 -> 1
	*sum += absM.m256i_u16[0];
}

inline void SoftXORWithSum4_m256(const __m256i & M1, const __m256i & M2, __m256i * out, int32_t * sum) {
	__m256i absM;
	SoftXOR_m256(M1, M2, &absM, out);
	absM = _mm256_cvtepi8_epi16(_mm256_extracti128_si256(absM, 0x0));
	absM = _mm256_hadd_epi16(absM, absM); // 4 -> 2
	absM = _mm256_hadd_epi16(absM, absM); // 2 -> 1
	*sum += absM.m256i_u16[0];
}

inline void SoftXORWithSum8_m256(const __m256i & M1, const __m256i & M2, __m256i * out, int32_t * sum) {
	__m256i absM;
	SoftXOR_m256(M1, M2, &absM, out);
	absM = _mm256_cvtepi8_epi16(_mm256_extracti128_si256(absM, 0x0));
	absM = _mm256_hadd_epi16(absM, absM); // 8 -> 4
	absM = _mm256_hadd_epi16(absM, absM); // 4 -> 2
	absM = _mm256_hadd_epi16(absM, absM); // 2 -> 1
	*sum += absM.m256i_u16[0];
}

inline void SoftXORWithSum16_m256(const __m256i & M1, const __m256i & M2, __m256i * out, int32_t * sum) {
	__m256i absM;
	SoftXOR_m256(M1, M2, &absM, out);
	absM = _mm256_cvtepi8_epi16(_mm256_extracti128_si256(absM, 0x0));
	absM = _mm256_hadd_epi16(absM, absM); // 16 -> 8
	absM = _mm256_hadd_epi16(absM, absM); // 8  -> 4
	absM = _mm256_hadd_epi16(absM, absM); // 4  -> 2
	absM = _mm256_hadd_epi16(absM, absM); // 2  -> 1
	*sum += absM.m256i_u16[0];
}

inline void SoftXORWithSum32_m256(const __m256i & M1, const __m256i & M2, __m256i * out, int32_t * sum) {
	__m256i absM;
	SoftXOR_m256(M1, M2, &absM, out);
	absM = _mm256_hadd_epi16(absM, absM); // 32 -> 16
	absM = _mm256_permute4x64_epi64(absM, 0x08);
	absM = _mm256_cvtepi8_epi16(_mm256_extracti128_si256(absM, 0x0));
	absM = _mm256_hadd_epi16(absM, absM); // 16 -> 8
	absM = _mm256_hadd_epi16(absM, absM); // 8  -> 4
	absM = _mm256_hadd_epi16(absM, absM); // 4  -> 2
	absM = _mm256_hadd_epi16(absM, absM); // 2  -> 1
	*sum += absM.m256i_u16[0];
}

int32_t SoftXORWithSum(int8_t *M_l, const int8_t *M_l_1, const int8_t *M_l_2, const uint32_t &N)
{
	CMP_BLOCK_COUNT(N);
	__m256i M1, M2, M, absM;
	int32_t Sum = 0;
	switch (N)
	{
	case 1:
		M1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_1));
		M2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_2));
		SoftXOR_m256(M1, M2, &absM, &M);
		*M_l = M.m256i_i8[0];
		return abs_value(M_l[0]);

	case 2:
		M1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_1));
		M2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_2));
		SoftXORWithSum2_m256(M1, M2, &M, &Sum);
		*reinterpret_cast<int16_t *>(M_l) = M.m256i_i16[0];
		return Sum;

	case 4:
		M1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_1));
		M2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_2));
		SoftXORWithSum4_m256(M1, M2, &M, &Sum);
		*reinterpret_cast<int32_t *>(M_l) = M.m256i_i32[0];
		return Sum;

	case 8:
		M1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_1));
		M2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_2));
		SoftXORWithSum8_m256(M1, M2, &M, &Sum);
		*reinterpret_cast<int64_t *>(M_l) = M.m256i_i64[0];
		return Sum;
	case 16:
		M1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_1));
		M2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_2));
		SoftXORWithSum8_m256(M1, M2, &M, &Sum);
		*reinterpret_cast<__m128i *>(M_l) = _mm256_extracti128_si256(M, 0x0);
		return Sum;

	default:
		for (uint32_t beta = 0; beta < N; beta += 32)
		{
			M1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_1 + beta));
			M2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(M_l_2 + beta));
			SoftXORWithSum32_m256(M1, M2, &M, &Sum);
			_mm256_store_si256(reinterpret_cast<__m256i *>(M_l + beta), M);
		}
		return Sum;
	}
}
#else
const __m128i Mask = _mm_set1_epi8(1);
const __m128i minValue = _mm_set1_epi8(-127);

/// Calculation of G function: 
/// g(lambda, beta) = (-1)^C_l[beta][0]*M_{l-1}[2*beta] + M_{l-1}[2*beta+1]
inline void SoftCombine_m128(const __m128i & M1, const __m128i & M2, const __m128i & C, __m128i * out) {
    //__m128i sign_mask = _mm_or_si128(C, Mask);       // 0 -> 1
    //__m128i M = _mm_sign_epi8(M1, sign_mask);        // negate elements
    //M = _mm_adds_epi8(M, M2);                        // summation
    __m128i sum1 = _mm_adds_epi8(M2, M1);
    __m128i sum2 = _mm_subs_epi8(M2, M1);
    __m128i M = _mm_blendv_epi8(sum1, sum2, C);
    *out = _mm_max_epi8(M, minValue);
    //__m128i correction = _mm_cmpeq_epi8(M, Mask1);
    //correction = _mm_and_si128(correction, Mask2);
    //*out = _mm_sub_epi8(M, correction);              // -128 -> -127
}

void SoftCombine(int8_t *M_l, const int8_t *M_l_1, const int8_t *M_l_2, const tBit *C_l, const uint32_t &N) {
    SUM_BLOCK_COUNT(N);
    __m128i C = _mm_loadu_si128(reinterpret_cast<const __m128i *>(C_l));
    __m128i M1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_1));
    __m128i M2 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_2));
    __m128i M; SoftCombine_m128(M1, M2, C, &M);
    switch (N)
    {
    case 1:
        *M_l = M.m128i_i8[0];
        break;
    case 2:
        *reinterpret_cast<int16_t *>(M_l) = M.m128i_i16[0];
        break;
    case 4:
        *reinterpret_cast<int32_t *>(M_l) = M.m128i_i32[0];
        break;
    case 8:
        _mm_storel_epi64(reinterpret_cast<__m128i *>(M_l), M);
        break;
    default:
        _mm_store_si128(reinterpret_cast<__m128i *>(M_l), M);
        for (uint32_t beta = 16; beta < N; beta += 16)
        {
            C = _mm_load_si128(reinterpret_cast<const __m128i *>(C_l + beta));
            M1 = _mm_load_si128(reinterpret_cast<const __m128i *>(M_l_1 + beta));
            M2 = _mm_load_si128(reinterpret_cast<const __m128i *>(M_l_2 + beta));
            SoftCombine_m128(M1, M2, C, &M);
            _mm_store_si128(reinterpret_cast<__m128i *>(M_l + beta), M);
        }
    }
}

/// Calculation of Q function:
/// Q(a, b) = sign(a)*sign(b)*min(a,b)
inline void SoftXOR_m128(const __m128i & M1, const __m128i & M2, __m128i * absM, __m128i * out) {
    // sign(a)*sign(b)*min(a,b) calculation
    // extract and remove signs
    __m128i M1_abs = _mm_abs_epi8(M1);
    __m128i M2_abs = _mm_abs_epi8(M2);
    // xor sign = mul sign
    __m128i signs = _mm_xor_si128(M1, M2);
    signs = _mm_or_si128(signs, Mask);
    // get minimum values
    *absM = _mm_min_epi8(M1_abs, M2_abs);
    // set result signs
    *out = _mm_sign_epi8(*absM, signs);
}

void SoftXOR(int8_t *M_l, const int8_t *M_l_1, const int8_t *M_l_2, const uint32_t &N) {
    CMP_BLOCK_COUNT(N);
    __m128i M1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_1));
    __m128i M2 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_2));
    __m128i M, absM; SoftXOR_m128(M1, M2, &absM, &M);
    switch (N)
    {
    case 1:
        *M_l = M.m128i_i8[0];
        break;
    case 2:
        *reinterpret_cast<int16_t *>(M_l) = M.m128i_i16[0];
        break;
    case 4:
        *reinterpret_cast<int32_t *>(M_l) = M.m128i_i32[0];
        break;
    case 8:
        _mm_storel_epi64(reinterpret_cast<__m128i *>(M_l), M);
        break;
    default:
        _mm_store_si128(reinterpret_cast<__m128i *>(M_l), M);
        for (uint32_t beta = 16; beta < N; beta += 16)
        {
            M1 = _mm_load_si128(reinterpret_cast<const __m128i *>(M_l_1 + beta));
            M2 = _mm_load_si128(reinterpret_cast<const __m128i *>(M_l_2 + beta));
            SoftXOR_m128(M1, M2, &absM, &M);
            _mm_store_si128(reinterpret_cast<__m128i *>(M_l + beta), M);
        }
    }
}

int32_t SumAbs(const int8_t *pLLR, unsigned N)
{
    __m128i sum_128 = _mm_setzero_si128();
    int32_t index = 0;
    __m128i llr;
    while (index < (N & (~0xF))) {
        llr = _mm_load_si128(reinterpret_cast<const __m128i*>(pLLR + index));
        llr = _mm_abs_epi8(llr);
        sum_128 = _mm_add_epi8(sum_128, llr);
        index += 16;
    }

    int sum = 0;
    for (int32_t i = 0; i < 16; ++i) {
        sum += sum_128.m128i_i8[i];
    }

    const int8_t *it = pLLR + index;
    for (int32_t i = index; i < N; ++i) {
        sum += abs(*it);
        ++it;
    }
    return sum;
}
//change LLR signs using a given mask
//if C[i]=1, L1[i]=-L[i], else L1[i]=L[i]
void ChangeLLRSigns(int8_t* pDestLLRs, const int8_t* pSrcLLRs, const tBit* pMask, unsigned  N)
{
	unsigned index = 0;
	while (index < (N & (~0xF)))
	{
		__m128i srcLLRs = _mm_loadu_si128(reinterpret_cast<const __m128i *>(pSrcLLRs + index));
		__m128i dstLLRs = _mm_loadu_si128(reinterpret_cast<const __m128i *>(pDestLLRs + index));
		__m128i mask = _mm_loadu_si128(reinterpret_cast<const __m128i *>(pMask + index));
		//mask = _mm_slli_epi64(mask, 7);
		mask = _mm_or_si128(mask, Mask);
		dstLLRs = _mm_sign_epi8(srcLLRs, mask);
		_mm_storeu_si128(reinterpret_cast<__m128i *>(pDestLLRs + index), dstLLRs);
		index += 16;
	}

	for (unsigned i = index; i < N; i++)
	{
		if (pMask[i])
		{
			pDestLLRs[i] = umin(pSrcLLRs[i]);
		}
		else
		{
			pDestLLRs[i] = pSrcLLRs[i];
		}
	}
};
inline void SoftXORWithSum2_m128(const __m128i & M1, const __m128i & M2, __m128i * out, int32_t * sum) {
	__m128i absM;
	SoftXOR_m128(M1, M2, &absM, out);
	absM = _mm_cvtepi8_epi16(absM);
	absM = _mm_hadd_epi16(absM, absM); // 2 -> 1
	*sum += absM.m128i_u16[0];
}

inline void SoftXORWithSum4_m128(const __m128i & M1, const __m128i & M2, __m128i * out, int32_t * sum) {
	__m128i absM;
	SoftXOR_m128(M1, M2, &absM, out);
	absM = _mm_cvtepi8_epi16(absM);
	absM = _mm_hadd_epi16(absM, absM); // 4 -> 2
	absM = _mm_hadd_epi16(absM, absM); // 2 -> 1
	*sum += absM.m128i_u16[0];
}

inline void SoftXORWithSum8_m128(const __m128i & M1, const __m128i & M2, __m128i * out, int32_t * sum) {
	__m128i absM;
	SoftXOR_m128(M1, M2, &absM, out);
	absM = _mm_cvtepi8_epi16(absM);
	absM = _mm_hadd_epi16(absM, absM); // 8 -> 4
	absM = _mm_hadd_epi16(absM, absM); // 4 -> 2
	absM = _mm_hadd_epi16(absM, absM); // 2 -> 1
	*sum += absM.m128i_u16[0];
}

inline void SoftXORWithSum16_m128(const __m128i & M1, const __m128i & M2, __m128i * out, int32_t * sum) {
	__m128i M;
	SoftXOR_m128(M1, M2, &M, out);
	// Warning! overflow.
	M = _mm_hadd_epi16(M, M); // 16 -> 8
	M = _mm_cvtepi8_epi16(M);
	M = _mm_hadd_epi16(M, M); // 8  -> 4
	M = _mm_hadd_epi16(M, M); // 4  -> 2
	M = _mm_hadd_epi16(M, M); // 2  -> 1
	*sum += M.m128i_u16[0];
}

int32_t SoftXORWithSum(int8_t *M_l, const int8_t *M_l_1, const int8_t *M_l_2, const uint32_t &N) 
{
	CMP_BLOCK_COUNT(N);
	__m128i M1, M2, M, absM;
	int32_t Sum = 0;
	switch (N)
	{
	case 1:
		M1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_1));
		M2 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_2));
		SoftXOR_m128(M1, M2, &absM, &M);
		*M_l = M.m128i_i8[0];
		return abs_value(M_l[0]);

	case 2:
		M1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_1));
		M2 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_2));
		SoftXORWithSum2_m128(M1, M2, &M, &Sum);
		*reinterpret_cast<int16_t *>(M_l) = M.m128i_i16[0];
		return Sum;

	case 4:
		M1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_1));
		M2 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_2));
		SoftXORWithSum4_m128(M1, M2, &M, &Sum);
		*reinterpret_cast<int32_t *>(M_l) = M.m128i_i32[0];
		return Sum;

	case 8:
		M1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_1));
		M2 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(M_l_2));
		SoftXORWithSum8_m128(M1, M2, &M, &Sum);
		_mm_storel_epi64(reinterpret_cast<__m128i *>(M_l), M);
		return Sum;

	default:
		for (uint32_t beta = 0; beta < N; beta += 16)
		{
			M1 = _mm_load_si128(reinterpret_cast<const __m128i *>(M_l_1 + beta));
			M2 = _mm_load_si128(reinterpret_cast<const __m128i *>(M_l_2 + beta));
			SoftXORWithSum16_m128(M1, M2, &M, &Sum);
			_mm_store_si128(reinterpret_cast<__m128i *>(M_l + beta), M);
		}
		return Sum;
	}
}
#endif
#else
#error not implemented
#endif
#else
void SoftCombine(MType *M_l, const MType *M_l_1, const MType *M_l_2,const tBit *C_l, const uint32_t &N) {
    for (uint32_t beta = 0; beta < N; ++beta)
    {
        SoftCombine2Elements(&M_l[beta], M_l_1[beta], M_l_2[beta], C_l[beta]);
    }
}

void SoftXOR(MType *M_l, const MType *M_l_1, const MType *M_l_2, const uint32_t &N) {
    for (uint32_t beta = 0; beta < N; beta++)
    {
        SoftXOR2Elements(&M_l[beta], M_l_1[beta], M_l_2[beta]);
    };
}

KType SoftXORWithSum(MType *M_l, const MType *M_l_1, const MType *M_l_2, const uint32_t &N) 
{
    KType Sum = 0;
    for (uint32_t beta = 0; beta < N; beta++)
    {
        SoftXOR2Elements(&M_l[beta], M_l_1[beta], M_l_2[beta]);
        Sum += abs_value(M_l[beta]);
    };
    return Sum;
}

//change LLR signs using a given mask
//if C[i]=1, L1[i]=-L[i], else L1[i]=L[i]
void ChangeLLRSigns(MType* pDestLLRs, const MType* pSrcLLRs, const tBit* pMask, unsigned  N)
{
    for (unsigned i = 0; i < N; i++)
    {
        pDestLLRs[i] = (pMask[i]) ? -pSrcLLRs[i] : pSrcLLRs[i];
    };
};

#if defined(USE_FLOAT)
float SumAbs(const float* L, unsigned N)
{
    float S = fabs(L[0]);
    SUM_BLOCK_COUNT(N - 1);
    for (unsigned i = 1; i < N; i++)
    {
        S += fabs(L[i]);
    };
    return S;
}
#elif defined(USE_INT32) || defined(USE_INT8)
int32_t SumAbs(const int32_t* L, unsigned N)
{
    int32_t S = abs(L[0]);;
    SUM_BLOCK_COUNT(N - 1);
    for (unsigned i = 1; i < N; i++)
    {
        S += abs(L[i]);
    };
    return S;
}

int32_t SumAbs(const int8_t* L, unsigned N)
{
    int32_t S = abs(L[0]);;
    SUM_BLOCK_COUNT(N - 1);
    for (unsigned i = 1; i < N; i++)
    {
        S += abs(L[i]);
    };
    return S;
}
#else
#error not implemented
#endif
#endif


void SoftXOREvenOdd(MType *M_l, const MType *M_l_1, const uint32_t &N) {
	const MType * M_l_2 = M_l_1 + N;
	for (unsigned i = 0; i < N; ++i) {
		SoftXOR2Elements(&M_l[i], M_l_1[2 * i], M_l_1[2 * i + 1]);
	}
}

//permuted softXor, induced by different derivative order
void SoftXOR_0312(MType *M_l, const MType *M_l_1, const uint32_t &N) {
	for (unsigned i = 0; i < N / 2; ++i) {
		SoftXOR2Elements(&M_l[2 * i], M_l_1[4 * i], M_l_1[4 * i + 3]);
		SoftXOR2Elements(&M_l[2 * i + 1], M_l_1[4 * i + 1], M_l_1[4 * i + 2]);
	}
}

//permuted softCombine, induced by different derivative order
void SoftCombineEvenOdd(MType *M_l, const MType *M_l_1, const tBit *C_l, const uint32_t &N) {
	const MType * M_l_2 = M_l_1 + N;
	for (unsigned i = 0; i < N; ++i) {
		SoftCombine2Elements(&M_l[i], M_l_1[2 * i], M_l_1[2 * i + 1], C_l[i]);
	}
}
