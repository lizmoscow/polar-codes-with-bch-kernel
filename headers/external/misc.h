#ifndef MISC_H
#define MISC_H
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <exception>
#include <immintrin.h>
#include <nmmintrin.h>
#include <csetjmp>
#include <algorithm>
#include <functional>

#include <complex>
#include "SeqConfig.h"



#ifdef DLL_EXPORTS
#pragma message("dll exports are on")
#define DLL_API __declspec(dllexport)
#elif defined(DLL_IMPORTS)
#define DLL_API __declspec(dllimport)
#else
#define DLL_API
#endif 

//maximal finite field extension
#define MAXFIELDEXTENSION 16
//one symbol of a codeword
typedef unsigned __int64 Word;
#define ONE  0x01ui64
#define ZERO 0x00ui64

typedef unsigned GFValue;
typedef double Float;
typedef std::complex<double> Complex;
typedef Complex tChannelSymbol;

struct tPairOfProbabilities {
	double P[2];
};
struct tPairOfBool
{
	tBit C[2];
};

typedef std::pair<KType, unsigned> tPairOfMTypeUnsigned;
typedef std::pair<double, unsigned> tPairOfDoubleUnsigned;



extern std::jmp_buf JUMP_BUFFER;

#define ERROR_MESSAGE_BUFFER_SIZE 1024
extern char ErrorMessageBuffer[ERROR_MESSAGE_BUFFER_SIZE];

class Exception :public std::exception
{
public:
	Exception(const char* Format, ...)
	{
		va_list args;
		va_start(args, Format);
		vsprintf_s(ErrorMessageBuffer, ERROR_MESSAGE_BUFFER_SIZE-1, Format, args);
		va_end(args);
		*(std::exception*)this=std::exception(ErrorMessageBuffer);
	};

	Exception() {
	}
};

//operation counters
extern unsigned long long SumCount;
extern unsigned long long CmpCount;
extern unsigned long long MulCount;
extern unsigned long long ExpCount;
extern unsigned long long CpyCount;
//extern unsigned long long LatencyCount;

#ifdef OPERATION_COUNTING
#define SUM_COUNT SumCount++
#define MUL_COUNT MulCount++
#define EXP_COUNT ExpCount++
#define CMP_COUNT CmpCount++
#define CPY_COUNT CpyCount++
//#define LATENCY_COUNT LatencyCount++
#define SUM_BLOCK_COUNT(N) SumCount+=N
#define MUL_BLOCK_COUNT(N) MulCount+=N
#define EXP_BLOCK_COUNT(N) ExpCount+=N
#define CMP_BLOCK_COUNT(N) CmpCount+=N
#define CPY_BLOCK_COUNT(N) CpyCount+=N
//#define LATENCY_BLOCK_COUNT(N) LatencyCount+=N
#else

#define SUM_COUNT 
#define CMP_COUNT 
#define EXP_COUNT
#define MUL_COUNT
#define CPY_COUNT
#define SUM_BLOCK_COUNT(N)
#define CMP_BLOCK_COUNT(N)
#define EXP_BLOCK_COUNT(N)
#define MUL_BLOCK_COUNT(N)
#define CPY_BLOCK_COUNT(N)

#endif

enum Counters {
	CNT_MEMORY_OUTER,
	CNT_MEMORY_CS,
	CNT_MEMORY_ERROR,
	CNT_UNKNOWN,
	CNT_END
};

extern Counters CurrentCounter;
extern unsigned long long pCounters[CNT_END];

#ifdef MEMORY_COUNTING
#define CNT_START(Counter) CurrentCounter = Counter
#define CNT_STOP() CurrentCounter = CNT_UNKNOWN
#define CNT_INFO(Counter) pCounters[Counter]
#define CNT_RESET(Counter) pCounters[Counter] = 0
#define CNT(Counter, N)  pCounters[Counter] += N
#else
#define CNT_START(Counter)
#define CNT_STOP()
#define CNT_INFO(Counter) 0
#define CNT_RESET(Counter) 
#define CNT(Counter, N)
#endif




///Gallager's function log tanh(x/2)
template<class F> F GallagerFunction(const F& x)
{
		if (x < 12)
		{
			if (x < 1E-5)
				return log(2.0) - log(x);
			else
				return -log(tanh(x / 2));
		}
		else
			return 2 / (exp(x) - 1);//use an approximate expression for very large argument values
}

#ifdef OPERATION_COUNTING
class CPairComparator :public std::greater<tPairOfMTypeUnsigned>
{

public:
	bool operator()(const tPairOfMTypeUnsigned& _Left, const tPairOfMTypeUnsigned& _Right) const
	{
		CMP_COUNT;
		return std::greater<tPairOfMTypeUnsigned>::operator()(_Left, _Right);
	};

};
#else 
typedef  std::greater<tPairOfMTypeUnsigned> CPairComparator;
#endif

/*
double pow(double x, unsigned y)
{
	return pow(x, (int)y);
};*/



//CRC-16
unsigned CRC16(const tBit* pSrc, unsigned K);

///compute binomial coefficient
unsigned long long Binomial(unsigned n, unsigned k);
//#define CMDLINE_PARSER(N,...) CMDLINE_NUMCHECK##N,__VA_ARGS__);PARSECMDLINE##N __VA_ARGS__; 
// inline void findMin(const double * a, int len, double *min, int *index); // !!!!

// negation
void neg(MType * value);

// unary minus
MType umin(const MType & value);

template <typename T> T abs_value(const T &a) 
{
#if defined(USE_FLOAT)
    return fabs(a);
#elif defined(USE_INT32) || defined(USE_INT8)
    return abs(a);
#else
#error not implemented
#endif
}


/*A stack  of capacity L is implemented as an array S of size L+1.
S[0] contains the number of elements in the stack
*/
inline void StackInit(const unsigned& capacity, unsigned* pStack)
{
    pStack[0] = capacity;
    pStack[capacity] = ~0u;
}

inline void Push(const unsigned& x, unsigned* pStack) {
	pStack[++pStack[0]] = x;
};

inline unsigned Pop(unsigned* pStack) {
	assert(pStack[0] > 0);
	if (pStack[pStack[0]] == ~0u) {
		//the stack was not properly initialized yet
		unsigned s = --pStack[0];
		if (s > 0)
			pStack[pStack[0]] = ~0u;
		return s;
	} else
		return pStack[pStack[0]--];
};

inline unsigned Peek(const unsigned *pStack) {
	return pStack[pStack[0]] != ~0u ? pStack[pStack[0]] : pStack[0] - 1;
}

inline unsigned IsInit(const unsigned *pStack) {
	return pStack[pStack[0]] != ~0u;
}

inline bool IsEmpty(const unsigned *pStack) {
	return pStack[0] == 0u;
}

 void *operator new(size_t Size);
 void *operator new[](size_t Size);

 void operator delete(void* ptr);
 void operator delete[](void* ptr);

int ceilLog2(unsigned a);

unsigned binomialCoeff(unsigned n, unsigned k); 
double GetSubsetCapacity(double Sigma, unsigned SubconstellationSize, const double* pSubConstellation);
double GetBPSKAWGNChannelEntropy(double Sigma);


#define NUMOFHERMITEINTEGRATIONPOINTS 10

// Coefficients for Hermite integration
extern double Hermite_wp[NUMOFHERMITEINTEGRATIONPOINTS];
extern double Hermite_xp[NUMOFHERMITEINTEGRATIONPOINTS];

#endif