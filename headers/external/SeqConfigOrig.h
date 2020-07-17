#ifndef _SEQCONFIGORIG_H_
#define _SEQCONFIGORIG_H_

#include <cstdint>

//#define CODE_GENERATOR_MODE
//#define SPECIALIZED_DECODER_MODE

#if defined(CODE_GENERATOR_MODE) && defined(SPECIALIZED_DECODER_MODE)
#error only code generation or specialized mode allowed in the same time
#endif

//#define CHECK_HARD_DECISIONS
//#define AVOID_COPYING

#define USE_FLOAT
//#define USE_INT32
//#define USE_INT8

#define AVX
#define IMMEDIATE_DF_EVALUATION
// #define AVX2
//#define USE_MEMORY_MANAGER

//#define EARLY_TERMINATION

//#define SCORE_STATISTICS
//#define SCORE_MEAN
//#define LLR_MEAN
//#define BITERRPROB_MEASURE
//#define BRANCH_STATISTICS
//#define FOLLOW_CORRECT_PATH

//#define NO_HEURISTIC

//add to a codeword a random vector after encoding
//#define RANDOMIZE_CODEWORD
#define USE_PARITY_BF

#ifdef AVX2
#define AVX
#endif

#ifndef __AVX__
#undef AVX
#endif

//#ifndef __AVX2__
//#undef AVX2
//#endif

#define USE_MBDS
#define MBDS_L 5

#define BEST_PATH_CACHING

// output phase at which the correct path is killed
//#define TRACK_CORRECT_PATH

// directive, that allows decoder not to exceed
// the maximal number of iterations by ajusting 
// the MaxBranchCount parameter
// L = max number of branches, that can be decreased
//#define ADAPTIVE_MAX_BRANCH_COUNT 

//#define ADAPTIVE_KILL // for increasing L parameter with killed paths restoring

//#define DISTRIBUTION_KILL // to get statistics
#ifdef ADAPTIVE_KILL
	#define SAVED_PATH_NUMBER_BOUND 100000
	#define GROUP_KILL_BOUND 20
	#define MAX_L_INCREASE_COEFF 4
	#define MAX_L_BOUND 2048
	#define SAVE_GAP 15
#endif

//define this to permute codeword symbols after encoding and before decoding
//#define REVERSE_BITS

//#define ZERO_CODEWORD
#define MINSUM_HEURISTIC


#define SYSTEMATIC_ENCODING 

#define OPERATION_COUNTING
//#define MEMORY_COUNTING

#define USE_S

//#define  PACKED_DF_EVALUATION

//#define TRACK_ZERO_CODEWORD
//#define TRACK_LIST_POSITION

//#define SUBTREE_KILLER

//#define REFCOUNTLOGGING

//#define OUTER_DECODER_USAGE_LOG

//#define DECODER_TRACE
#ifdef EXTERNAL_DEFINITIONS
#undef USE_FLOAT
#undef MINSUM_HEURISTIC
#undef SYSTEMATIC_ENCODING

#ifdef EXTERNAL_USE_FLOAT
#define USE_FLOAT
#endif
#ifdef EXTERNAL_MIN_SUM_HEURISTIC
#define MINSUM_HEURISTIC
#endif
#ifdef EXTERNAL_SYSTEMATIC_ENCODING
#define SYSTEMATIC_ENCODING
#endif
#endif


#ifdef TRACK_LIST_POSITION
#define TRACK_ZERO_CODEWORD
#endif

//#define FORCE_ZERO_CODEWORD
#if defined(SCORE_STATISTICS)||defined(SCORE_MEAN)||defined(BITERRPROB_MEASURE)||defined(LLR_MEAN)
#define FORCE_ZERO_CODEWORD
#define ZERO_CODEWORD
#endif

#ifdef TRACK_ZERO_CODEWORD
#define ZERO_CODEWORD
#endif


#ifdef USE_S

//#define HEURISTIC_QUALITY_STATISTIC
#ifdef HEURISTIC_QUALITY_STATISTIC 
#define ZERO_CODEWORD
#endif

#else
#ifdef USE_MBDS
#error MBDS is not compatible with log-probability domain decoding
#endif
#endif

// Visual Studio configuration
#ifdef _USE_FLOAT_
#undef USE_FLOAT
#undef USE_INT32
#undef USE_INT8
#define USE_FLOAT
#endif
#ifdef _USE_INT8_
#undef USE_FLOAT
#undef USE_INT32
#undef USE_INT8
#define USE_INT8
#endif

typedef unsigned char tBit;
/**  Metric type, discretization
* and upper bound (MBDS mapping): 
* T' = UPPER_BOUND - T, where
* T is a classic path metric.
* BIT_1 - one bit in codeword
* BIT_0 - zero bit in codeword
*/
#if defined(USE_FLOAT)
typedef float MType;
typedef float KType;
#define MTYPE_UPPER_BOUND 100000
#define KTYPE_UPPER_BOUND 100000
#define BIT_0 0x00
#define BIT_1 0x80
#define BOOL2BIT( b ) b << 7
#elif defined(USE_INT32)
typedef int32_t MType;
typedef int32_t KType;
#define NORMALIZE_MODEM
#define MTYPE_UPPER_BOUND 0x3FFFFFFF
#define KTYPE_UPPER_BOUND 0x3FFFFFFF
#define BIT_0 0x00
#define BIT_1 0x01
#elif defined(USE_INT8)
typedef int8_t MType;
typedef int32_t KType;  // TODO: can use int16?
#define NORMALIZE_MODEM
#define MTYPE_UPPER_BOUND 0x3F
#define KTYPE_UPPER_BOUND 0x3FFFFFFF
#define BIT_0 0x00
#define BIT_1 0x80
#define BOOL2BIT( b ) b << 7
#else
#error not implemented
#endif

#define ALIGN 32

// Guard
#if (defined(USE_FLOAT) && (defined(USE_INT32) || defined(USE_INT8))) \
    || (defined(USE_INT32) && defined(USE_INT8))
#error one and only one USE_x directive should be defined
#endif


#endif  // _SEQCONFIGORIG_H_