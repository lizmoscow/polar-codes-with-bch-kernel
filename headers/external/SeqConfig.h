#ifndef SEQCONFIG_H
#define SEQCONFIG_H

#define DECODER_TYPE 0  // Original SeqConfig.h for compatibility
//#define DECODER_TYPE 1  // SEQUENTIAL DECODER
//#define DECODER_TYPE 2  // BLOCK SEQUENTIAL DECODER
//#define DECODER_TYPE 3  // BLOCK SEQUENTIAL DECODER X
//#define DECODER_TYPE 4  // BLIND SEQUENTIAL DECODER X
//#define DECODER_TYPE 17 // MULTILEVEL SEQUENTIAL DECODER

#if defined(_SEQUENTIAL_DECODER_)
#undef DECODER_TYPE
#define DECODER_TYPE 1
#endif
#if defined(_BLOCK_SEQUENTIAL_DECODER_)
#undef DECODER_TYPE
#define DECODER_TYPE 2
#endif
#if defined(_BLOCK_SEQUENTIAL_DECODER_X_)
#undef DECODER_TYPE
#define DECODER_TYPE 3
#endif
#if defined(_BLIND_SEQUENTIAL_DECODER_X_)
#undef DECODER_TYPE
#define DECODER_TYPE 4
#endif
#if defined(_MULTILEVEL_SEQUENTIAL_DECODER_)
#undef DECODER_TYPE
#define DECODER_TYPE 17
#endif

#if   DECODER_TYPE == 1
#include "SeqDecConfig.h"
#elif DECODER_TYPE == 2
#include "BlockSeqConfig.h"
#elif DECODER_TYPE == 3 
#include "BSDXConfig.h"
#elif DECODER_TYPE == 4 
#include "BlindSeqConfig.h"
#elif DECODER_TYPE == 0
#include "SeqConfigOrig.h"
#elif DECODER_TYPE == 17
#include "MLSeqConfig.h"
#else
#error No configuration file was specified
#endif

#endif // SEQCONFIG_H