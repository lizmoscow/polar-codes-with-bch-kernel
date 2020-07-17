#pragma once

#include "misc.h"

// Memory bulk for aligned array blocks allocation
struct CBlockMemoryAllocator {
	unsigned char *pBulk, *pBulkStart, *pBulkPos;
	const unsigned char *pBulkEnd;
	unsigned m_MemoryLimit;

	CBlockMemoryAllocator() :
		pBulk(nullptr), pBulkStart(nullptr), pBulkPos(nullptr), pBulkEnd(nullptr),m_MemoryLimit(0) {
	}

	~CBlockMemoryAllocator() {
#ifndef SPECIALIZED_DECODER_MODE
		_aligned_free(pBulk);
#endif
	}

	inline void Init(unsigned MemoryLimit) {
		pBulkStart = pBulk = static_cast<unsigned char *>(_aligned_malloc(MemoryLimit, ALIGN));
		pBulkEnd = pBulk + MemoryLimit;
		m_MemoryLimit = MemoryLimit;
		Reset();
	}

	inline void InitStatic(unsigned MemoryLimit, unsigned char * pBulk) {
		this->pBulk = pBulkStart = pBulk;
		pBulkEnd = pBulk + MemoryLimit;
		m_MemoryLimit = MemoryLimit;
		Reset();
	}

	template<typename T>
	inline T * Allocate(const unsigned &Size) {
		unsigned char *pOldPos = pBulkPos;
		size_t Mask = ALIGN - 1u;
		pBulkPos = reinterpret_cast<unsigned char *>((reinterpret_cast<size_t>(pBulkPos)+Mask) & ~Mask);
		size_t BlockSizeInBytes = sizeof(T) * Size;
		T *pRes = reinterpret_cast<T *>(pBulkPos);
		pBulkPos += BlockSizeInBytes;
		if (pBulkPos >= pBulkEnd) {
			CNT(CNT_MEMORY_ERROR, 1);
			longjmp(JUMP_BUFFER, 1);
		}

		CNT(CNT_MEMORY_CS, reinterpret_cast<size_t>(pBulkPos)-reinterpret_cast<size_t>(pOldPos));

		return pRes;
	}

	template<class T>
	inline void Free(T *pArray) {
	}

	inline void Reset() {
		pBulkPos = pBulkStart;
		CNT_RESET(CNT_MEMORY_CS);
	}
};

struct CSimpleMemoryAllocator {
	inline void Init(unsigned MemoryLimit = 0) {
		Reset();
	}

	inline void InitStatic(unsigned MemoryLimit, unsigned char * pBulk) {
	}

	inline void Reset() {
	}

	template<typename T>
	inline T * Allocate(const unsigned &Size) {
		unsigned long long BlockSizeInBytes = Size * sizeof(T);
		CNT(CNT_MEMORY_CS, BlockSizeInBytes);
		return static_cast<T *>(_aligned_malloc(BlockSizeInBytes, ALIGN));
	}

	

	template<typename T>
	inline void Free(T *pArray) {
		_aligned_free(pArray);
	}

	template<typename T>
	inline void FreeStatic(T *pArray) {
		Free(pArray);
	}
};