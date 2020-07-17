#pragma once

#include <memory>

#include "SeqConfig.h"
#include "misc.h"
#include "MemoryAllocator.h"

typedef tBit TV_TYPE_C;
typedef MType TV_TYPE_S;
typedef tPairOfProbabilities TV_TYPE_P;

typedef unsigned long long tDFMask;

class DLL_API CTVMemoryEngine
{

protected:
#ifdef USE_MEMORY_MANAGER
    CBlockMemoryAllocator
#else
    CSimpleMemoryAllocator
#endif 
        m_Memory;
    unsigned m_EmptyIndex;
    unsigned m_UpperLayer;

    /// This contains m+1 stacks of size L
    unsigned* m_pInactiveArrayIndices;

public:

	tDFMask*  m_pDFMask;
    unsigned* m_pArrayReferenceCount;
    /// Stack of size L
    unsigned* m_pInactivePathIndices;
    unsigned* m_pPathIndex2ArrayIndex;
    unsigned  m_MaxNumOfPaths;
    /// The phase of each path
    unsigned* m_pPathPhases;
    unsigned* m_pBlocks;
    unsigned* m_pBlocksMax;
    unsigned  m_LogLength;

    CTVMemoryEngine() {};
    CTVMemoryEngine(unsigned MaxNumOfPaths);
    virtual ~CTVMemoryEngine();

    /// Initialize data structures for given max num of paths
    void InitMemoryEngine(unsigned MaxNumOfPaths);

    /// Initialize the Tal-Vardy  data structures for specific number of layers. 
    /// Only layers up to UpperLayerBoundary will be actually initialized
    void Init(unsigned LogLength, unsigned UpperLayerBoundary, unsigned MemorySize);

    /// cleanup all datastructures
    void Cleanup();

	// return pointer to pathIndex2ArrayIndex for path number l (pointer, because there is an index for each layer)
    inline unsigned * GetIndexArrayPointer(unsigned l) const
    {
        return m_pPathIndex2ArrayIndex + l * (m_UpperLayer + 1);
    }

	// return path index l by pointer to pathIndex2ArrayIndex (calculate difference between pointers)
    inline unsigned GetPathIndex(const unsigned *pArrayIndex) const
    {
        return static_cast<unsigned>((pArrayIndex - m_pPathIndex2ArrayIndex) / (m_UpperLayer + 1));
    }

	// check if pointer to layer lambda of this path is used more than one time
	// if so, then pIndexArray[lambda] is assigned to a new unused index for this layer
    inline unsigned UpdateIndex(unsigned lambda, unsigned *pIndexArray)
    {
        unsigned Index = pIndexArray[lambda];

        if (m_pArrayReferenceCount[Index] > 1) {
            m_pArrayReferenceCount[Index]--;
            Index = GetNewArrayIndex(lambda);
            m_pArrayReferenceCount[Index] = 1;
            pIndexArray[lambda] = Index;
        }

        return Index;
    }

    /// Get an index of an unused array at layer lambda
    unsigned GetNewArrayIndex(unsigned lambda)
    {
        unsigned* pStack = m_pInactiveArrayIndices + lambda*(m_MaxNumOfPaths + 1);
        assert(pStack[0] != ~0u);
        if (!pStack[0])
            //out of memory
            longjmp(JUMP_BUFFER, 1);
        if (pStack[pStack[0]] == ~0u) {
            //the stack was not properly initialized yet. The memory has not been allocated too
            unsigned s = --pStack[0];
            if (s > 0u)
                pStack[pStack[0]] = ~0u;
            s = s * (m_UpperLayer + 1) + lambda;
            AllocateEntry(lambda, s);
            return s;
        } else
            return pStack[pStack[0]--];
    }

	// allocate memory for new C array at layer lambda for path with index s1
    virtual void AllocateEntry(unsigned lambda, unsigned s1) = 0;

    unsigned AssignInitialPath();
    virtual void KillPath(unsigned l);
    unsigned ClonePath(unsigned l);

    unsigned GetMaxNumOfPaths() const
    {
        return m_MaxNumOfPaths;
    }

    unsigned GetNumOfUnusedPaths() const
    {
        return m_pInactivePathIndices[0];
    }
};