#include "TVMemoryEngine.h"

// assume that instance was initialized via default constructor
void CTVMemoryEngine::InitMemoryEngine(unsigned MaxNumOfPaths)
{
    m_MaxNumOfPaths = MaxNumOfPaths;
	m_pPathPhases = new unsigned[m_MaxNumOfPaths];
	m_pInactivePathIndices = new unsigned[m_MaxNumOfPaths + 1]; //this is a stack
	m_pDFMask = new tDFMask[MaxNumOfPaths];
}

CTVMemoryEngine::CTVMemoryEngine(unsigned MaxNumOfPaths) :
m_MaxNumOfPaths(MaxNumOfPaths), m_pPathIndex2ArrayIndex(nullptr),
m_pInactiveArrayIndices(nullptr), m_pArrayReferenceCount(nullptr) {
    InitMemoryEngine(MaxNumOfPaths);
    //m_pPathPhases = new unsigned[m_MaxNumOfPaths];
    //m_pInactivePathIndices = new unsigned[m_MaxNumOfPaths + 1]; //this is a stack
    //m_pDFMask = new tDFMask[MaxNumOfPaths];
}

CTVMemoryEngine::~CTVMemoryEngine() {
#ifndef SPECIALIZED_DECODER_MODE
	delete[] m_pPathPhases;
	delete[] m_pInactivePathIndices;

	delete[] m_pDFMask;
	delete[]m_pPathIndex2ArrayIndex;
	delete[]m_pArrayReferenceCount;
	delete[]m_pInactiveArrayIndices;
#endif
}

void CTVMemoryEngine::Init(unsigned LogLength, unsigned UpperLayerBoundary, unsigned MemorySize) {
	m_Memory.Init(MemorySize);

	m_LogLength = LogLength;
	m_UpperLayer = UpperLayerBoundary;

	m_pPathIndex2ArrayIndex = new unsigned[(m_UpperLayer + 1) * m_MaxNumOfPaths];
	m_pInactiveArrayIndices = new unsigned[(m_UpperLayer + 1) * (m_MaxNumOfPaths + 1)];
	m_EmptyIndex = (m_UpperLayer + 1) * m_MaxNumOfPaths;
	m_pArrayReferenceCount = new unsigned[m_EmptyIndex + 1];
	for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) {
		for (unsigned s = 0; s < m_MaxNumOfPaths; s++)
			m_pPathIndex2ArrayIndex[s * (m_UpperLayer + 1) + lambda] = m_EmptyIndex;
	}

#ifdef MEMORY_COUNTING
	m_pBlocks = new unsigned[m_LogLength + 1];
	m_pBlocksMax = new unsigned[m_LogLength + 1];
	memset(m_pBlocks, 0, sizeof(unsigned) * (m_LogLength + 1));
	memset(m_pBlocksMax, 0, sizeof(unsigned) * (m_LogLength + 1));
#endif

	Cleanup();
}

void CTVMemoryEngine::Cleanup() {
#ifdef MEMORY_COUNTING
	for (unsigned i = 0; i < m_LogLength + 1; ++i) {
		if (m_pBlocks[i] > m_pBlocksMax[i])
			m_pBlocksMax[i] = m_pBlocks[i];
	}
#ifdef USE_MEMORY_MANAGER
	memset(m_pBlocks, 0, sizeof(unsigned) * (m_LogLength + 1));
#endif
#endif
	m_Memory.Reset();

	//final initialization
	unsigned *pInactiveArray = m_pInactiveArrayIndices + 0 * (m_MaxNumOfPaths + 1);

	for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) {
		pInactiveArray[0] = m_MaxNumOfPaths;//mark the stack as empty
		pInactiveArray[m_MaxNumOfPaths] = ~0;//force on-demand initialization of the stack
		pInactiveArray += m_MaxNumOfPaths + 1;
	}

	m_pInactivePathIndices[0] = m_MaxNumOfPaths;//mark the stack as empty
	m_pInactivePathIndices[m_MaxNumOfPaths] = ~0; //force on-demand initialization of the stack

	m_pArrayReferenceCount[m_EmptyIndex] = (~0u) / 2u;
}

unsigned CTVMemoryEngine::AssignInitialPath() {
	unsigned l = Pop(m_pInactivePathIndices);

	unsigned *pIndexArray = GetIndexArrayPointer(l);
	for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++)
		pIndexArray[lambda] = m_EmptyIndex;

	m_pDFMask[l] = 0;
	return l;
}

void CTVMemoryEngine::KillPath(unsigned l) {
	//mark the path index as inactive
	Push(l, m_pInactivePathIndices);

	//disassociate arrays with path index
	const unsigned* pArrayIndex = GetIndexArrayPointer(l);
	unsigned* pInactiveArray = m_pInactiveArrayIndices + 0 * (m_MaxNumOfPaths + 1);

	for (unsigned lambda = 0; lambda <= m_UpperLayer; lambda++) 
	{
		unsigned s = pArrayIndex[lambda];
		assert(m_pArrayReferenceCount[s] > 0);
		if (/*s != EMPTY_INDEX && */!--m_pArrayReferenceCount[s]) 
		{
			assert(s != m_EmptyIndex);
			Push(s, pInactiveArray);
		}
		pInactiveArray += (m_MaxNumOfPaths + 1);
	}
}

unsigned CTVMemoryEngine::ClonePath(unsigned l) {
	unsigned l1 = Pop(m_pInactivePathIndices);
	//make l1 reference same arrays as l
#ifdef LAYER_0_MEMORY_SAVING
	const unsigned Z = 1;
#else
	const unsigned Z = 0;
#endif
	const unsigned *pIndexArray = GetIndexArrayPointer(l);
	unsigned *pIndexArray1 = GetIndexArrayPointer(l1);

	for (unsigned lambda = Z; lambda <= m_UpperLayer; lambda++) 
	{
		unsigned s = pIndexArray[lambda];
		pIndexArray1[lambda] = s;
		//if (s != EMPTY_INDEX)
			m_pArrayReferenceCount[s]++;
	}

#ifdef PROCESS_ONLY_UNFROZEN
	memcpy(m_pLayerPhases + l1 * m_LogLength, m_pLayerPhases + l * m_LogLength, sizeof(unsigned) * m_LogLength);
#endif

	m_pDFMask[l1] = m_pDFMask[l];
	return l1;
}
