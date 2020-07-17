#pragma once
#include <map>
#include "MixedKernelEncoder.h"
#include "KernProc.h"
#include "TVMemoryEngine.h"


//#define ENABLE_SIMPLIFIED_ENGINES

//#define ENABLE_LOGGING

//#define CALL_ACCOUNTING

class CListKernelEngine :public CMixedKernelEncoder, public CTVMemoryEngine
{
protected:
	///LLR calculators
	const CKernProcLLR*** m_pppLLREngines;
//#ifdef ENABLE_SIMPLIFIED_ENGINES
/*	///simplified LLR evaluation engines 
	const CKernProcLLR** m_ppSimplifiedLLREngines;*/
	///for each branch on each layer show provide a 
	unsigned** m_ppProcessorID;
	unsigned* m_pMaxEngineID;
	/*///state variable for simplified engine
	void** m_ppSimplifiedState;
	//tDFMask* m_pDFMaskCorrection;
#endif*/
#ifdef CALL_ACCOUNTING
	typedef std::map<const CKernProcLLR*, unsigned*> tAccountMap;
	tAccountMap m_CallAccountTable;
	unsigned m_NumOfDecoderCalls;
	///unpuncture and unshorten the code
	virtual void LoadLLRs(MType* pDest ///the LLR vector as required by the decoder
		, const MType* pSrc
	)
	{
		m_NumOfDecoderCalls++;
		CMixedKernelEncoder::LoadLLRs(pDest, pSrc);
	}

#endif

	///temporary storage for LLR calculators
	void*** m_pppEngineTemp;
	///state variables for each path
	void*** m_pppEngineState;
	///correction masks for evaluation of the dynamic frozen symbols
	tDFMask* m_pDFCorrection;
	///the indices of the bits where the values of dynamic frozen symbols are accumulated
	unsigned* m_pDFValueBits;
	///codewords of outer codes
	tBit** m_ppArrayPointer_C;
	///LLRs 
	MType** m_ppArrayPointer_S;
	///path scores
	KType* m_pR;
	///lengths of outer codes
	unsigned* m_pOuterLengths;
	///allocate a new array at layer lambda
	void AllocateEntry(unsigned lambda, unsigned s1) ;
	///copy the data if the array is cloned
	unsigned UpdateIndex(unsigned lambda ///layer ID
				        ,unsigned *pIndexArray ///index array for the path
						,unsigned CurPathPhase ///current path phase
						);

public:
	CListKernelEngine(std::istream& SpecFile,unsigned MaxNumOfPaths);
	~CListKernelEngine();
	//Tal&Vardy functions
	tBit* GetArrayPointer_C(
		unsigned  lambda      //layer id
		,unsigned* pIndexArray //the pointer obtained from GetIndexArrayPointer
		,unsigned PathPhase 
	)
	{
		unsigned Index = UpdateIndex(lambda, pIndexArray,PathPhase);

		tBit* pC = m_ppArrayPointer_C[Index];
		return pC;
	}
	tBit* GetArrayPointer_C(
		unsigned lambda //layer id
		,unsigned l      //the pointer obtained from GetIndexArrayPointer
		, unsigned PathPhase ///phase of the path
	)
	{
		unsigned* pIndexArray = GetIndexArrayPointer(l);
		return GetArrayPointer_C(lambda, pIndexArray, PathPhase);
	};

	const tBit* CGetArrayPointer_C(
		unsigned lambda, //layer id
		unsigned l      //path index
	) const
	{
		unsigned *pIndexArray = GetIndexArrayPointer(l);
		return CGetArrayPointer_C(lambda, pIndexArray);
	}
	unsigned CExtractArrayIndex_C(
		unsigned lambda, //layer id
		unsigned l       //path index
	) const
	{
		unsigned *pIndexArray = GetIndexArrayPointer(l);
		unsigned Index = pIndexArray[lambda];
		m_pArrayReferenceCount[Index]++;
		return Index;
	}
	const tBit* CListKernelEngine::CGetArrayPointer_C(
		unsigned        lambda,             //layer id
		const unsigned* pIndexArrayPointer //the pointer obtained from GetIndexArrayPointer
	) const
	{
		unsigned s = pIndexArrayPointer[lambda];
		tBit* pC = m_ppArrayPointer_C[s];
		return pC;
	};
	
	///get a pointer to the state variable. This function must be called only after GetArrayPointer_C for the previous layer
	void* GetStatePointer(unsigned lambda ///layer ID
		, const unsigned* pIndexArrayPointer //the pointer obtained from GetIndexArrayPointer
		,unsigned EngineID ///ID of the engine
	)
	{
		return m_pppEngineState[EngineID][lambda - 1];
	}


	MType* GetArrayPointer_S(unsigned ArrayIndex)
	{
		return m_ppArrayPointer_S[ArrayIndex];
	}

	MType* GetArrayPointer_S(
		unsigned  lambda     //layer id
		,unsigned* pIndexArray //path index array;
		,unsigned PathPhase ///phase of the path
	)
	{
		unsigned s = UpdateIndex(lambda, pIndexArray,PathPhase);
		return m_ppArrayPointer_S[s];
	}

	MType* GetArrayPointer_S(
		unsigned lambda //layer id
		,unsigned l       //path index
		, unsigned PathPhase ///phase of the path
	)
	{
		unsigned* pIndexArray = GetIndexArrayPointer(l);
		return GetArrayPointer_S(lambda, pIndexArray,PathPhase);
	}

	const MType* CGetArrayPointer_S(
		unsigned lambda, //layer id
		unsigned l       //path index
	) const
	{
		unsigned* pIndexArray = GetIndexArrayPointer(l);
		return CGetArrayPointer_S(lambda, pIndexArray);
	}

	const MType* CGetArrayPointer_S(
		unsigned        lambda,            //layer id
		const unsigned* pIndexArrayPointer //the pointer obtained from GetIndexArrayPointer
	) const
	{
		unsigned s = pIndexArrayPointer[lambda];
		return m_ppArrayPointer_S[s];
	};


	///update the codewords of outer codes
	void IterativelyUpdateC(
		unsigned l      ///path ID
		,unsigned phi     ///phase
		, unsigned Layers2Skip=0  ///number of layers to be skipped
	);
	///compute the LLR for the input symbol of a polarizing transformation
	const MType* IterativelyCalcS(
		unsigned l      ///path ID
		,unsigned phi    ///phase
		, unsigned Layers2Skip=0 ///number of last layers to be skipped
	);

	inline const tBit * CExtractCodeword(
		unsigned l //path index
	) const
	{
		return CGetArrayPointer_C(0, l);
	}
	

	//void MakeOuterPermutationC(unsigned l, unsigned lambda, unsigned phi);
	//initialize the Tal-Vardy  data structures for specific number of layers. Only layers up to UpperLayerBoundary will be actually initialized
	//void Init(unsigned NumOfLayers, unsigned MemoryLimit);

};