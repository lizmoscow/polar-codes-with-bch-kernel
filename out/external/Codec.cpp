#include <algorithm>
#include "Codec.h"


using namespace std;

//construct a generator matrix in packed form
//pGenMatrix must have size >=m_CodeDimension*ceil(GetLength/(8*sizeof(Word)) if GetExtension()==1
//Otherwise, >=m_CodeDimension*GetLength
Word* CBinaryEncoder::GetGeneratorMatrix()
{
    tBit* pEncodingBuffer=new tBit[GetDimension()];
    memset(pEncodingBuffer,0,sizeof(tBit)*GetDimension());
    tBit* pTemp=new tBit[GetLength()];
    unsigned RowSize=GetLength()/(8*sizeof(Word))+((GetLength()%(8*sizeof(Word)))?1:0);
    Word* pGenMatrix=new Word[RowSize*GetDimension()];
    memset(pGenMatrix,0,sizeof(Word)*RowSize*GetDimension());
    Word* pDest=pGenMatrix;
    for(unsigned i=0;i<GetDimension();i++,pDest+=RowSize)
    {
        pEncodingBuffer[i]=1;
        Encode(pEncodingBuffer,pTemp);
        for(unsigned j=0;j<GetLength();j++)
            if (pTemp[j])
                pDest[j/(8*sizeof(Word))]|=ONE<<(j%(8*sizeof(Word)));
        pEncodingBuffer[i]=0;
    };
    delete[]pTemp;
    delete[]pEncodingBuffer;
    return pGenMatrix;
};

//constructs check matrix in packed form 
//pCheckMatrix must have size >=m_CodeLength*ceil((m_CodeLength-m_CodeDimension)/(8*sizeof(Word)) 
Word* CBlockCodec::GetCheckMatrix()
{
    //get the code generator matrix
    unsigned RowSize=m_Length/(8*sizeof(Word))+((m_Length%(8*sizeof(Word)))?1:0);
    Word* pGenMatrix=GetGeneratorMatrix();
    unsigned* pColumnPermutation=new unsigned[m_Length];
    for(unsigned i=0;i<m_Length;i++)
        pColumnPermutation[i]=i;
    //perform Gaussian elimination
    Word* pCurRow=pGenMatrix;
    for (unsigned i=0;i<m_Dimension;i++,pCurRow+=RowSize)
    {
        unsigned Column;
        for (Column=i;Column<m_Length;Column++)
        {
            bool Success=false;
            //look for a row having 1 at this position
            if (!(pCurRow[pColumnPermutation[Column]/(8*sizeof(Word))]&(ONE<<(pColumnPermutation[Column]%(8*sizeof(Word))))))
            {
                //try other rows
                Word* pSearchRow=pCurRow+RowSize;
                for (unsigned j=i+1;j<m_Dimension;j++,pSearchRow+=RowSize)
                {
                    if (pSearchRow[pColumnPermutation[Column]/(8*sizeof(Word))]&(ONE<<(pColumnPermutation[Column]%(8*sizeof(Word)))))
                    {
                        //add this row to the original one
                        Success=true;
                        for (unsigned k=0;k<RowSize;k++)
                            pCurRow[k]^=pSearchRow[k];
                        break;
                    }
                };
                if (Success)
                    //exchange columns
                    swap(pColumnPermutation[i],pColumnPermutation[Column]);

            } else
            {
                swap(pColumnPermutation[i],pColumnPermutation[Column]);
                break;
            };
        };
        if (Column==m_Length)
            throw exception("Generator matrix is not full rank!");
        //eliminate 1's outsize the diagonal
        Word* pRow=pGenMatrix;
        for (unsigned j=0;j<m_Dimension;j++,pRow+=RowSize)
        {
            if (j==i) continue;
            if (pRow[pColumnPermutation[i]/(8*sizeof(Word))]&(ONE<<(pColumnPermutation[i]%(8*sizeof(Word)))))
            {
                //add this row to the original one
                for (unsigned k=0;k<RowSize;k++)
                    pRow[k]^=pCurRow[k];
            };
        };
    };
    unsigned Redundancy=m_Length-m_Dimension;
    unsigned ColumnSize=Redundancy/(8*sizeof(Word))+((Redundancy%(8*sizeof(Word)))?1:0);
    Word* pCheckMatrix=new Word[m_Length*ColumnSize];
    memset(pCheckMatrix,0,sizeof(Word)*m_Length*ColumnSize);
    for(unsigned i=0;i<Redundancy;i++)
        pCheckMatrix[pColumnPermutation[m_Dimension+i]*ColumnSize+i/(8*sizeof(Word))]=ONE<<(i%(8*sizeof(Word)));
    for (unsigned i=0;i<m_Dimension;i++)
    {
        for (unsigned j=0;j<Redundancy;j++)
        {
            if (pGenMatrix[i*RowSize+pColumnPermutation[m_Dimension+j]/(8*sizeof(Word))]&(ONE<<(pColumnPermutation[m_Dimension+j]%(8*sizeof(Word)))))
                pCheckMatrix[pColumnPermutation[i]*ColumnSize+j/(8*sizeof(Word))]|=ONE<<(j%(8*sizeof(Word)));
        }
    };
    delete[]pColumnPermutation;
    delete[]pGenMatrix;
    return pCheckMatrix;
};



CBinaryLinearEncoder::CBinaryLinearEncoder(unsigned Length, unsigned Dimension, const tBit* pGenMatrix, unsigned MinDist) :m_MinDist(MinDist), CBinaryEncoder(Length,Dimension)
{
	m_pGeneratorMatrix=new tBit[Dimension*Length];
	memcpy(m_pGeneratorMatrix,pGenMatrix,Dimension*Length*sizeof(tBit));
	m_Dimension=Dimension;
	m_Length=Length;
	m_pInformationSet=new unsigned [GetLength()];
	for(unsigned i=0;i<GetLength();i++)
		m_pInformationSet[i]=i;
	unsigned RowSize0=GetLength()/(8*sizeof(Word))+((GetLength()%(8*sizeof(Word)))?1:0);
	unsigned RowSize1=GetDimension()/(8*sizeof(Word))+((GetDimension()%(8*sizeof(Word)))?1:0);
	//get the generator matrix in the binary format
	Word* pGeneratorMatrix=GetGeneratorMatrix();
	m_pMessageExtractionMatrix=new Word[GetDimension()*RowSize1];
	memset(m_pMessageExtractionMatrix,0,sizeof(Word)*GetDimension()*RowSize1);
	for(unsigned i=0;i<GetDimension();i++)
		m_pMessageExtractionMatrix[i*RowSize1+(i/(8*sizeof(Word)))]|=ONE<<(i%(8*sizeof(Word)));

	//perform Gaussian elimination
	Word* pCurRow0=pGeneratorMatrix;
	Word* pCurRow1=m_pMessageExtractionMatrix;
	for (unsigned i=0;i<GetDimension();i++,pCurRow0+=RowSize0,pCurRow1+=RowSize1)
	{
		unsigned Column;
		for (Column=i;Column<m_Length;Column++)
		{
			bool Success=false;
			//look for a row having 1 at this position
			if (!(pCurRow0[m_pInformationSet[Column]/(8*sizeof(Word))]&(ONE<<(m_pInformationSet[Column]%(8*sizeof(Word))))))
			{
				//try other rows
				Word* pSearchRow0=pCurRow0+RowSize0;
				Word* pSearchRow1=pCurRow1+RowSize1;
				for (unsigned j=i+1;j<GetDimension();j++,pSearchRow0+=RowSize0,pSearchRow1+=RowSize1)
				{
					if (pSearchRow0[m_pInformationSet[Column]/(8*sizeof(Word))]&(ONE<<(m_pInformationSet[Column]%(8*sizeof(Word)))))
					{
						//add this row to the original one
						Success=true;
						for (unsigned k=0;k<RowSize0;k++)
						{
							pCurRow0[k]^=pSearchRow0[k];
						};
						for (unsigned k=0;k<RowSize1;k++)
						{
							pCurRow1[k]^=pSearchRow1[k];
						};
						break;
					}
				};
				if (Success)
				{
					//exchange columns
					unsigned Temp=m_pInformationSet[i];
					m_pInformationSet[i]=m_pInformationSet[Column];
					m_pInformationSet[Column]=Temp;break;
				}

			} else
			{
				unsigned Temp=m_pInformationSet[i];
				m_pInformationSet[i]=m_pInformationSet[Column];
				m_pInformationSet[Column]=Temp;
				break;
			};
		};
		if (Column==GetLength())
		{
			throw exception("Generator matrix is singular");//the matrix is not full rank
		};
		//eliminate 1's below the diagonal
		Word* pRow0=pCurRow0+RowSize0;
		Word* pRow1=pCurRow1+RowSize1;
		for (unsigned j=i+1;j<GetDimension();j++,pRow0+=RowSize0,pRow1+=RowSize1)
		{
			if (pRow0[m_pInformationSet[i]/(8*sizeof(Word))]&(ONE<<(m_pInformationSet[i]%(8*sizeof(Word)))))
            {
                //add this row to the original one
				for (unsigned k=0;k<RowSize0;k++)
                    pRow0[k]^=pCurRow0[k];
				for (unsigned k=0;k<RowSize1;k++)
                    pRow1[k]^=pCurRow1[k];
           };
        };
    };
    //reverse pass
	pCurRow0=pGeneratorMatrix+RowSize0;
	pCurRow1=m_pMessageExtractionMatrix+RowSize1;
	for (unsigned i=1;i<GetDimension();i++,pCurRow0+=RowSize0,pCurRow1+=RowSize1)
    {
		Word* pDestRow0=pGeneratorMatrix;
		Word* pDestRow1=m_pMessageExtractionMatrix;
		for (unsigned j=0;j<i;j++,pDestRow0+=RowSize0,pDestRow1+=RowSize1)
        {
			if (pDestRow0[m_pInformationSet[i]/(8*sizeof(Word))]&(ONE<<(m_pInformationSet[i]%(8*sizeof(Word)))))
            {
				for (unsigned k=0;k<RowSize0;k++)
                    pDestRow0[k]^=pCurRow0[k];
				for (unsigned k=0;k<RowSize1;k++)
                    pDestRow1[k]^=pCurRow1[k];
            }
        }
    }
	delete[]pGeneratorMatrix;

};

//given a codeword, extract information bits from it
void CBinaryLinearEncoder::ExtractInformationBits(const tBit* pCodeword,tBit* pInfBits)
{
	memset(pInfBits,0,sizeof(tBit)*GetDimension());
	const Word* pRow=m_pMessageExtractionMatrix;
	unsigned RowSize=GetDimension()/(8*sizeof(Word))+((GetDimension()%(8*sizeof(Word)))?1:0);
	for(unsigned i=0;i<GetDimension();i++,pRow+=RowSize)
	{
		if (pCodeword[m_pInformationSet[i]])
		{
			for(unsigned j=0;j<GetDimension();j++)
				if (pRow[j/(8*sizeof(Word))]&(ONE<<(j%(8*sizeof(Word)))))
					pInfBits[j]^=1;
		};
	};
};



CBinaryLinearEncoder::~CBinaryLinearEncoder()
{
	delete[]m_pInformationSet;
	delete[]m_pMessageExtractionMatrix;
	delete[]m_pGeneratorMatrix;
};
//encoder of a binary linear block code given by a generator matrix
//encode m_CodeDimension bits; return m_CodeLength on success
void CBinaryLinearEncoder::Encode(const tBit* pSrc,tBit* pEncoded)
{
	for(unsigned j=0;j<m_Length;j++)
	{
		tBit R=0;
		for(unsigned i=0;i<m_Dimension;i++)
		{
			if (pSrc[i]) R^=m_pGeneratorMatrix[i*m_Length+j];
		};
		pEncoded[j]=R;
	};
};



