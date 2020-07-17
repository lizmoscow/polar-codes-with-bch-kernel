/** \file RayleighBlockChannel.h
*  \brief Rayleigh block fading channel model in real domain
*
*  \author Stanislav Rets
*/

#pragma once

#include "Channel.h"

namespace nsSimulation
{
	template<>
	struct CChannelData<eChannelType::RAYLEIGH_BLOCK_REAL> : CChannelDataTemplate<eChannelTag::TAG_FADING>
	{
		static const eChannelType channel_type = eChannelType::RAYLEIGH_BLOCK_REAL;

		const unsigned BlockSize; /// number of symbols per block

		/**
		* \brief initialize block fading channel output
		* \param NumOfBlocks number of blocks per symbol array
		* \param BlockSize number of symbols per block
		*/
		CChannelData(unsigned NumOfBlocks, unsigned BlockSize)
			: CChannelDataTemplate(channel_type, NumOfBlocks * BlockSize), BlockSize(BlockSize) {}
	};

	template<>
	class CChannel<eChannelType::RAYLEIGH_BLOCK_REAL> 
		: public CChannelTemplate<eChannelType::RAYLEIGH_BLOCK_REAL>
		, public CChannelAWGNStdDev
	{
	private:
		const unsigned m_BlockSize; /// number of symbols per block
	public:

		/**
		* \brief initialize block fading channel
		* \param AWGNStdDev standard deviation of AWGN
		* \param BlockSize number of symbols per block
		*/
		explicit CChannel(double AWGNStdDev, unsigned BlockSize)
			: m_BlockSize(BlockSize)
		{
			SetAWGNStdDev(AWGNStdDev);
		}

		/**
		* \brief Transmit data through the channel
		* \param pRNG pointer to the rng source
		* \param pIn pointer to the array of symbols to be transmitted
		* \param NumOfSymbols number of symbols in the array
		* \param OutData  channel output, should be allocated by allocate_data(NumOfSymbols)
		*/
		virtual void apply(gsl_rng *pRNG, const tChannelSymbol *pIn, unsigned NumOfSymbols, CChannelDataInterface &OutData) override
		{
			channel_data_type &Data = dynamic_cast<channel_data_type &>(OutData);
			Data.SetAWGNStdDev(*this);

			const unsigned NumOfBlocks = NumOfSymbols / m_BlockSize + (NumOfSymbols % m_BlockSize > 0);

			for (unsigned blk = 0; blk < NumOfBlocks; ++blk) 
				std::fill_n(Data.pFadingFactors + blk * m_BlockSize, m_BlockSize, gsl_ran_rayleigh(pRNG, M_SQRT1_2));

			for (unsigned SymbolID = 0; SymbolID < NumOfSymbols; ++SymbolID)
			{
				const Complex additive_noise = gsl_ran_gaussian(pRNG, m_AWGNStdDev);
				Data.pSymbols[SymbolID] = Data.pFadingFactors[SymbolID] * pIn[SymbolID] + additive_noise;
			}
		}

		/**
		* \brief Check if ml-error accured, returns true if decoded sequence is closer to the received data than the correct one
		* \param ChannelData data received from the channel output
		* \param pCorrectArray pointer to the array of transmitted symbols
		* \param pDecodedArray pointer to the array of decoded symbols
		* \param NumOfSymbols number of symbols in both arrays
		* \return True if maximum likelihood error occured, false - otherwise
		*/
		bool check_ml_error(
			const CChannelDataInterface &ChannelData
			, const tChannelSymbol *pCorrectArray
			, const tChannelSymbol *pDecodedArray
			, unsigned NumOfSymbols) override
		{
			const channel_data_type &Data = dynamic_cast<const channel_data_type &>(ChannelData);

			// compare the euclidean distance
			double MLCorr = 0, CurCorr = 0;
			for (unsigned i = 0; i < NumOfSymbols; i++)
			{
				MLCorr -= std::norm(pCorrectArray[i] * Data.pFadingFactors[i] - Data.pSymbols[i]);
				CurCorr -= std::norm(pDecodedArray[i] * Data.pFadingFactors[i] - Data.pSymbols[i]);
			}

			return CurCorr >= MLCorr;
		}

		/**
		* \brief Allocate container for the channel data. Should be deallocated with delete operator.
		* \param NumOfSymbols number of symbols in the transmitted block
		* \return pointer to the allocated data
		*/
		CChannelDataInterface * allocate_data(unsigned NumOfSymbols) override
		{
			const unsigned NumOfBlocks = NumOfSymbols / m_BlockSize + (NumOfSymbols % m_BlockSize > 0);
			return new channel_data_type(NumOfBlocks, m_BlockSize);
		}
	};
}