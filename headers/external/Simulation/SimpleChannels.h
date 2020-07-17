/** \file SimpleChannels.h
*  \brief Some simple channel models in real domain.
*         AWGN channel, Rayleigh fading channel
*
*  The covariance matrix is given by \f$ a_{i,j}=\rho^{|i-j|} \f$. The real domain output is obtained by transforming the complex domain one.
*
*  \author Stanislav Rets
*/

#pragma once

#include "ChannelTemplates.h"

namespace nsSimulation 
{
	template<>
	struct CChannelData<eChannelType::AWGN_REAL> : CChannelDataTemplate<eChannelTag::TAG_ADDITIVE>
	{
		static const eChannelType channel_type = eChannelType::AWGN_REAL;

		/**
		 * \brief 
		 * \param NumOfSymbols number of symbols transmitted through the channel
		 */
		CChannelData(unsigned NumOfSymbols) : CChannelDataTemplate(channel_type, NumOfSymbols) {}
	};

	template<>
	class CChannel<eChannelType::AWGN_REAL> : public CChannelTemplate<eChannelType::AWGN_REAL>
											, public CChannelAWGNStdDev
	{
	public:
		
		/**
		 * \brief initialize the channel
		 * \param AWGNStdDev standard deviation of AWGN
		 */
		CChannel(double AWGNStdDev)
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

			for (unsigned i = 0; i < NumOfSymbols; ++i) {
				const Complex additive_noise = gsl_ran_gaussian(pRNG, m_AWGNStdDev);
				OutData.pSymbols[i] = pIn[i] + additive_noise;
			}
		}

		/**
		* \brief Allocate container for the channel data. Should be deallocated with delete operator.
		* \param NumOfSymbols number of symbols in the transmitted block
		* \return pointer to the allocated data
		*/
		virtual CChannelDataInterface * allocate_data(unsigned NumOfSymbols = 1) override
		{
			auto pData = new CChannelData<channel_type>(NumOfSymbols);
			return pData;
		}

		/**
		* \brief Check if ml-error accured, returns true if decoded sequence is closer to the received data than the correct one
		* \param Data data received from the channel output
		* \param pCorrectArray pointer to the array of transmitted symbols
		* \param pDecodedArray pointer to the array of decoded symbols
		* \param NumOfSymbols number of symbols in both arrays
		* \return True if maximum likelihood error occured, false - otherwise
		*/
		virtual bool check_ml_error(
			const CChannelDataInterface &Data
			, const tChannelSymbol *pCorrectArray
			, const tChannelSymbol *pDecodedArray
			, unsigned NumOfSymbols) override
		{
			// compare the euclidean distance
			double MLCorr = 0, CurCorr = 0;
			for (unsigned i = 0; i < NumOfSymbols; i++)
			{
				MLCorr -= std::norm(pCorrectArray[i] - Data.pSymbols[i]);
				CurCorr -= std::norm(pDecodedArray[i] - Data.pSymbols[i]);
			}

			return CurCorr >= MLCorr;
		}
	};

	template<>
	struct CChannelData<eChannelType::RAYLEIGH_REAL> : CChannelDataTemplate<eChannelTag::TAG_FADING>
	{
		static const eChannelType channel_type = eChannelType::RAYLEIGH_REAL;
		
		/**
		* \brief
		* \param NumOfSymbols number of symbols transmitted through the channel
		*/
		CChannelData(unsigned NumOfSymbols)
			: CChannelDataTemplate(channel_type, NumOfSymbols) {}
		~CChannelData() {}
	};

	template<>
	class CChannel<eChannelType::RAYLEIGH_REAL> 
		: public CChannelTemplate<eChannelType::RAYLEIGH_REAL>
		, public CChannelAWGNStdDev
	{
	public:

		/**
		* \brief initialize the channel
		* \param AWGNStdDev standard deviation of AWGN
		*/
		CChannel(double AWGNStdDev) { SetAWGNStdDev(AWGNStdDev); }

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

			for (unsigned i = 0; i < NumOfSymbols; ++i)
			{
				Data.pFadingFactors[i] = gsl_ran_rayleigh(pRNG, M_SQRT1_2);
				const Complex additive_noise = gsl_ran_gaussian(pRNG, m_AWGNStdDev);
				Data.pSymbols[i] = Data.pFadingFactors[i] * pIn[i] + additive_noise;
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
		CChannelDataInterface * allocate_data(unsigned NumOfSymbols = 1) override
		{
			auto pData = new channel_data_type(NumOfSymbols);
			return pData;
		}
	};
}