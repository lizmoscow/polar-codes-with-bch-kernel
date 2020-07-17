/** \file RayleighCorrelatedChannel.h
*  \brief Rayleigh correlated fading channel
*  
*  The covariance matrix is given by \f$ a_{i,j}=\rho^{|i-j|} \f$. The real domain output is obtained by transforming the complex domain one.
*
*  \author Stanislav Rets
*/

#pragma once

#include "ChannelTemplates.h"
#include <gsl/gsl_linalg.h>

namespace nsSimulation
{
	template<>
	struct CChannelData<eChannelType::RAYLEIGH_CORRELATED_REAL> : CChannelDataTemplate<eChannelTag::TAG_FADING>
	{
		static const eChannelType channel_type = eChannelType::RAYLEIGH_CORRELATED_REAL;

		double Rho;
		const double *pCovarienceMatrix;

		/**
		 * \brief 
		 * \param NumOfSymbols number of symbols transmitted through the channel
		 */
		CChannelData(unsigned NumOfSymbols) : CChannelDataTemplate(channel_type, NumOfSymbols), Rho(0),
											pCovarienceMatrix(nullptr)
		{
		}
	};

	template<>
	class CChannel<eChannelType::RAYLEIGH_CORRELATED_REAL> 
		: public CChannelTemplate<eChannelType::RAYLEIGH_CORRELATED_REAL>
		, public CChannelAWGNStdDev
	{
	private:
		double m_Rho; /// covariance factor
		double *m_pCovarianceMatrix; /// covariance matrix
		gsl_matrix *m_pCholeskyFactor; /// cholesky decomposition of the covariance matrix
		gsl_vector *m_pMean; /// mean values of the transmitted symbols, should be zero.
		gsl_vector *m_pGeneratedFadingFactors1, *m_pGeneratedFadingFactors2; /// generated fading factors for symbols of QAM
		unsigned m_NumOfSymbols; /// number of real domain symbols transmitted through the channel
		unsigned m_NumOfSymbolsQAM; /// corresponding number of QAM-symbls = m_NumOfSymbols / 2
	public:

		/**
		* \brief Initialize rayleigh correlated fading channel
		* \param AWGNStdDev standard deviation of AWGN
		* \param Rho correlation coefficient. the covariance matrix is given by \f$ a_{i,j}=\rho^{|i-j|} \f$
		* \param NumOfSymbols number of symbols transmitted through the channel
		*/
		CChannel(double AWGNStdDev, double Rho, unsigned NumOfSymbols)
		{
			SetAWGNStdDev(AWGNStdDev);
			m_Rho = Rho;
			m_NumOfSymbols = NumOfSymbols;
			m_NumOfSymbolsQAM = m_NumOfSymbols / 2 + (m_NumOfSymbols % 2 > 0);

			m_pCovarianceMatrix = new double[m_NumOfSymbolsQAM * m_NumOfSymbolsQAM];
			for (unsigned i = 0; i < m_NumOfSymbolsQAM; ++i)
			{
				for (unsigned j = 0; j < m_NumOfSymbolsQAM; ++j)
					m_pCovarianceMatrix[i * m_NumOfSymbolsQAM + j] = std::pow(m_Rho, std::abs(int(i) - int(j)));
			}

			m_pCholeskyFactor = gsl_matrix_alloc(m_NumOfSymbolsQAM, m_NumOfSymbolsQAM);
			for (unsigned i = 0; i < m_NumOfSymbolsQAM; ++i)
				std::copy_n(m_pCovarianceMatrix + i * m_NumOfSymbolsQAM, m_NumOfSymbolsQAM, m_pCholeskyFactor->data + i * m_pCholeskyFactor->tda);
			m_pMean = gsl_vector_alloc(m_NumOfSymbolsQAM);
			m_pGeneratedFadingFactors1 = gsl_vector_alloc(m_NumOfSymbolsQAM);
			m_pGeneratedFadingFactors2 = gsl_vector_alloc(m_NumOfSymbolsQAM);

			gsl_vector_set_zero(m_pMean);
			auto status = gsl_linalg_cholesky_decomp(m_pCholeskyFactor);
		}

		~CChannel()
		{
			gsl_matrix_free(m_pCholeskyFactor);
			gsl_vector_free(m_pMean);
			gsl_vector_free(m_pGeneratedFadingFactors1);
			gsl_vector_free(m_pGeneratedFadingFactors2);

			delete[] m_pCovarianceMatrix;
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
			assert(m_NumOfSymbols == NumOfSymbols);

			channel_data_type &Data = dynamic_cast<channel_data_type &>(OutData);
			Data.SetAWGNStdDev(*this);
			Data.Rho = m_Rho;

			// Generate correlated fading factors for QAM
			auto status = gsl_ran_multivariate_gaussian(pRNG, m_pMean, m_pCholeskyFactor, m_pGeneratedFadingFactors1);
			status |= gsl_ran_multivariate_gaussian(pRNG, m_pMean, m_pCholeskyFactor, m_pGeneratedFadingFactors2);

			for (unsigned i = 0; i < m_NumOfSymbolsQAM; ++i)
			{
				const Complex additive_noise = Complex(
					gsl_ran_gaussian(pRNG, m_AWGNStdDev / M_SQRT2)
					, gsl_ran_gaussian(pRNG, m_AWGNStdDev / M_SQRT2)
				);
				tChannelSymbol In;
				transform::PAM_into_QAM(pIn + 2 * i, In);

				const Complex Fading = Complex(
					gsl_vector_get(m_pGeneratedFadingFactors1, i) / M_SQRT2
					, gsl_vector_get(m_pGeneratedFadingFactors2, i) / M_SQRT2
				);

				const Complex Noisy = In * Fading + additive_noise;

				Complex NoisyPAM[2], FadingPAM[2];
				// Split QAM into two PAM
				transform::QAM_into_PAM(Noisy, Fading, NoisyPAM, FadingPAM);

				Data.pSymbols[2 * i] = NoisyPAM[0];
				Data.pFadingFactors[2 * i] = FadingPAM[0];

				if (NumOfSymbols % 2 == 0 || i != m_NumOfSymbolsQAM - 1) {
					Data.pSymbols[2 * i + 1] = NoisyPAM[1];
					Data.pFadingFactors[2 * i + 1] = FadingPAM[1];
				}
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
			pData->pCovarienceMatrix = m_pCovarianceMatrix;
			return pData;
		}
	};
}