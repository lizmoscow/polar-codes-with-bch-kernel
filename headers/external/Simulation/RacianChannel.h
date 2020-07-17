/** \file RacianChannel.h
*  \brief Racian channel model in real domain
*
*  \author Stanislav Rets
*/

#pragma once
#include "ChannelTemplates.h"

#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>

namespace nsSimulation
{
	/**
	 * \brief Racian channel output data
	 */
	template<>
	struct CChannelData<eChannelType::RACIAN_REAL> : CChannelDataTemplate<eChannelTag::TAG_FADING>
	{
		static const eChannelType channel_type = eChannelType::RACIAN_REAL;
		double KFactor; /// K-factor, see <a href="https://en.wikipedia.org/wiki/Rice_distribution#Characterization">Shape Parameter</a> 

		CChannelData(unsigned NumOfSymbols) : CChannelDataTemplate(channel_type, NumOfSymbols), KFactor(0) {}
	};

	template<>
	class CChannel<eChannelType::RACIAN_REAL> 
		: public CChannelTemplate<eChannelType::RACIAN_REAL>
		,  public CChannelAWGNStdDev
	{
	private:
		static void gsl_error(const char * reason, const char * file,
			int line, int gsl_errno)
		{
		}

		static double Variance(double RiceOmega, void *params)
		{
			const double RiceK = *reinterpret_cast<double *>(params);

			double s2 = RiceOmega / 2 / (RiceK + 1), v2 = RiceOmega - 2 * s2;
			if (s2 == 0)
				return -1;

			double x = -v2 / 2 / s2;
			gsl_sf_result I0, I1;

			auto status = gsl_sf_bessel_I0_e(-x / 2, &I0);
			status |= gsl_sf_bessel_I1_e(-x / 2, &I1);
			if (status != GSL_SUCCESS)
				return -1;

			const double L = std::exp(x / 2) * ((1-x) * I0.val - x * I1.val);
			const double Mean2 = s2 * M_PI / 2 * L * L;

			const double Var = 2 * s2 + v2 - M_PI * s2 / 2 * L * L;
			return Var + Mean2 - 1;
		}

	public:
		double m_RiceK; /// K-factor, see <a href="https://en.wikipedia.org/wiki/Rice_distribution#Characterization">shape parameter</a> 
		double m_RiceOmega; /// <a href="https://en.wikipedia.org/wiki/Rice_distribution#Characterization">scale parameter</a> 
		double m_RiceV, m_RiceSigma;

		/**
		 * \brief Initialize racian channel
		 * \param AWGNStdDev standard deviation of the additive gaussian noise
		 * \param RiceK K-factor,the ratio of the power contributions by line-of-sight path to the remaining multipaths
		 */
		CChannel(double AWGNStdDev, double RiceK) 
			: m_RiceK(RiceK)
		{
			SetAWGNStdDev(AWGNStdDev);

			// override GSL default error handle to manually handle errors
			const auto old_handler = gsl_set_error_handler(&gsl_error);

			// We are searching for the suitable parameters for racian distribution which normalize everage signal energy for the given K-factor
			gsl_function F;
			F.function = &Variance;
			F.params = reinterpret_cast<void *>(&m_RiceK);

			gsl_root_fsolver *s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
			gsl_root_fsolver_set(s, &F, 0, 1000);

			int status;
			unsigned iter = 0;
			double r;
			do
			{
				++iter;
				status = gsl_root_fsolver_iterate(s);
				r = gsl_root_fsolver_root(s);
				double x_lo = gsl_root_fsolver_x_lower(s);
				double x_hi = gsl_root_fsolver_x_upper(s);
				status = gsl_root_test_interval(x_lo, x_hi, 0, 0.000001);
			} while (status == GSL_CONTINUE && iter < 1000);

			if (status != GSL_SUCCESS)
				throw Exception("gsl failed to find parameters for the rice distribution");

			gsl_set_error_handler(old_handler);
			m_RiceOmega = r;
			m_RiceSigma = std::sqrt(m_RiceOmega / 2 / (RiceK + 1));
			if (m_RiceOmega - 2 * m_RiceSigma * m_RiceSigma < 0)
				m_RiceV = 0;
			else
				m_RiceV = std::sqrt(m_RiceOmega - 2 * m_RiceSigma * m_RiceSigma);
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
			Data.KFactor = m_RiceK;

			for (unsigned i = 0; i < NumOfSymbols; ++i)
			{
				const unsigned P = gsl_ran_poisson(pRNG, m_RiceV * m_RiceV / 2 / m_RiceSigma / m_RiceSigma);
				const double RiceValue = m_RiceSigma * sqrt(gsl_ran_chisq(pRNG, 2 * P + 2));
				Data.pFadingFactors[i] = RiceValue;
				const Complex additive_noise = gsl_ran_gaussian(pRNG, m_AWGNStdDev);
				Data.pSymbols[i] = Data.pFadingFactors[i] * pIn[i] + additive_noise;
			}
		}

		/**
		* \brief Check if ml-error accured, returns true if decoded sequence is closer to the received data than the correct one
		* \param Data data received from the channel output
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
		* \brief Allocate container for the channel data. Should be manually deallocated with delete operator.
		* \param NumOfSymbols number of symbols in the transmitted block
		* \return
		*/
		CChannelDataInterface * allocate_data(unsigned NumOfSymbols = 1) override { return new channel_data_type(NumOfSymbols); }
	};
}
