/** \file NextSimulator.h
*  \brief Monte-Carlo simulation
*
*  This contains implememntation of Monte-Carlo simulation algorithm for the given modulated decoder-encoder pair and the channel.
*  It also contains adapters for the old interfaces based on #CBinaryCodec.
*
*  \author Stanislav Rets
*/

#pragma once

#include <gsl/gsl_rng.h>
#include <chrono>
#include <memory>
#include <iterator>

#include "Codec.h"
#include "Simulation/ChannelTemplates.h"
#include "Modem.h"
#include "NextSimulator.h"


namespace nsSimulation {
	/**
	 * \brief Transform CBinaryEncoder with given CModem into the CModulatedEncoder
	 */
	class DLL_API CModulatedEncoderAdapter : public CModulatedEncoder
	{
	protected:
		CBinaryEncoder &m_Encoder;
		CModem &m_Modem;
		/// temporary codeword array for the #m_Encoder
		tBit *m_pBinaryCodeword; 

	public:
		/**
		 * \brief Construct the adapter with the given binary encoder and modem
		 * \param Encoder Reference to the encoder
		 * \param Modem Reference to the modem
		 */
		CModulatedEncoderAdapter(CBinaryEncoder &Encoder, CModem &Modem);

		virtual ~CModulatedEncoderAdapter();

		/**
		 * \brief Encode the data and obtain the corresponding channel symbols
		 * \param pSrc pointer to the array of #m_Encoder.GetDimension() bits
		 * \param pEncoded pointer to the output modulated codeword of length std::ceil(#m_Encoder.GetLength() / #m_Modem.GetBitsPerSymbol())
		 */
		virtual void Encode(const tBit *pSrc, tChannelSymbol *pEncoded) override;
	};

	/**
	 * \brief Transform CBinarySoftDecoder with given CModem into the CModulatedDecoder
	 */
	class DLL_API CModulatedDecoderAdapter : public CModulatedDecoder
	{
	protected:
		CBinarySoftDecoder &m_Decoder;
		CModem &m_Modem;
		/// Temporary array to store LLRs returned from the #m_Modem 
		MType *m_pLLRs;
		/// Temporary array for the codeword list returned by the #m_Decoder.Decode()
		tBit *m_pBinaryCodewordList;
		/// true if the the last bits of the codeword does עשו fully fill up a channel symbol
		const bool m_ArrayRoundingRequired;
		/// Temporary array for the last bits of the codeword
		tBit *m_pRoundingArray;

		/**
		 * \brief Modulate the binary codeword
		 * \param pBinaryCodeword pointer to the binary codeword
		 * \param pCodeword pointer to the codeword array of #GetLength() channel symbols 
		 */
		void ModulateBinaryCodeword(const tBit *pBinaryCodeword, tChannelSymbol *pCodeword) const;
	public:

		/**
		 * \brief Construct the adapter with the given binary decoder and modem
		 * \param Decoder reference to the decoder
		 * \param Modem reference to the modem
		 */
		CModulatedDecoderAdapter(CBinarySoftDecoder &Decoder, CModem &Modem);

		virtual ~CModulatedDecoderAdapter();
		
		/**
		* \brief Soft-input modulated decoder. Produces a list of information vectors
		* \param pChannelData   channel output with noised symbols and additional data
		* \param pInfVectorList output array. Must be sufficiently big to accomodate GetMaxListSize() vectors of length GetDimension()
		* \param pCodewordList output array. If non-zero, must be sufficiently big to accomodate GetMaxListSize() vectors of length GetLength()
		* \return the number of vectors obtained, or -1 in case of error
		*/
		virtual int Decode(
			const CChannelDataInterface *pChannelData
			, tBit *pInfVectorList, tChannelSymbol *pCodewordList) override;

		/**
		* \brief Adjust the decoder to the current channel
		* \param Channel current channel
		*/
		virtual void SetChannel(nsSimulation::CChannelInterface &Channel) override
		{
			auto AWGNStdDev = dynamic_cast<nsSimulation::CChannelAWGNStdDev &>(Channel).GetAWGNStdDev();
			m_Decoder.SetAWGNStdDev(AWGNStdDev);
		}
	};

	class CNextSimulator;

	/**
	* \brief this function will be called after each simulation iteration. If it returns false, simulations will be terminated
	*/
	typedef  bool(*tSimulatorHook)(CNextSimulator *pEngine ///pointer to the simulator engine
		, bool DecodingFailed ///true if the last decoding attempt has failed
		);

	/**
	 * \brief Simulator for the given modulated encoder and decoder
	 */
	class CNextSimulator
	{
	public:
		/// is used for the Ctrl+C handling
		static bool KeepRunning;
		/// pointer to the random number generator
		gsl_rng* m_pRNG;

		/// pointer to the current modulated encoder
		CModulatedEncoder *m_pEncoder;
		/// pointer to the current modulated decoder
		CModulatedDecoder *m_pDecoder;
		/// true if adapters were allocated in the constructor
		bool m_IsAdapterUsed;
		/// the channel to be used
		CChannelInterface &m_Channel;
		/// pointer to the data object allocated by #m_Channel.allocate_data()
		CChannelDataInterface *m_pChannelData;

		/// the data being transmitted
		tBit *m_pData;
		/// the encoded symbols
		tChannelSymbol *m_pEncoded;
		/// decoded information symbols
		tBit *m_pDecodedData;
		/// decoded codeword
		tChannelSymbol *m_pDecodedCW;

		/// number of bit errors in the decoded data
		unsigned long long m_NumOfBitErrors;
		/// number of codeword errors
		unsigned long long m_NumOfBlockErrors;
		/// number of most-likelihood errors
		unsigned long long m_NumOfMLErrors;
		/// number of undetected errors
		unsigned long long m_NumOfUndetectedError;
		/// number of iterations already done
		unsigned long long m_IterationsDone;
		/// overall decoding time
		unsigned long long m_DecodingTime;

		unsigned long long m_MemoryUsage;
		unsigned long long m_MemoryUsageSuccess;
		unsigned long long m_MemoryUsageStatic;

		
		/**
		 * \brief Run one simulation iteration
		 * \param It iteration ID
		 * \param Indeed true if the data should be indeed coded and decoded
		 * \param pData do this for some specific data
		 * \return true if decoding is correct
		 */
		bool Iterate(unsigned long long It, bool Indeed, const tBit* pData = nullptr);
		
		/**
		 * \brief Initialize simulator by given modulated coder, decoder and channel.
		 * Note that the simulator overrides a random generator of the channel.
		 * \param pEncoder Pointer to a modulated encoder
		 * \param pDecoder Pointer to a modulated decoder
		 * \param Channel Reference to a channel model
		 */
		CNextSimulator(CModulatedEncoder *pEncoder, CModulatedDecoder *pDecoder, CChannelInterface &Channel);
		

		/**
		 * \brief Constructor for backward compatibility with old interfaces.
		 * Initialize simulator by given binary codec, modem and channel model.
		 * Note that the simulator overrides the random generator of the channel.
		 * \param Codec The codec to be used. Must represent both binary coder and decoder
		 * \param Modem The modem to be used
		 * \param Channel The Reference to the channel model
		 */
		CNextSimulator(CBinaryCodec &Codec, CModem &Modem, CChannelInterface &Channel);

		virtual ~CNextSimulator();

		/**
		 * \brief Reset the current counters in the simulator
		 */
		void ResetCounters();

		
		/**
		 * \brief run simulations 
		 * \param MaxIterations stop after this number of iterations
		 * \param MaxErrors stop after this number of errors
		 * \param Hook the function to be called after each iteration 
		 */
		virtual void Simulate(unsigned long long MaxIterations, unsigned long long MaxErrors, tSimulatorHook Hook = DefaultHook);


		/**
		 * \brief Advance to and do the real encoding and decoding only for the specified iteration
		 * \param IterationID iteration id
		 * \param pData do this for some specific data pattern. can be #nullptr
		 */
		void RunOneIteration(unsigned IterationID, const tBit* pData = 0);

		/**
		 * \brief Get current bit error reate
		 */
		virtual double BER() const
		{
			return double(m_NumOfBitErrors) / (m_IterationsDone*m_pEncoder->GetDimension());
		}

		/**
		 * \brief Get current frame error rate
		 */
		virtual double FER() const
		{
			return double(m_NumOfBlockErrors) / m_IterationsDone;
		}

		/**
		 * \brief Get most-likelihood frame error rate
		 */
		virtual double MLFER() const
		{
			return double(m_NumOfMLErrors) / m_IterationsDone;
		}

		/**
		* \brief Get undetected frame error rate
		*/
		virtual double UndetectedFER() const
		{
			return double(m_NumOfUndetectedError) / (m_IterationsDone);
		}

		/**
		 * \brief Get average decoding time
		 */
		double GetAvgDecodingTime() const;

		/**
		 * \brief Set seed for the random number generator
		 */
		void SetRngSeed(unsigned long Seed = std::chrono::system_clock::now().time_since_epoch().count()) const;

		/**
		 * \brief Print two data array into the stdout and compare them
		 * \tparam T data type
		 * \param pData0 pointer to the first array
		 * \param pData1 pointer to the second array
		 * \param N size of arrays
		 */
		template<typename T>
		static void printData(const T* pData0, const T* pData1, unsigned N)
		{
			std::copy_n(pData0, N, std::ostream_iterator<T>(std::cout, ""));
			std::cout << std::endl;
			std::copy_n(pData1, N, std::ostream_iterator<T>(std::cout, ""));
			std::cout << std::endl;
			const bool match = std::equal(pData0, pData0 + N, pData1);
			std::cout << "equals? " << (match ? "YES" : "NO") << std::endl;
		}

		/**
		 * \brief Get the reference to the current channel
		 */
		const CChannelInterface & Channel() const
		{
			return m_Channel;
		}

		/// set it to true if erroneous iteration numbers should be printed
		bool m_PrintErrorIterations;
		/// set it to true if zero data should be used
		bool m_ZeroData;
		/// set it to true and data will contain no zeroes
		bool m_OneData;
		/// set it to true to force assertion failure on non-ML decoding error
		bool m_AssertNonMLError;
		/// set it to true if BER is not needed
		bool m_NoBER;
		/// set it to true if decoded data should be printed to the stdout
		bool m_PrintDecodedData;
		/// set it to true if decoded codeword should be printed to the stdout
		bool m_PrintCodeword;

		static bool EmptyHook(CNextSimulator *pEngine, bool Correct)
		{
			return true;
		}

		static bool DefaultHook(CNextSimulator* pSim, bool Correct) {
			if (!(pSim->m_IterationsDone % 1000))
				std::cout << pSim->m_IterationsDone << ' ' << std::scientific << pSim->FER() << std::endl;
			return true;
		}
	};
}
