/** \file Channel.h
*  \brief File with channel types and basic interfaces
*
*  This contains implememnted channel types and channel tags. 
*  The latter define interfaces for different channel families like additive channels and multipath ones.
*
*  \author Stanislav Rets
*/

#pragma once

#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

#include "misc.h"

namespace nsSimulation {
	
	/**
	 * \brief Supported channel global tags
	 */
	enum class eChannelTag
	{
		TAG_ADDITIVE,        /// Additive channel tag
		TAG_FADING,          /// Fading channel tag
		TAG_MULTIPATH_FADING /// Multipath channel tag
	};

	/**
	 * \brief List of existing channel implementations
	 */
	enum class eChannelType
	{
		AWGN_REAL,	              /// AWGN channel, real
		RAYLEIGH_REAL,	          /// Rayleigh fading channel, real
		RAYLEIGH_BLOCK_REAL,      /// Block rayleigh fading channel, real
		RAYLEIGH_CORRELATED_REAL, /// Correlated Rayleigh fading channel, real
		RACIAN_REAL,              /// Racian fading channel, real
	};


	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	/**
	* \brief This template should be used to implement new types of channel. 
	* Each specialization must inherit CChannelInterface or CChannelWithDataInterface
	* \tparam Type Channel type listed in #eChannelType
	*/
	template<eChannelType Type>
	class CChannel;

	/**
	 * \brief Specialization of this template contains data returned by a channel.
	 * Each specilization must inherit #CChannelDataInterface
	 * \tparam Type Channel type listed in #eChannelType
	 */
	template<eChannelType Type>
	struct CChannelData;


	/**
	 * \brief Templates which provide basic channel interfaces given by a tag
	 * \tparam Tag channel tag which specifies the channel interface
	 */
	template<eChannelTag Tag>
	struct CChannelDataTemplate;


	template<eChannelType Tag>
	class CChannelTemplate;

	/**
	 * \brief Interface for cast operations for a channel and a channel data
	 * \tparam TypeParamType Type of the template variable of the given \t Type
	 * \tparam Type Type with one template varaible
	 * \tparam HasArbitraryTypeCast True if arbitrary type casts are required
	 */
	template<typename TypeParamType, template<TypeParamType Value> typename Type, bool HasArbitraryTypeCast = true>
	struct CCastTemplate
	{
		virtual ~CCastTemplate() {}

		// @{
		/* 
		 * \brief Safe cast to \t SomeType. These casts are disabled if \p X == false using SFINAE 
		 * \tparam SomeType 
		 * \tparam E 
		 * \return 
		 */
		template<typename SomeType, typename = typename std::enable_if<HasArbitraryTypeCast>::type> SomeType * cast() { return dynamic_cast<SomeType *>(this); }
		template<typename SomeType, typename = typename std::enable_if<HasArbitraryTypeCast>::type> const SomeType * castc() const { return dynamic_cast<const SomeType *>(this); }
		// @}

		// @{
		/*
		* \brief Unsafe \v this pointer cast into the pointer to \tparam Type with the given template argument \tparam Value
		*/
		template<TypeParamType Value> Type<Value> * cast() { return dynamic_cast<Type<Value> *>(this); }
		template<TypeParamType Value> const Type<Value> * castc() const { return dynamic_cast<const Type<Value> *>(this); }
		// @}

		// @{
		/*
		* \brief Same as cast/castc functions but a bit faster. If the cast does not exist it leads to UB
		*/
		template<typename SomeType, typename = typename std::enable_if<HasArbitraryTypeCast>::type> SomeType * unsafe_cast() { return static_cast<SomeType *>(this); }
		template<typename SomeType, typename = typename std::enable_if<HasArbitraryTypeCast>::type> const SomeType * unsafe_castc() const { return static_cast<const SomeType *>(this); }
		template<TypeParamType Value> Type<Value> * unsafe_cast() { return static_cast<Type<Value> *>(this); }
		template<TypeParamType Value> const Type<Value> * unsafe_castc() const { return static_cast<const Type<Value> *>(this); }
		// @}
	};

	/**
	* \brief Base class for the channel output. If possible your class should iherit one of #CChannelDataTemplate insteed of this one
	*/
	struct CChannelDataInterface : CCastTemplate<eChannelTag, CChannelDataTemplate>, CCastTemplate<eChannelType, CChannelData, false>
	{
	private:
		const eChannelType m_Type; /// channel type
		const eChannelTag  m_Tag;  /// channel tag
	public:
		const unsigned NumOfSymbols; /// number of symbols per transmitted block
		tChannelSymbol * const pSymbols; /// array of \p NumOfSymbols symbols

		CChannelDataInterface(eChannelType Type, eChannelTag Tag, unsigned NumOfSymbols) 
		: m_Type(Type), m_Tag(Tag), NumOfSymbols(NumOfSymbols), pSymbols(new tChannelSymbol[NumOfSymbols]) {}
		
		virtual ~CChannelDataInterface() { delete[] pSymbols; }
		eChannelType type() const { return m_Type; }
		eChannelTag tag() const { return m_Tag; }
	};

	/**
	 * \brief Channel interface, its implementations should be listened in #tChannelType
	 */
	class CChannelInterface : public CCastTemplate<eChannelType, CChannel>
	{
	private:
		const eChannelType m_Type; /// channel type
		const eChannelTag  m_Tag;  /// channel tag
	public:

		CChannelInterface(eChannelType type, eChannelTag tag) : m_Type(type), m_Tag(tag) {}

		virtual ~CChannelInterface() {}

		/**
		 * \brief Allocate container for the channel data. Should be deallocated with delete operator.
		 * \param NumOfSymbols number of symbols in the transmitted block
		 * \return pointer to the allocated data
		 */
		virtual CChannelDataInterface * allocate_data(unsigned NumOfSymbols = 1) = 0;

		/**
		 * \brief Check if ml-error accured, returns true if decoded sequence is closer to the received data than the correct one
		 * \param Data data received from the channel output
		 * \param pCorrectArray pointer to the array of transmitted symbols
		 * \param pDecodedArray pointer to the array of decoded symbols
		 * \param NumOfSymbols number of symbols in both arrays
		 * \return True if maximum likelihood error occured, false - otherwise
		 */
		virtual bool check_ml_error(const CChannelDataInterface &Data
			, const tChannelSymbol *pCorrectArray
			, const tChannelSymbol *pDecodedArray, unsigned NumOfSymbols) = 0;

		/**
		 * \brief Transmit data through the channel
		 * \param pRNG pointer to the rng source
		 * \param pIn pointer to the array of symbols to be transmitted
		 * \param NumOfSymbols number of symbols in the array
		 * \param OutData  channel output, should be allocated by allocate_data(NumOfSymbols)
		 */
		virtual void apply(gsl_rng *pRNG, const tChannelSymbol *pIn, unsigned NumOfSymbols, CChannelDataInterface &OutData) = 0;


		/**
		 * \brief Get channel type
		 */
		eChannelType type() const { return m_Type; }


		/**
		 * \brief Get channel tag
		 */
		eChannelTag tag() const { return m_Tag; }

		/**
		 * \brief Get name of the channel. Default: the class name
		 */
		virtual std::string name() const { return typeid(*this).name(); }
	};

	template<eChannelType Type>
	class CChannelTemplate : public CChannelInterface
	{
	public:
		typedef CChannelData<Type> channel_data_type;
		static const eChannelType channel_type = Type;
		static const eChannelTag channel_tag = channel_data_type::channel_tag;

		CChannelTemplate() : CChannelInterface(channel_type, channel_tag) {}
	};

	namespace transform
	{
		/**
		 * \brief Decompose one QAM symbol and the given fading factor into two PAM symbols with corresponding fading factors.
		 *  Energy of the resulting symbols is normailzed.
		 * \param QAMSymbol Symbol of QAM
		 * \param QAMFading Fading factor of the symbol
		 * \param pSymbolsOut Pointer to the resulting array of two PAM symbols
		 * \param pFadingOut Pointer to the resulting array of two fading factors of the OAM symbols
		 */
		inline void QAM_into_PAM(const tChannelSymbol &QAMSymbol, const Complex &QAMFading, 
			tChannelSymbol *pSymbolsOut, Complex *pFadingOut)
		{
			double abs_h = abs(QAMFading);
			Complex y = QAMSymbol * std::conj(QAMFading / abs_h) * M_SQRT2;
			if (pSymbolsOut) 
			{
				pSymbolsOut[0] = y.real();
				pSymbolsOut[1] = y.imag();
			}
			if (pFadingOut)
			{
				pFadingOut[0] = abs_h;
				pFadingOut[1] = abs_h;
			}
		}

		/**
		 * \brief Merge two PAM symbols into one QAM and normalize their energy
		 * \param pPAMSymbols pointer to the array of 2 PAM symbols
		 * \param QAMSymbolOut reference to the result
		 */
		inline void PAM_into_QAM(const tChannelSymbol *pPAMSymbols, tChannelSymbol &QAMSymbolOut)
		{
			QAMSymbolOut = tChannelSymbol(pPAMSymbols[0].real(), pPAMSymbols[1].real()) / M_SQRT2;
		}
	}
}
