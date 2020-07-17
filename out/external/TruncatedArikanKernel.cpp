#pragma once
#include "TruncatedArikanKernel.h"
#include "LinAlg.h"
#include "SoftProcessing.h"
#include "IdentityProcessor.h"
/**

Modified Kernel:
100
110
011


*/

///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CGKernel::Multiply(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	/*XOR(pDest, pSrc, pSrc + Stride, Stride);
	XOR(pDest + Stride, pSrc, pSrc + 2 * Stride, Stride);
	XOR(pDest + 2 * Stride, pDest + Stride, pSrc + Stride, Stride);*/
	XOR(pDest, pSrc, pSrc + Stride, Stride);
	XOR(pDest + Stride, pSrc + Stride, pSrc + 2 * Stride, Stride);
	memcpy(pDest + 2 * Stride, pSrc + 2 * Stride, sizeof(tBit)*Stride);
	
}


///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^T, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CGKernel::MultiplyTransposed(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	/*XOR(pDest + Stride, pSrc, pSrc + 2 * Stride, Stride);
	XOR(pDest + 2 * Stride, pSrc + Stride, pSrc + 2 * Stride, Stride);
	XOR(pDest, pDest + Stride, pSrc + Stride, Stride);*/
	memcpy(pDest, pSrc, sizeof(tBit)*Stride);
	XOR(pDest + Stride, pSrc, pSrc + Stride, Stride);
	XOR(pDest + 2 * Stride, pSrc + Stride, pSrc + 2 * Stride, Stride);
}


///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^{-1}, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CGKernel::InverseMultiply(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	/*XOR(pDest + Stride, pSrc + Stride, pSrc + 2 * Stride, Stride);
	XOR(pDest + 2 * Stride, pSrc, pSrc + 2 * Stride, Stride);
	XOR(pDest, pDest + Stride, pSrc, Stride);*/

	memcpy(pDest + 2 * Stride, pSrc + 2 * Stride, sizeof(tBit)*Stride);
	XOR(pDest + Stride, pSrc + Stride, pSrc + 2 * Stride, Stride);
	XOR(pDest, pSrc, pDest + Stride, Stride);
}
CGKernel GKernel;
const tBit G[3 * 3] = {

	/*1, 1, 1,
	1, 0, 1,
	0, 1, 1*/
	1, 0, 0,
	1,1,0,
	0,1,1
};
const tBit* CGKernel::GetKernel()const
{
	return G;
}
///allocate an LLR computation engine
CKernProcLLR* CGKernel::GetProcessor(unsigned ProcessorID)const
{
	switch (ProcessorID)
	{
	case 0:return new CGKernelProcessor();
	case 1:	return new CIdentityProcessor(3);
	default:throw Exception("Invalid processor ID %d",ProcessorID);
	}
	
}
///get kernel name
const char* CGKernel::Name()const
{
	return "G";
};
const unsigned NumOfGPoints = 79 + 2;
double GPhysCap[NumOfGPoints] = { 0,0.00144125,0.00161692,0.00181396,0.00203499,0.0022829,0.00256096,0.00287283,0.00322258,0.00361481,0.00405465,0.00454784,0.0051008,0.00572073,0.00641567,0.0071946,0.00806759,0.00904583,0.0101419,0.0113697,0.0127448,0.0142846,0.0160083,0.0179375,0.020096,0.0225102,0.0252093,0.0282258,0.0315953,0.0353573,0.039555,0.0442359,0.0494517,0.0552591,0.0617194,0.0688989,0.0768691,0.0857063,0.0954918,0.106312,0.118255,0.131416,0.145889,0.161771,0.179157,0.19814,0.218805,0.241229,0.265479,0.2916,0.319619,0.349531,0.381298,0.414841,0.450032,0.486688,0.524572,0.563394,0.602826,0.64251,0.682077,0.721132,0.759224,0.795783,0.830095,0.861356,0.888889,0.912429,0.932287,0.94918,0.963732,0.975962,0.985271,0.9913,0.994683,0.996608,0.997968,0.99905,0.999697,0.999914,1 };
double GCapacities[NumOfGPoints * 3] = { 0,1.18854591955383e-005,1.33341451396009e-005,1.4959061621744e-005,1.67818148193085e-005,1.88262375004297e-005,2.11192961536206e-005,2.36911734540976e-005,2.65754331964321e-005,2.98100098904588e-005,3.34372087612761e-005,3.75043654798521e-005,4.20644234273038e-005,4.71767583581555e-005,5.29076731982925e-005,5.93312227082183e-005,6.65304504779411e-005,7.45976363259504e-005,8.36365228900119e-005,9.37617383628874e-005,0.000105101682813735,0.000117799847649322,0.000132014568215046,0.000147923971774479,0.000165724328182853,0.000185633348540091,0.000207891834517317,0.000232767801673148,0.000260554830127175,0.000291578661866023,0.000326195551416837,0.000364797213826825,0.00040781,0.000506799,0.000610253,0.000756405,0.000943053,0.00123752,0.0015833,0.0020914,0.00271007,0.0035995,0.00473678,0.00625239,0.00825361,0.0108616,0.0142686,0.0186309,0.0242931,0.0314866,0.0405817,0.051994,0.0661443,0.0835528,0.104775,0.130298,0.160551,0.19604,0.236944,0.283315,0.334997,0.391484,0.451912,0.515093,0.57954,0.64367,0.705592,0.763491,0.815913,0.861588,0.899816,0.930455,0.95393,0.971038,0.982786,0.990376,0.994963,0.997553,0.998901,0.999553,1,
0,0.00135706,0.00149397,0.00166538,0.00184779,0.00206226,0.00228423,0.00256247,0.00283135,0.00313461,0.00351766,0.00394657,0.00439017,0.00493225,0.00551259,0.00618433,0.00697534,0.00785856,0.00879427,0.00989343,0.0111451,0.0125329,0.0141343,0.0159458,0.0179656,0.0202841,0.022921,0.0258211,0.0291947,0.0329705,0.0372656,0.0421115,0.0476285,0.0538247,0.0608451,0.0687647,0.0777123,0.087815,0.099169,0.111892,0.126201,0.142199,0.160085,0.179985,0.202039,0.226399,0.253223,0.282571,0.314532,0.349031,0.386118,0.425631,0.467323,0.51091,0.555953,0.601901,0.648122,0.693894,0.738409,0.780885,0.820516,0.856644,0.888676,0.916214,0.939152,0.957504,0.971622,0.98194,0.989132,0.993844,0.996749,0.998418,0.999304,0.999726,0.999913,0.99998,0.999998,1,1,1,1,
0,0.00303316,0.00339764,0.0037943,0.00423372,0.00471935,0.00525758,0.00589811,0.00658087,0.00733536,0.00823135,0.00920118,0.0103004,0.0115248,0.0129121,0.0144304,0.0161657,0.0180806,0.0202412,0.0226712,0.0253597,0.0283738,0.0317311,0.0354887,0.0397008,0.0443777,0.0495915,0.0553887,0.0618421,0.0690118,0.0769914,0.0858062,0.0955898,0.106414,0.118329,0.131498,0.145952,0.16182,0.179186,0.198136,0.218812,0.241206,0.265465,0.291547,0.319561,0.349436,0.381204,0.414734,0.449908,0.486589,0.524542,0.563475,0.603049,0.642839,0.682448,0.721352,0.759038,0.795031,0.828746,0.859808,0.887803,0.912417,0.933573,0.951155,0.965327,0.976339,0.984575,0.990419,0.994386,0.996909,0.998436,0.999295,0.999724,0.999912,0.999977,0.999998,1,1,1,1,1
};

double GBEC[3 * 3] = {
	3,	3,	1,
	0,	2,	1,
	0,	1,	1
};

CAWGNCapacityInterpolator GInterp(3, NumOfGPoints, GCapacities, GPhysCap);
CBECCapacityEvaluator GBECEval(3, GBEC);

CGKernel::CGKernel():CBinaryKernel(3)
{
	memset(m_ppCapacityEvaluators, 0, sizeof(m_ppCapacityEvaluators));
	m_ppCapacityEvaluators[ec_AWGN] = &GInterp;
	m_ppCapacityEvaluators[ec_BEC] = &GBECEval;
};


///compute the capacity of a subchannel induced by the kernel for a given physical channel capacity
double CGKernel::GetSubchannelCapacity(eChannel  Channel  ///channel model to be used
	, unsigned SubchannelID ///ID of the subchannel to be analyzed
	, double PhysicalCapacity ///capacity of the underlying channel
)const
{
	if (m_ppCapacityEvaluators[Channel])
		return m_ppCapacityEvaluators[Channel]->GetSubchannelCapacity(PhysicalCapacity, SubchannelID);
	else throw Exception("G Kernel  does not support capacity evaluation for channel model %d", (unsigned)Channel);
}

///allocate state variables if necessary
void* CGKernelProcessor::GetState(unsigned Stride  ///enable processing of data blocks of size Stride
)const
{
	return new MType[Stride];
}
///deallocate state variables
void CGKernelProcessor::FreeState(void * State)const
{
	delete[]State;
}
///copy the state variable
void CGKernelProcessor::CopyState(unsigned phase ///the number of known kernel input symbols
	, unsigned Stride ///the state variable is assumed to be configured for processing of Stride LLRs
	, void* pDest ///destination state variable 
	, const void* pSrc ///source state variable
)const
{
	if (phase == 0)
		memcpy(pDest, pSrc, sizeof(MType)*Stride);
}


void CGKernelProcessor::GetLLRs(unsigned Stride ///number of LLRs to be computed 
	, unsigned phase ///local kernel phase
	, const tBit* pKnownInputSymbols ///known input symbols of the kernel. This is a vector of size Stride*m_pKernel->size
	, const MType* pChannelLLRs ///the LLRs for kernel output symbols. This is a vector of size Stride*m_pKernel->size
	, MType* pLLRs ///output: the LLRs for kernel input symbol phase. This is a vector of size Stride
	, void* pState ///state variable
	, void* pTemp  ///temporary variable
)const
{
	switch (phase)
	{
	case 0:
	{
		MType* pL1L2 = (MType*)pState;
		SoftXOR(pL1L2, pChannelLLRs + Stride, pChannelLLRs + 2 * Stride, Stride);
		SoftXOR(pLLRs, pChannelLLRs, pL1L2, Stride);
		break;
	}
	case 1:
	{
		MType* pL1L2 = (MType*)pState;
		SoftCombine(pLLRs, pChannelLLRs, pL1L2, pKnownInputSymbols, Stride);
		break;
	}
	case 2:
	{
		SoftCombine(pLLRs, pChannelLLRs + Stride, pChannelLLRs + 2 * Stride, pKnownInputSymbols + Stride, Stride);
		//ChangeLLRSigns(pLLRs, pLLRs, pKnownInputSymbols, Stride);
		break;
	}
	default: throw Exception("Invalid phase");
	};
}




/**

Shortened Arikan kernel

1 0 0 0 0
1 1 0 0 0
1 0 1 0 0
1 0 0 0 1
1 1 1 1 0


*/

///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CA5Kernel::Multiply(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	XOR(pDest, pSrc, pSrc + Stride, Stride);
	XOR(pDest+2*Stride, pSrc+2*Stride, pSrc + 4*Stride, Stride);
	XOR(pDest + Stride, pSrc + Stride, pSrc + 4 * Stride, Stride);
	XOR(pDest, pDest + 2 * Stride, Stride);
	XOR(pDest, pSrc + 3 * Stride, Stride);
	memcpy(pDest+4*Stride, pSrc + 3 * Stride, sizeof(tBit)*Stride);
	memcpy(pDest + 3 * Stride, pSrc + 4 * Stride, sizeof(tBit)*Stride);
}


///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^T, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CA5Kernel::MultiplyTransposed(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	memcpy(pDest , pSrc , sizeof(tBit)*Stride);
	XOR(pDest+Stride, pSrc, pSrc + Stride, Stride);
	XOR(pDest + 4*Stride, pSrc+2*Stride, pSrc + 3*Stride, Stride);
	XOR(pDest + 4 * Stride, pDest + Stride, Stride);
	XOR(pDest + 2 * Stride, pSrc, pSrc + 2 * Stride, Stride);
	XOR(pDest+3*Stride, pSrc,pSrc + 4 * Stride, Stride);
}


///compute (y_i,y_{i+d},y_{i+2d},...,y_{i+(l-1)d})=(x_i,x_{i+d},x_{i+2d},...,x_{i+(l-1)d})F_l^{-1}, 0\leq i<d
///where l is kernel dimension (m_Size), and d is stride
void CA5Kernel::InverseMultiply(unsigned Stride ///d
	, const tBit* pSrc ///input data (x)
	, tBit* pDest
)const
{
	XOR(pDest, pSrc, pSrc + Stride, Stride);
	XOR(pDest + 2 * Stride, pSrc + 2 * Stride, pSrc + 3 * Stride, Stride);
	XOR(pDest + Stride, pSrc + Stride, pSrc + 3 * Stride, Stride);
	XOR(pDest, pDest + 2 * Stride, Stride);
	XOR(pDest, pSrc + 4 * Stride, Stride);
	memcpy(pDest + 4 * Stride, pSrc + 3 * Stride, sizeof(tBit)*Stride);
	memcpy(pDest + 3 * Stride, pSrc + 4 * Stride, sizeof(tBit)*Stride);
}
CA5Kernel A5Kernel;
const tBit A5[5 * 5] = {

	1,0,0,0,0,
 1,1,0,0,0,
 1,0,1,0,0,
 1,1,1,1,0,
 1,0,0,0,1
};
const tBit* CA5Kernel::GetKernel()const
{
	return A5;
}
///allocate an LLR computation engine
CKernProcLLR* CA5Kernel::GetProcessor(unsigned ProcessorID)const
{
	switch (ProcessorID)
	{
	case 0:return new CA5KernelProcessor();
	case 1:	return new CIdentityProcessor(5);
	default: throw Exception("Invalid processor ID %d", ProcessorID);
	}
	
}
///get kernel name
const char* CA5Kernel::Name()const
{
	return "A5";
};
const unsigned NumOfA5Points = 81 + 2;
double A5PhysCap[NumOfA5Points] = { 0,0.00144125,0.00161692,0.00181396,0.00203499,0.0022829,0.00256096,0.00287283,0.00322258,0.00361481,0.00405465,0.00454784,0.0051008,0.00572073,0.00641567,0.0071946,0.00806759,0.00904583,0.0101419,0.0113697,0.0127448,0.0142846,0.0160083,0.0179375,0.020096,0.0225102,0.0252093,0.0282258,0.0315953,0.0353573,0.039555,0.0442359,0.0494517,0.0552591,0.0617194,0.0688989,0.0768691,0.0857063,0.0954918,0.106312,0.118255,0.131416,0.145889,0.161771,0.179157,0.19814,0.218805,0.241229,0.265479,0.2916,0.319619,0.349531,0.381298,0.414841,0.450032,0.486688,0.524572,0.563394,0.602826,0.64251,0.682077,0.721132,0.759224,0.795783,0.830095,0.861356,0.888889,0.912429,0.932287,0.94918,0.963732,0.975962,0.985271,0.9913,0.994683,0.996608,0.997968,0.99905,0.999697,0.999914,0.999962,0.999978,1 };
double A5Capacities[NumOfA5Points * 5] = { 0,3.62139835451587e-006,4.06280064345797e-006,4.55789887886044e-006,5.11327627924111e-006,5.7361944864002e-006,6.43486996009087e-006,7.21849910480752e-006,8.09730852336219e-006,9.08285653834346e-006,1.01880331921164e-005,1.14272612611285e-005,1.28166721434273e-005,1.43743571265426e-005,1.61205181482164e-005,1.80777190642844e-005,2.02712625505004e-005,2.27292654829005e-005,2.54833373610855e-005,2.85684044207036e-005,3.20235890710383e-005,3.58926119236201e-005,4.02237164118622e-005,4.50711763983545e-005,5.04947936390987e-005,5.65609028550378e-005,6.3342874267821e-005,7.09223699392153e-005,7.93888412353411e-005,8.88415389697305e-005,9.93890108675631e-005,0.00011115060917296,0.000124256239381101,0.000138848370381286,0.000155081029385399,0.000173120806999447,0.000193147359759317,0.00021535240505925,0.000239940223687604,0.000267127911094738,0.000297136834285013,0.000330206200282434,0.00036657220089642,0.000406478565972862,0.000450164,0.000637257,0.000927104,0.00142351,0.00217533,0.00336043,0.00513261,0.00776091,0.0116218,0.0171213,0.0249157,0.0357103,0.050389,0.0699638,0.0954911,0.127946,0.16811,0.216636,0.273571,0.338382,0.409827,0.486065,0.564419,0.641852,0.715465,0.782326,0.840168,0.887733,0.924928,0.952479,0.971599,0.984064,0.991655,0.995964,0.998203,0.99928,0.999745,0.999923,1,
0,4.08221370320396e-005,4.5797835080552e-005,5.13788195598533e-005,5.7639299662675e-005,6.46611320939763e-005,7.25369367240744e-005,8.13703798298383e-005,9.12767405770756e-005,0.000102386309914857,0.000114844390575514,0.000128813563004192,0.000144475668047201,0.000162034639363955,0.000181718202874134,0.000203780709169618,0.000228507382132394,0.000256215168658134,0.000287260386168425,0.000322036739922415,0.000360985236458587,0.000404598715453858,0.000453421,0.000556667,0.000666332,0.000756,0.000886399,0.00108076,0.00128136,0.00155949,0.00190499,0.00235949,0.00297196,0.00363309,0.00454582,0.00566451,0.0070964,0.00889821,0.0111035,0.0139161,0.017399,0.0217349,0.0270902,0.0337345,0.0418952,0.0519213,0.0641458,0.0788974,0.0967017,0.117955,0.143053,0.172439,0.206532,0.245624,0.28984,0.339154,0.393105,0.451093,0.512223,0.575171,0.638274,0.699865,0.75817,0.811453,0.858337,0.897884,0.929702,0.954076,0.971662,0.983618,0.991196,0.99568,0.998089,0.999241,0.999741,0.999925,0.999992,0.999998,1,1,1,1,1,
0,5.41148513153381e-005,6.07107617615241e-005,6.81090551201879e-005,7.64081049631916e-005,8.57164225969023e-005,9.61567872503232e-005,0.000107866621546743,0.00012099874244703,0.000135725866909417,0.000152240611889496,0.000170758498113407,0.000191520578379377,0.000214797192274203,0.00024089021900314,0.000270136832106388,0.000302915131533814,0.000339645270059896,0.000380799591018233,0.0004269,0.000502822,0.000580531,0.000676894,0.000842796,0.000998451,0.00123913,0.0015232,0.00189222,0.00232599,0.00286806,0.00360421,0.00446565,0.00556789,0.0069212,0.0086144,0.0107128,0.0133565,0.0166113,0.0206023,0.0255642,0.0315778,0.0389165,0.0477773,0.0585279,0.0714242,0.0867928,0.105011,0.126434,0.151428,0.180322,0.213447,0.250948,0.292905,0.339293,0.389759,0.443815,0.500632,0.559223,0.618329,0.676529,0.732491,0.784723,0.831892,0.873072,0.907697,0.935597,0.957102,0.972873,0.983809,0.990913,0.995229,0.997686,0.998963,0.999583,0.999858,0.999961,0.999997,1,1,1,1,1,1,
0,0.00129979,0.00144831,0.00158817,0.00174193,0.00192738,0.00216269,0.00242558,0.00272227,0.00302515,0.00338048,0.00381022,0.00426568,0.0048148,0.00543637,0.00612942,0.00690279,0.00775749,0.00877504,0.0099432,0.0112579,0.0127738,0.0144626,0.0164309,0.0186603,0.0212264,0.0241295,0.0274506,0.0312598,0.0356173,0.0405567,0.0462522,0.0527436,0.0600851,0.0685235,0.078075,0.0889334,0.101211,0.115138,0.13081,0.14843,0.168245,0.190266,0.214754,0.241827,0.271583,0.304017,0.339194,0.376978,0.417179,0.459597,0.503763,0.549268,0.595477,0.641763,0.687309,0.731439,0.773375,0.812459,0.848036,0.879709,0.907208,0.930428,0.949429,0.964518,0.976063,0.98456,0.990511,0.994499,0.997016,0.998478,0.999298,0.999704,0.999893,0.999974,0.999997,1,1,1,1,1,1,1,
0,0.00578175,0.00648002,0.0072416,0.0081188,0.00910687,0.0101932,0.0114098,0.0127936,0.0143142,0.0160119,0.017945,0.020091,0.0225204,0.025215,0.0282102,0.0315829,0.0353631,0.039569,0.0442348,0.0494347,0.0552507,0.0617136,0.0689024,0.0768849,0.0857295,0.0955311,0.106368,0.118328,0.131508,0.146023,0.161919,0.179357,0.198371,0.219092,0.241554,0.265878,0.292033,0.32014,0.350104,0.381938,0.415542,0.450831,0.487582,0.525581,0.56457,0.604197,0.644016,0.683596,0.722484,0.76009,0.795984,0.829615,0.860534,0.88838,0.912868,0.933891,0.951398,0.965492,0.976459,0.984641,0.990472,0.994421,0.996928,0.998438,0.999256,0.999677,0.999895,0.999969,0.999993,1,1,1,1,1,1,1,1,1,1,1,1,1,
};

double A5BEC[5 * 5] = {
	5,	10,	10,	5,	1,
	0,	6,	9,	5,	1,
	0,	3,	8,	5,	1,
	0,	1,	3,	4,	1,
	0,	0,	0,	1,	1
};

CAWGNCapacityInterpolator A5Interp(5, NumOfA5Points, A5Capacities, A5PhysCap);
CBECCapacityEvaluator A5BECEval(5, A5BEC);

CA5Kernel::CA5Kernel() :CBinaryKernel(5)
{
	memset(m_ppCapacityEvaluators, 0, sizeof(m_ppCapacityEvaluators));
	m_ppCapacityEvaluators[ec_AWGN] = &A5Interp;
	m_ppCapacityEvaluators[ec_BEC] = &A5BECEval;
};


///compute the capacity of a subchannel induced by the kernel for a given physical channel capacity
double CA5Kernel::GetSubchannelCapacity(eChannel  Channel  ///channel model to be used
	, unsigned SubchannelID ///ID of the subchannel to be analyzed
	, double PhysicalCapacity ///capacity of the underlying channel
)const
{
	if (m_ppCapacityEvaluators[Channel])
		return m_ppCapacityEvaluators[Channel]->GetSubchannelCapacity(PhysicalCapacity, SubchannelID);
	else throw Exception("G Kernel  does not support capacity evaluation for channel model %d", (unsigned)Channel);
}


///allocate state variables if necessary
void* CA5KernelProcessor::GetState(unsigned Stride  ///enable processing of data blocks of size Stride
)const
{
	return _aligned_malloc(Stride*(3 * sizeof(MType) + sizeof(tBit)),ALIGNMENT);
}
///deallocate state variables
void CA5KernelProcessor::FreeState(void * State)const
{
	_aligned_free(State);
}
///copy the state variable
void CA5KernelProcessor::CopyState(unsigned phase ///the number of known kernel input symbols
	, unsigned Stride ///the state variable is assumed to be configured for processing of Stride LLRs
	, void* pDest ///destination state variable 
	, const void* pSrc ///source state variable
)const
{
	memcpy(pDest, pSrc, Stride*(3 * sizeof(MType) + sizeof(tBit)));
}

///compute LLRs for a given kernel phase
void CA5KernelProcessor::GetLLRs(unsigned Stride ///number of LLRs to be computed 
	, unsigned phase ///local kernel phase
	, const tBit* pKnownInputSymbols ///known input symbols of the kernel. This is a vector of size Stride*m_pKernel->size
	, const MType* pChannelLLRs ///the LLRs for kernel output symbols. This is a vector of size Stride*m_pKernel->size
	, MType* pLLRs ///output: the LLRs for kernel input symbol phase. This is a vector of size Stride
	, void* pState ///state variable
	, void* pTemp  ///temporary variable
)const
{
	MType* S = (MType*)pState;
	MType* pL04 = S;
	MType* pL13 = S + Stride;
	MType* pL024 = S + 2 * Stride;
	tBit* pC = (tBit*)(S + 3 * Stride);
	switch (phase)
	{
	case 0:
	{
		SoftXOR(pL04, pChannelLLRs, pChannelLLRs + 4 * Stride, Stride);
		SoftXOR(pL024, (const MType*)pL04, pChannelLLRs + 2 * Stride, Stride);
		SoftXOR(pL13, pChannelLLRs + 1 * Stride, pChannelLLRs + 3 * Stride, Stride);
		SoftXOR(pLLRs, (const MType*)pL024, pL13, Stride);
		break;
	};
	case 1:
	{
		SoftCombine(pLLRs, pL024, pL13, pKnownInputSymbols, Stride);
		break;
	};
	case 2:
	{
		XOR(pC, pKnownInputSymbols, pKnownInputSymbols + Stride, Stride);
		SoftCombine(pL024, pL04, pChannelLLRs + 2 * Stride, pC, Stride);
		SoftCombine(pL13, pChannelLLRs + 1 * Stride, pChannelLLRs + 3 * Stride, pKnownInputSymbols+Stride, Stride);
		SoftXOR(pLLRs, (const MType*)pL024, pL13, Stride);
		break;
	}
	case 3:
	{
		//combine L1,L2,L3
		SoftCombine(pL13, pChannelLLRs + Stride, pChannelLLRs + 3 * Stride, pKnownInputSymbols + Stride, Stride);
		SoftCombine(pL13, pChannelLLRs + 2 * Stride, pL13, pKnownInputSymbols + 2 * Stride, Stride);
		SoftXOR(pL04, pChannelLLRs, pL13, Stride);
		XOR(pC, pKnownInputSymbols + 2 * Stride, Stride);
		SoftCombine(pLLRs, pL04, pChannelLLRs + 4 * Stride, pC, Stride);
		break;
	}
	case 4:
	{
		XOR(pC, pKnownInputSymbols + 3 * Stride, Stride);
		SoftCombine(pLLRs, pChannelLLRs, pL13, pC, Stride);
		break;
	}
	default:
		throw Exception("Invalid phase");
	}
}


