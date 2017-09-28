#include <stdio.h>
#include <string.h>
#include "ippi.h"
#include "ipps.h"
#include "ippcore.h"

// This is a sample code for ipp wavelet for image


void func_wavelet()
{
	IppiWTFwdSpec_32f_C1R* pSpecFwd;
	IppiWTInvSpec_32f_C1R* pSpecInv;
	int specSizeFwd, specSizeInv;

	static const int initLength = 8;

	Ipp32f pSrc[initLength * initLength] = { 
		0.0,  0.0,  0.0, 11, 11,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0, 11, 11,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0, 11, 11,  0.0,  0.0,  0.0,
		11, 11, 11, 11, 11, 11, 11, 11,
		11, 11, 11, 11, 11, 11, 11, 11,
		0.0,  0.0,  0.0, 11, 11,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0, 11, 11,  0.0,  0.0,  0.0,
		0.0,  0.0,  0.0, 11, 11,  0.0,  0.0,  0.0 };
	for (int i = 0; i < initLength * initLength; i++)
	{
		if (i % initLength == 0)
			printf("\n");
		printf("%-*.2f", 10, pSrc[i]);
	}

	// DWT Forward

	int srcStep = initLength * sizeof(Ipp32f);
	IppiSize roiSize = { initLength, initLength };


	// forward, filter length

	static const Ipp32f pFwdFltL[4] =
	{ -1.294095225509215e-001f,     2.241438680418574e-001f,     8.365163037374690e-001f,     4.829629131446903e-001f };
	static const Ipp32f pFwdFltH[4] =
	{ -4.829629131446903e-001f,     8.365163037374690e-001f,    -2.241438680418574e-001f,    -1.294095225509215e-001f };
	
	static const int fwdFltLenL = 4;
	static const int fwdFltLenH = 4;

	static const int fwdFltOffsH = 1;
	static const int fwdFltOffsL = 1;

	static const int add = 4;
	static const int addoffset = 2;

	// forward, added length
	static const int totalWidthWithBorders = initLength + add;
	static const int totalHeightWithBorders = initLength + add;
	Ipp32f pSrcB[totalWidthWithBorders * totalHeightWithBorders];
	int srcStepB = totalWidthWithBorders * sizeof(Ipp32f);
	IppiSize roiSizeB = { totalWidthWithBorders, totalHeightWithBorders };

	// forward, dst
	Ipp32f pDetailXDst[4 * 4];
	Ipp32f pDetailYDst[4 * 4];
	Ipp32f pDetailXYDst[4 * 4];
	Ipp32f pApproxDst[4 * 4];
	IppiSize dstRoiSize = { 4, 4 };
	Ipp32f pDstInv[8 * 8];

	//adds border to the source image
	ippiCopyWrapBorder_32s_C1R((Ipp32s*)pSrc, srcStep, roiSize, (Ipp32s*)pSrcB, srcStepB, roiSizeB, addoffset, addoffset);

	printf("\n\nAfter Adding Border, pSrcB: \n");
	for (int i = 0; i < totalWidthWithBorders * totalHeightWithBorders; i++)
	{
		if (i % totalWidthWithBorders == 0)
			printf("\n");
		printf("%-*.2f", 10, pSrcB[i]);
	}

	//performs forward  wavelet transform 
	int bufSizeFwd, bufSizeInv;
	Ipp8u* pBufferFwd;
	Ipp8u* pBufferInv;

	ippiWTFwdGetSize_32f(1, fwdFltLenL, fwdFltOffsL, fwdFltLenH, fwdFltOffsH, &specSizeFwd, &bufSizeFwd);

	pSpecFwd = (IppiWTFwdSpec_32f_C1R*)ippMalloc(specSizeFwd);
	pBufferFwd = (Ipp8u*)ippMalloc(bufSizeFwd);

	ippiWTFwdInit_32f_C1R(pSpecFwd, pFwdFltL, fwdFltLenL, fwdFltOffsL, pFwdFltH, fwdFltLenH, fwdFltOffsH);

	int approxStep, detailXStep, detailYStep, detailXYStep;
	approxStep = detailXStep = detailYStep = detailXYStep = 4 * sizeof(Ipp32f);
	ippiWTFwd_32f_C1R(
		pSrcB + addoffset * roiSizeB.width + addoffset, srcStepB,
		pApproxDst, approxStep,
		pDetailXDst, detailXStep,
		pDetailYDst, detailYStep,
		pDetailXYDst, detailXYStep,
		dstRoiSize, pSpecFwd, pBufferFwd);


	// DWT Inverse

	static const Ipp32f pInvFltL[4] =
	{ 4.829629131446903e-001f,     8.365163037374690e-001f,     2.241438680418574e-001f,    -1.294095225509215e-001f };
	static const Ipp32f pInvFltH[4] =
	{ -1.294095225509215e-001f,    -2.241438680418574e-001f,     8.365163037374690e-001f,    -4.829629131446903e-001f };

	static const int invFltLenL = 4;
	static const int invFltLenH = 4;

	static const int invFltOffsH = 2; 
	static const int invFltOffsL = 2; 

	Ipp32f pAppB[10 * 10];
	Ipp32f pXB[10 * 10];
	Ipp32f pYB[10 * 10];
	Ipp32f pXYB[10 * 10];

	// inverse, init length
	int stepDstInv = initLength * sizeof(Ipp32f);

	ippiWTInvGetSize_32f(1, invFltLenL, invFltOffsL, invFltLenH, invFltOffsH, &specSizeInv, &bufSizeInv);
	pSpecInv = (IppiWTInvSpec_32f_C1R*)ippMalloc(specSizeInv);
	pBufferInv = (Ipp8u*)ippMalloc(bufSizeInv);
	ippiWTInvInit_32f_C1R(pSpecInv, pInvFltL, invFltLenL, invFltOffsL, pInvFltH, invFltLenH, invFltOffsH);

	//adds border to four images obtained after ippiWTFwd    
	ippiCopyWrapBorder_32s_C1R((Ipp32s*)pApproxDst, approxStep, dstRoiSize, (Ipp32s*)pAppB, 10 * sizeof(Ipp32f), { 10,10 }, 3, 3);
	ippiCopyWrapBorder_32s_C1R((Ipp32s*)pDetailXDst, detailXStep, dstRoiSize, (Ipp32s*)pXB, 10 * sizeof(Ipp32f), { 10,10 }, 3, 3); 
	ippiCopyWrapBorder_32s_C1R((Ipp32s*)pDetailYDst, detailYStep, dstRoiSize, (Ipp32s*)pYB, 10 * sizeof(Ipp32f), { 10,10 }, 3, 3);
	ippiCopyWrapBorder_32s_C1R((Ipp32s*)pDetailXYDst, detailXYStep, dstRoiSize, (Ipp32s*)pXYB, 10 * sizeof(Ipp32f), { 10,10 }, 3, 3);

	printf("\n\npApproxDst: \n");
	for (int i = 0; i < 4 * 4; i++)
	{
		if (i % 4 == 0)
			printf("\n");
		printf("%-*.2f", 10, pApproxDst[i]);
	}
	printf("\n\nAfter Adding Border, pApproxDst: \n");
	for (int i = 0; i < 10 * 10; i++)
	{
	if (i % 10 == 0)
	printf("\n");
	printf("%-*.2f", 10, pAppB[i]);
	}

	printf("\n\npDetailXDst: \n");
	for (int i = 0; i < 4 * 4; i++)
	{
	if (i % 4 == 0)
	printf("\n");
	printf("%-*.2f", 10, pDetailXDst[i]);
	}
	printf("\n\nAfter Adding Border, pDetailXDst: \n");
	for (int i = 0; i < 10 * 10; i++)
	{
	if (i % 10 == 0)
	printf("\n");
	printf("%-*.2f", 10, pXB[i]);
	}

	printf("\n\npDetailYDst: \n");
	for (int i = 0; i < 4 * 4; i++)
	{
	if (i % 4 == 0)
	printf("\n");
	printf("%-*.2f", 10, pDetailYDst[i]);
	}
	printf("\n\nAfter Adding Border, pDetailYDst: \n");
	for (int i = 0; i < 10 * 10; i++)
	{
	if (i % 10 == 0)
	printf("\n");
	printf("%-*.2f", 10, pYB[i]);
	}


	// performs inverse wavelet transform   
	ippiWTInv_32f_C1R(
		pAppB + 33, 10 * sizeof(Ipp32f),
		pXB + 33, 10 * sizeof(Ipp32f),
		pYB + 33, 10 * sizeof(Ipp32f),
		pXYB + 33, 10 * sizeof(Ipp32f),
		{ 4,4 }, pDstInv, stepDstInv, pSpecInv, pBufferInv);

	printf("\n\nAfter WTFinv: \n");
	for (int i = 0; i < 8 * 8; i++)
	{
		if (i % 8 == 0)
			printf("\n");
		printf("%-*.2f", 10, pDstInv[i]);
	}

	ippFree(pSpecFwd);
	ippFree(pSpecInv);
	ippFree(pBufferFwd);
	ippFree(pBufferInv);
}

int main()
{
	func_wavelet();
	printf("\n \nend \n");
	return 0;
}
