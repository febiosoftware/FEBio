// FEElementTraits.cpp: implementation of the FEElementTraits class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEElementTraits.h"
#include "FEElement.h"
#include "FEException.h"

//*****************************************************************************
//                          F E H E X E L E M E N T
//*****************************************************************************

void FEHexElementTraits::init()
{
	int n;
	
	// integration point coordinates
	const double a = 1.0 / sqrt(3.0);
	gr[0] = -a; gs[0] = -a; gt[0] = -a; gw[0] = 1;
	gr[1] =  a; gs[1] = -a; gt[1] = -a; gw[1] = 1;
	gr[2] =  a; gs[2] =  a; gt[2] = -a; gw[2] = 1;
	gr[3] = -a; gs[3] =  a; gt[3] = -a; gw[3] = 1;
	gr[4] = -a; gs[4] = -a; gt[4] =  a; gw[4] = 1;
	gr[5] =  a; gs[5] = -a; gt[5] =  a; gw[5] = 1;
	gr[6] =  a; gs[6] =  a; gt[6] =  a; gw[6] = 1;
	gr[7] = -a; gs[7] =  a; gt[7] =  a; gw[7] = 1;
	
	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 - gt[n]);
		H[n][1] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 - gt[n]);
		H[n][2] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 - gt[n]);
		H[n][3] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 - gt[n]);
		H[n][4] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 + gt[n]);
		H[n][5] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 + gt[n]);
		H[n][6] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 + gt[n]);
		H[n][7] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 + gt[n]);
	}

	Hi = H.inverse();
	
	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.125*(1 - gs[n])*(1 - gt[n]);
		Gr[n][1] =  0.125*(1 - gs[n])*(1 - gt[n]);
		Gr[n][2] =  0.125*(1 + gs[n])*(1 - gt[n]);
		Gr[n][3] = -0.125*(1 + gs[n])*(1 - gt[n]);
		Gr[n][4] = -0.125*(1 - gs[n])*(1 + gt[n]);
		Gr[n][5] =  0.125*(1 - gs[n])*(1 + gt[n]);
		Gr[n][6] =  0.125*(1 + gs[n])*(1 + gt[n]);
		Gr[n][7] = -0.125*(1 + gs[n])*(1 + gt[n]);
		
		Gs[n][0] = -0.125*(1 - gr[n])*(1 - gt[n]);
		Gs[n][1] = -0.125*(1 + gr[n])*(1 - gt[n]);
		Gs[n][2] =  0.125*(1 + gr[n])*(1 - gt[n]);
		Gs[n][3] =  0.125*(1 - gr[n])*(1 - gt[n]);
		Gs[n][4] = -0.125*(1 - gr[n])*(1 + gt[n]);
		Gs[n][5] = -0.125*(1 + gr[n])*(1 + gt[n]);
		Gs[n][6] =  0.125*(1 + gr[n])*(1 + gt[n]);
		Gs[n][7] =  0.125*(1 - gr[n])*(1 + gt[n]);
		
		Gt[n][0] = -0.125*(1 - gr[n])*(1 - gs[n]);
		Gt[n][1] = -0.125*(1 + gr[n])*(1 - gs[n]);
		Gt[n][2] = -0.125*(1 + gr[n])*(1 + gs[n]);
		Gt[n][3] = -0.125*(1 - gr[n])*(1 + gs[n]);
		Gt[n][4] =  0.125*(1 - gr[n])*(1 - gs[n]);
		Gt[n][5] =  0.125*(1 + gr[n])*(1 - gs[n]);
		Gt[n][6] =  0.125*(1 + gr[n])*(1 + gs[n]);
		Gt[n][7] =  0.125*(1 - gr[n])*(1 + gs[n]);
	}
	
	// calculate local second derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Grr[n][0] = 0.0;
		Grr[n][1] = 0.0;
		Grr[n][2] = 0.0;
		Grr[n][3] = 0.0;
		Grr[n][4] = 0.0;
		Grr[n][5] = 0.0;
		Grr[n][6] = 0.0;
		Grr[n][7] = 0.0;
		
		Gsr[n][0] = 0.125*(1 - gt[n]);
		Gsr[n][1] = -0.125*(1 - gt[n]);
		Gsr[n][2] =  0.125*(1 - gt[n]);
		Gsr[n][3] = -0.125*(1 - gt[n]);
		Gsr[n][4] =  0.125*(1 + gt[n]);
		Gsr[n][5] = -0.125*(1 + gt[n]);
		Gsr[n][6] =  0.125*(1 + gt[n]);
		Gsr[n][7] = -0.125*(1 + gt[n]);
		
		Gtr[n][0] =  0.125*(1 - gs[n]);
		Gtr[n][1] = -0.125*(1 - gs[n]);
		Gtr[n][2] = -0.125*(1 + gs[n]);
		Gtr[n][3] =  0.125*(1 + gs[n]);
		Gtr[n][4] = -0.125*(1 - gs[n]);
		Gtr[n][5] =  0.125*(1 - gs[n]);
		Gtr[n][6] =  0.125*(1 + gs[n]);
		Gtr[n][7] = -0.125*(1 + gs[n]);
		
		Grs[n][0] =  0.125*(1 - gt[n]);
		Grs[n][1] = -0.125*(1 - gt[n]);
		Grs[n][2] =  0.125*(1 - gt[n]);
		Grs[n][3] = -0.125*(1 - gt[n]);
		Grs[n][4] =  0.125*(1 + gt[n]);
		Grs[n][5] = -0.125*(1 + gt[n]);
		Grs[n][6] =  0.125*(1 + gt[n]);
		Grs[n][7] = -0.125*(1 + gt[n]);
		
		Gss[n][0] = 0.0;
		Gss[n][1] = 0.0;
		Gss[n][2] = 0.0;
		Gss[n][3] = 0.0;
		Gss[n][4] = 0.0;
		Gss[n][5] = 0.0;
		Gss[n][6] = 0.0;
		Gss[n][7] = 0.0;
		
		Gts[n][0] =  0.125*(1 - gr[n]);
		Gts[n][1] =  0.125*(1 + gr[n]);
		Gts[n][2] = -0.125*(1 + gr[n]);
		Gts[n][3] = -0.125*(1 - gr[n]);
		Gts[n][4] = -0.125*(1 - gr[n]);
		Gts[n][5] = -0.125*(1 + gr[n]);
		Gts[n][6] =  0.125*(1 + gr[n]);
		Gts[n][7] =  0.125*(1 - gr[n]);
		
		Grt[n][0] =  0.125*(1 - gs[n]);
		Grt[n][1] = -0.125*(1 - gs[n]);
		Grt[n][2] = -0.125*(1 + gs[n]);
		Grt[n][3] =  0.125*(1 + gs[n]);
		Grt[n][4] = -0.125*(1 - gs[n]);
		Grt[n][5] =  0.125*(1 - gs[n]);
		Grt[n][6] =  0.125*(1 + gs[n]);
		Grt[n][7] = -0.125*(1 + gs[n]);
		
		Gst[n][0] =  0.125*(1 - gr[n]);
		Gst[n][1] =  0.125*(1 + gr[n]);
		Gst[n][2] = -0.125*(1 + gr[n]);
		Gst[n][3] = -0.125*(1 - gr[n]);
		Gst[n][4] = -0.125*(1 - gr[n]);
		Gst[n][5] = -0.125*(1 + gr[n]);
		Gst[n][6] =  0.125*(1 + gr[n]);
		Gst[n][7] =  0.125*(1 - gr[n]);
		
		Gtt[n][0] = 0.0;
		Gtt[n][1] = 0.0;
		Gtt[n][2] = 0.0;
		Gtt[n][3] = 0.0;
		Gtt[n][4] = 0.0;
		Gtt[n][5] = 0.0;
		Gtt[n][6] = 0.0;
		Gtt[n][7] = 0.0;
	}
	
}

//*****************************************************************************
//                          F E H E X 2 0 E L E M E N T
//*****************************************************************************

void FEHex20ElementTraits::init()
{
	int n;
	
	// integration point coordinates
	const double a = 0.774596669241483;
	const double w1 = 5.0 / 9.0;
	const double w2 = 8.9 / 9.0;
	gr[ 0] = -a; gs[ 0] = -a; gt[ 0] = -a; gw[ 0] = w1*w1*w1;
	gr[ 1] =  0; gs[ 1] = -a; gt[ 1] = -a; gw[ 1] = w2*w1*w1;
	gr[ 2] =  a; gs[ 2] = -a; gt[ 2] = -a; gw[ 2] = w1*w1*w1;
	gr[ 3] = -a; gs[ 3] =  0; gt[ 3] = -a; gw[ 3] = w1*w2*w1;
	gr[ 4] =  0; gs[ 4] =  0; gt[ 4] = -a; gw[ 4] = w2*w2*w1;
	gr[ 5] =  a; gs[ 5] =  0; gt[ 5] = -a; gw[ 5] = w1*w2*w1;
	gr[ 6] = -a; gs[ 6] =  a; gt[ 6] = -a; gw[ 6] = w1*w1*w1;
	gr[ 7] =  0; gs[ 7] =  a; gt[ 7] = -a; gw[ 7] = w2*w1*w1;
	gr[ 8] =  a; gs[ 8] =  a; gt[ 8] = -a; gw[ 8] = w1*w1*w1;
	gr[ 9] = -a; gs[ 9] = -a; gt[ 9] =  0; gw[ 9] = w1*w1*w2;
	gr[10] =  0; gs[10] = -a; gt[10] =  0; gw[10] = w2*w1*w2;
	gr[11] =  a; gs[11] = -a; gt[11] =  0; gw[11] = w1*w1*w2;
	gr[12] = -a; gs[12] =  0; gt[12] =  0; gw[12] = w1*w2*w2;
	gr[13] =  0; gs[13] =  0; gt[13] =  0; gw[13] = w2*w2*w2;
	gr[14] =  a; gs[14] =  0; gt[14] =  0; gw[14] = w1*w2*w2;
	gr[15] = -a; gs[15] =  a; gt[15] =  0; gw[15] = w1*w1*w2;
	gr[16] =  0; gs[16] =  a; gt[16] =  0; gw[16] = w2*w1*w2;
	gr[17] =  a; gs[17] =  a; gt[17] =  0; gw[17] = w1*w1*w2;
	gr[18] = -a; gs[18] = -a; gt[18] =  a; gw[18] = w1*w1*w1;
	gr[19] =  0; gs[19] = -a; gt[19] =  a; gw[19] = w2*w1*w1;
	gr[20] =  a; gs[20] = -a; gt[20] =  a; gw[20] = w1*w1*w1;
	gr[21] = -a; gs[21] =  0; gt[21] =  a; gw[21] = w1*w2*w1;
	gr[22] =  0; gs[22] =  0; gt[22] =  a; gw[22] = w2*w2*w1;
	gr[23] =  a; gs[23] =  0; gt[23] =  a; gw[23] = w1*w2*w1;
	gr[24] = -a; gs[24] =  a; gt[24] =  a; gw[24] = w1*w1*w1;
	gr[25] =  0; gs[25] =  a; gt[25] =  a; gw[25] = w2*w1*w1;
	gr[26] =  a; gs[26] =  a; gt[26] =  a; gw[26] = w1*w1*w1;
	
	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		H[n][ 8] = 0.25*(1 - gr[n]*gr[n])*(1 - gs[n])*(1 - gt[n]);
		H[n][ 9] = 0.25*(1 - gs[n]*gs[n])*(1 + gr[n])*(1 - gt[n]);
		H[n][10] = 0.25*(1 - gr[n]*gr[n])*(1 + gs[n])*(1 - gt[n]);
		H[n][11] = 0.25*(1 - gs[n]*gs[n])*(1 - gr[n])*(1 - gt[n]);
		H[n][12] = 0.25*(1 - gr[n]*gr[n])*(1 - gs[n])*(1 + gt[n]);
		H[n][13] = 0.25*(1 - gs[n]*gs[n])*(1 + gr[n])*(1 + gt[n]);
		H[n][14] = 0.25*(1 - gr[n]*gr[n])*(1 + gs[n])*(1 + gt[n]);
		H[n][15] = 0.25*(1 - gs[n]*gs[n])*(1 - gr[n])*(1 + gt[n]);
		H[n][16] = 0.25*(1 - gt[n]*gt[n])*(1 - gr[n])*(1 - gs[n]);
		H[n][17] = 0.25*(1 - gt[n]*gt[n])*(1 + gr[n])*(1 - gs[n]);
		H[n][18] = 0.25*(1 - gt[n]*gt[n])*(1 + gr[n])*(1 + gs[n]);
		H[n][19] = 0.25*(1 - gt[n]*gt[n])*(1 - gr[n])*(1 + gs[n]);

		H[n][0] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 - gt[n]) - 0.5*(H[n][ 8] + H[n][11] + H[n][16]);
		H[n][1] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 - gt[n]) - 0.5*(H[n][ 8] + H[n][ 9] + H[n][17]);
		H[n][2] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 - gt[n]) - 0.5*(H[n][ 9] + H[n][10] + H[n][18]);
		H[n][3] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 - gt[n]) - 0.5*(H[n][10] + H[n][11] + H[n][19]);
		H[n][4] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 + gt[n]) - 0.5*(H[n][12] + H[n][15] + H[n][16]);
		H[n][5] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 + gt[n]) - 0.5*(H[n][12] + H[n][13] + H[n][17]);
		H[n][6] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 + gt[n]) - 0.5*(H[n][13] + H[n][14] + H[n][18]);
		H[n][7] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 + gt[n]) - 0.5*(H[n][14] + H[n][15] + H[n][19]);
	}

//	Hi = H.inverse();
	
	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][ 8] = -0.5*gr[n]*(1 - gs[n])*(1 - gt[n]);
		Gr[n][ 9] =  0.25*(1 - gs[n]*gs[n])*(1 - gt[n]);
		Gr[n][10] = -0.5*gr[n]*(1 + gs[n])*(1 - gt[n]);
		Gr[n][11] = -0.25*(1 - gs[n]*gs[n])*(1 - gt[n]);
		Gr[n][12] = -0.5*gr[n]*(1 - gs[n])*(1 + gt[n]);
		Gr[n][13] = 0.25*(1 - gs[n]*gs[n])*(1 + gt[n]);
		Gr[n][14] = -0.5*gr[n]*(1 + gs[n])*(1 + gt[n]);
		Gr[n][15] = -0.25*(1 - gs[n]*gs[n])*(1 + gt[n]);
		Gr[n][16] = -0.25*(1 - gt[n]*gt[n])*(1 - gs[n]);
		Gr[n][17] =  0.25*(1 - gt[n]*gt[n])*(1 - gs[n]);
		Gr[n][18] =  0.25*(1 - gt[n]*gt[n])*(1 + gs[n]);
		Gr[n][19] = -0.25*(1 - gt[n]*gt[n])*(1 + gs[n]);

		Gr[n][0] = -0.125*(1 - gs[n])*(1 - gt[n]) - 0.5*(Gr[n][ 8] + Gr[n][11] + Gr[n][16]);
		Gr[n][1] =  0.125*(1 - gs[n])*(1 - gt[n]) - 0.5*(Gr[n][ 8] + Gr[n][ 9] + Gr[n][17]);
		Gr[n][2] =  0.125*(1 + gs[n])*(1 - gt[n]) - 0.5*(Gr[n][ 9] + Gr[n][10] + Gr[n][18]);
		Gr[n][3] = -0.125*(1 + gs[n])*(1 - gt[n]) - 0.5*(Gr[n][10] + Gr[n][11] + Gr[n][19]);
		Gr[n][4] = -0.125*(1 - gs[n])*(1 + gt[n]) - 0.5*(Gr[n][12] + Gr[n][15] + Gr[n][16]);
		Gr[n][5] =  0.125*(1 - gs[n])*(1 + gt[n]) - 0.5*(Gr[n][12] + Gr[n][13] + Gr[n][17]);
		Gr[n][6] =  0.125*(1 + gs[n])*(1 + gt[n]) - 0.5*(Gr[n][13] + Gr[n][14] + Gr[n][18]);
		Gr[n][7] = -0.125*(1 + gs[n])*(1 + gt[n]) - 0.5*(Gr[n][14] + Gr[n][15] + Gr[n][19]);
		
		Gs[n][ 8] = -0.25*(1 - gr[n]*gr[n])*(1 - gt[n]);
		Gs[n][ 9] = -0.5*gs[n]*(1 + gr[n])*(1 - gt[n]);
		Gs[n][10] = 0.25*(1 - gr[n]*gr[n])*(1 - gt[n]);
		Gs[n][11] = -0.5*gs[n]*(1 - gr[n])*(1 - gt[n]);
		Gs[n][12] = -0.25*(1 - gr[n]*gr[n])*(1 + gt[n]);
		Gs[n][13] = -0.5*gs[n]*(1 + gr[n])*(1 + gt[n]);
		Gs[n][14] = 0.25*(1 - gr[n]*gr[n])*(1 + gt[n]);
		Gs[n][15] = -0.5*gs[n]*(1 - gr[n])*(1 + gt[n]);
		Gs[n][16] = -0.25*(1 - gt[n]*gt[n])*(1 - gr[n]);
		Gs[n][17] = -0.25*(1 - gt[n]*gt[n])*(1 + gr[n]);
		Gs[n][18] =  0.25*(1 - gt[n]*gt[n])*(1 + gr[n]);
		Gs[n][19] =  0.25*(1 - gt[n]*gt[n])*(1 - gr[n]);

		Gs[n][0] = -0.125*(1 - gr[n])*(1 - gt[n]) - 0.5*(Gs[n][ 8] + Gs[n][11] + Gs[n][16]);
		Gs[n][1] = -0.125*(1 + gr[n])*(1 - gt[n]) - 0.5*(Gs[n][ 8] + Gs[n][ 9] + Gs[n][17]);
		Gs[n][2] =  0.125*(1 + gr[n])*(1 - gt[n]) - 0.5*(Gs[n][ 9] + Gs[n][10] + Gs[n][18]);
		Gs[n][3] =  0.125*(1 - gr[n])*(1 - gt[n]) - 0.5*(Gs[n][10] + Gs[n][11] + Gs[n][19]);
		Gs[n][4] = -0.125*(1 - gr[n])*(1 + gt[n]) - 0.5*(Gs[n][12] + Gs[n][15] + Gs[n][16]);
		Gs[n][5] = -0.125*(1 + gr[n])*(1 + gt[n]) - 0.5*(Gs[n][12] + Gs[n][13] + Gs[n][17]);
		Gs[n][6] =  0.125*(1 + gr[n])*(1 + gt[n]) - 0.5*(Gs[n][13] + Gs[n][14] + Gs[n][18]);
		Gs[n][7] =  0.125*(1 - gr[n])*(1 + gt[n]) - 0.5*(Gs[n][14] + Gs[n][15] + Gs[n][19]);

		Gt[n][ 8] = -0.25*(1 - gr[n]*gr[n])*(1 - gs[n]);
		Gt[n][ 9] = -0.25*(1 - gs[n]*gs[n])*(1 + gr[n]);
		Gt[n][10] = -0.25*(1 - gr[n]*gr[n])*(1 + gs[n]);
		Gt[n][11] = -0.25*(1 - gs[n]*gs[n])*(1 - gr[n]);
		Gt[n][12] =  0.25*(1 - gr[n]*gr[n])*(1 - gs[n]);
		Gt[n][13] =  0.25*(1 - gs[n]*gs[n])*(1 + gr[n]);
		Gt[n][14] =  0.25*(1 - gr[n]*gr[n])*(1 + gs[n]);
		Gt[n][15] =  0.25*(1 - gs[n]*gs[n])*(1 - gr[n]);
		Gt[n][16] = -0.5*gt[n]*(1 - gr[n])*(1 - gs[n]);
		Gt[n][17] = -0.5*gt[n]*(1 + gr[n])*(1 - gs[n]);
		Gt[n][18] = -0.5*gt[n]*(1 + gr[n])*(1 + gs[n]);
		Gt[n][19] = -0.5*gt[n]*(1 - gr[n])*(1 + gs[n]);
		
		Gt[n][0] = -0.125*(1 - gr[n])*(1 - gs[n]) - 0.5*(Gt[n][ 8] + Gt[n][11] + Gt[n][16]);
		Gt[n][1] = -0.125*(1 + gr[n])*(1 - gs[n]) - 0.5*(Gt[n][ 8] + Gt[n][ 9] + Gt[n][17]);
		Gt[n][2] = -0.125*(1 + gr[n])*(1 + gs[n]) - 0.5*(Gt[n][ 9] + Gt[n][10] + Gt[n][18]);
		Gt[n][3] = -0.125*(1 - gr[n])*(1 + gs[n]) - 0.5*(Gt[n][10] + Gt[n][11] + Gt[n][19]);
		Gt[n][4] =  0.125*(1 - gr[n])*(1 - gs[n]) - 0.5*(Gt[n][12] + Gt[n][15] + Gt[n][16]);
		Gt[n][5] =  0.125*(1 + gr[n])*(1 - gs[n]) - 0.5*(Gt[n][12] + Gt[n][13] + Gt[n][17]);
		Gt[n][6] =  0.125*(1 + gr[n])*(1 + gs[n]) - 0.5*(Gt[n][13] + Gt[n][14] + Gt[n][18]);
		Gt[n][7] =  0.125*(1 - gr[n])*(1 + gs[n]) - 0.5*(Gt[n][14] + Gt[n][15] + Gt[n][19]);
	}
	
	// TODO: calculate local second derivatives of shape functions at gauss points (need for biphasic problems)
}

//*****************************************************************************
//                          F E R I H E X E L E M E N T
//*****************************************************************************

void FERIHexElementTraits::init()
{
	int n;
	
	// This is for a six point integration rule
	// integration point coordinates
	const double a = 8.0 / 6.0;
	gr[0] = -1; gs[0] = 0; gt[0] = 0; gw[0] = a;
	gr[1] =  1; gs[1] = 0; gt[1] = 0; gw[1] = a;
	gr[2] =  0; gs[2] =-1; gt[2] = 0; gw[2] = a;
	gr[3] =  0; gs[3] = 1; gt[3] = 0; gw[3] = a;
	gr[4] =  0; gs[4] = 0; gt[4] =-1; gw[4] = a;
	gr[5] =  0; gs[5] = 0; gt[5] = 1; gw[5] = a;
	
	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 - gt[n]);
		H[n][1] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 - gt[n]);
		H[n][2] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 - gt[n]);
		H[n][3] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 - gt[n]);
		H[n][4] = 0.125*(1 - gr[n])*(1 - gs[n])*(1 + gt[n]);
		H[n][5] = 0.125*(1 + gr[n])*(1 - gs[n])*(1 + gt[n]);
		H[n][6] = 0.125*(1 + gr[n])*(1 + gs[n])*(1 + gt[n]);
		H[n][7] = 0.125*(1 - gr[n])*(1 + gs[n])*(1 + gt[n]);
	}
	
//	Hi = H.inverse();

	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.125*(1 - gs[n])*(1 - gt[n]);
		Gr[n][1] =  0.125*(1 - gs[n])*(1 - gt[n]);
		Gr[n][2] =  0.125*(1 + gs[n])*(1 - gt[n]);
		Gr[n][3] = -0.125*(1 + gs[n])*(1 - gt[n]);
		Gr[n][4] = -0.125*(1 - gs[n])*(1 + gt[n]);
		Gr[n][5] =  0.125*(1 - gs[n])*(1 + gt[n]);
		Gr[n][6] =  0.125*(1 + gs[n])*(1 + gt[n]);
		Gr[n][7] = -0.125*(1 + gs[n])*(1 + gt[n]);
		
		Gs[n][0] = -0.125*(1 - gr[n])*(1 - gt[n]);
		Gs[n][1] = -0.125*(1 + gr[n])*(1 - gt[n]);
		Gs[n][2] =  0.125*(1 + gr[n])*(1 - gt[n]);
		Gs[n][3] =  0.125*(1 - gr[n])*(1 - gt[n]);
		Gs[n][4] = -0.125*(1 - gr[n])*(1 + gt[n]);
		Gs[n][5] = -0.125*(1 + gr[n])*(1 + gt[n]);
		Gs[n][6] =  0.125*(1 + gr[n])*(1 + gt[n]);
		Gs[n][7] =  0.125*(1 - gr[n])*(1 + gt[n]);
		
		Gt[n][0] = -0.125*(1 - gr[n])*(1 - gs[n]);
		Gt[n][1] = -0.125*(1 + gr[n])*(1 - gs[n]);
		Gt[n][2] = -0.125*(1 + gr[n])*(1 + gs[n]);
		Gt[n][3] = -0.125*(1 - gr[n])*(1 + gs[n]);
		Gt[n][4] =  0.125*(1 - gr[n])*(1 - gs[n]);
		Gt[n][5] =  0.125*(1 + gr[n])*(1 - gs[n]);
		Gt[n][6] =  0.125*(1 + gr[n])*(1 + gs[n]);
		Gt[n][7] =  0.125*(1 - gr[n])*(1 + gs[n]);
	}
	
	// calculate local second derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Grr[n][0] = 0.0;
		Grr[n][1] = 0.0;
		Grr[n][2] = 0.0;
		Grr[n][3] = 0.0;
		Grr[n][4] = 0.0;
		Grr[n][5] = 0.0;
		Grr[n][6] = 0.0;
		Grr[n][7] = 0.0;
		
		Gsr[n][0] = 0.125*(1 - gt[n]);
		Gsr[n][1] = -0.125*(1 - gt[n]);
		Gsr[n][2] =  0.125*(1 - gt[n]);
		Gsr[n][3] = -0.125*(1 - gt[n]);
		Gsr[n][4] =  0.125*(1 + gt[n]);
		Gsr[n][5] = -0.125*(1 + gt[n]);
		Gsr[n][6] =  0.125*(1 + gt[n]);
		Gsr[n][7] = -0.125*(1 + gt[n]);
		
		Gtr[n][0] =  0.125*(1 - gs[n]);
		Gtr[n][1] = -0.125*(1 - gs[n]);
		Gtr[n][2] = -0.125*(1 + gs[n]);
		Gtr[n][3] =  0.125*(1 + gs[n]);
		Gtr[n][4] = -0.125*(1 - gs[n]);
		Gtr[n][5] =  0.125*(1 - gs[n]);
		Gtr[n][6] =  0.125*(1 + gs[n]);
		Gtr[n][7] = -0.125*(1 + gs[n]);
		
		Grs[n][0] =  0.125*(1 - gt[n]);
		Grs[n][1] = -0.125*(1 - gt[n]);
		Grs[n][2] =  0.125*(1 - gt[n]);
		Grs[n][3] = -0.125*(1 - gt[n]);
		Grs[n][4] =  0.125*(1 + gt[n]);
		Grs[n][5] = -0.125*(1 + gt[n]);
		Grs[n][6] =  0.125*(1 + gt[n]);
		Grs[n][7] = -0.125*(1 + gt[n]);
		
		Gss[n][0] = 0.0;
		Gss[n][1] = 0.0;
		Gss[n][2] = 0.0;
		Gss[n][3] = 0.0;
		Gss[n][4] = 0.0;
		Gss[n][5] = 0.0;
		Gss[n][6] = 0.0;
		Gss[n][7] = 0.0;
		
		Gts[n][0] =  0.125*(1 - gr[n]);
		Gts[n][1] =  0.125*(1 + gr[n]);
		Gts[n][2] = -0.125*(1 + gr[n]);
		Gts[n][3] = -0.125*(1 - gr[n]);
		Gts[n][4] = -0.125*(1 - gr[n]);
		Gts[n][5] = -0.125*(1 + gr[n]);
		Gts[n][6] =  0.125*(1 + gr[n]);
		Gts[n][7] =  0.125*(1 - gr[n]);
		
		Grt[n][0] =  0.125*(1 - gs[n]);
		Grt[n][1] = -0.125*(1 - gs[n]);
		Grt[n][2] = -0.125*(1 + gs[n]);
		Grt[n][3] =  0.125*(1 + gs[n]);
		Grt[n][4] = -0.125*(1 - gs[n]);
		Grt[n][5] =  0.125*(1 - gs[n]);
		Grt[n][6] =  0.125*(1 + gs[n]);
		Grt[n][7] = -0.125*(1 + gs[n]);
		
		Gst[n][0] =  0.125*(1 - gr[n]);
		Gst[n][1] =  0.125*(1 + gr[n]);
		Gst[n][2] = -0.125*(1 + gr[n]);
		Gst[n][3] = -0.125*(1 - gr[n]);
		Gst[n][4] = -0.125*(1 - gr[n]);
		Gst[n][5] = -0.125*(1 + gr[n]);
		Gst[n][6] =  0.125*(1 + gr[n]);
		Gst[n][7] =  0.125*(1 - gr[n]);
		
		Gtt[n][0] = 0.0;
		Gtt[n][1] = 0.0;
		Gtt[n][2] = 0.0;
		Gtt[n][3] = 0.0;
		Gtt[n][4] = 0.0;
		Gtt[n][5] = 0.0;
		Gtt[n][6] = 0.0;
		Gtt[n][7] = 0.0;
	}
	
}

//*****************************************************************************
//                          F E U D F H E X E L E M E N T
//*****************************************************************************

void FEUDFHexElementTraits::init()
{
	// single gauss-point integration rule
	gr[0] = 0; gs[0] = 0; gt[0] = 0; gw[0] = 8.0;
	
	// calculate shape function values at gauss points
	H[0][0] = 0.125*(1 - gr[0])*(1 - gs[0])*(1 - gt[0]);
	H[0][1] = 0.125*(1 + gr[0])*(1 - gs[0])*(1 - gt[0]);
	H[0][2] = 0.125*(1 + gr[0])*(1 + gs[0])*(1 - gt[0]);
	H[0][3] = 0.125*(1 - gr[0])*(1 + gs[0])*(1 - gt[0]);
	H[0][4] = 0.125*(1 - gr[0])*(1 - gs[0])*(1 + gt[0]);
	H[0][5] = 0.125*(1 + gr[0])*(1 - gs[0])*(1 + gt[0]);
	H[0][6] = 0.125*(1 + gr[0])*(1 + gs[0])*(1 + gt[0]);
	H[0][7] = 0.125*(1 - gr[0])*(1 + gs[0])*(1 + gt[0]);

//	Hi = H.inverse();
	
	// calculate local derivatives of shape functions at gauss points
	Gr[0][0] = -0.125*(1 - gs[0])*(1 - gt[0]);
	Gr[0][1] =  0.125*(1 - gs[0])*(1 - gt[0]);
	Gr[0][2] =  0.125*(1 + gs[0])*(1 - gt[0]);
	Gr[0][3] = -0.125*(1 + gs[0])*(1 - gt[0]);
	Gr[0][4] = -0.125*(1 - gs[0])*(1 + gt[0]);
	Gr[0][5] =  0.125*(1 - gs[0])*(1 + gt[0]);
	Gr[0][6] =  0.125*(1 + gs[0])*(1 + gt[0]);
	Gr[0][7] = -0.125*(1 + gs[0])*(1 + gt[0]);
	
	Gs[0][0] = -0.125*(1 - gr[0])*(1 - gt[0]);
	Gs[0][1] = -0.125*(1 + gr[0])*(1 - gt[0]);
	Gs[0][2] =  0.125*(1 + gr[0])*(1 - gt[0]);
	Gs[0][3] =  0.125*(1 - gr[0])*(1 - gt[0]);
	Gs[0][4] = -0.125*(1 - gr[0])*(1 + gt[0]);
	Gs[0][5] = -0.125*(1 + gr[0])*(1 + gt[0]);
	Gs[0][6] =  0.125*(1 + gr[0])*(1 + gt[0]);
	Gs[0][7] =  0.125*(1 - gr[0])*(1 + gt[0]);
	
	Gt[0][0] = -0.125*(1 - gr[0])*(1 - gs[0]);
	Gt[0][1] = -0.125*(1 + gr[0])*(1 - gs[0]);
	Gt[0][2] = -0.125*(1 + gr[0])*(1 + gs[0]);
	Gt[0][3] = -0.125*(1 - gr[0])*(1 + gs[0]);
	Gt[0][4] =  0.125*(1 - gr[0])*(1 - gs[0]);
	Gt[0][5] =  0.125*(1 + gr[0])*(1 - gs[0]);
	Gt[0][6] =  0.125*(1 + gr[0])*(1 + gs[0]);
	Gt[0][7] =  0.125*(1 - gr[0])*(1 + gs[0]);
	
	// calculate local second derivatives of shape functions at gauss points
	Grr[0][0] = 0.0;
	Grr[0][1] = 0.0;
	Grr[0][2] = 0.0;
	Grr[0][3] = 0.0;
	Grr[0][4] = 0.0;
	Grr[0][5] = 0.0;
	Grr[0][6] = 0.0;
	Grr[0][7] = 0.0;
	
	Gsr[0][0] = 0.125*(1 - gt[0]);
	Gsr[0][1] = -0.125*(1 - gt[0]);
	Gsr[0][2] =  0.125*(1 - gt[0]);
	Gsr[0][3] = -0.125*(1 - gt[0]);
	Gsr[0][4] =  0.125*(1 + gt[0]);
	Gsr[0][5] = -0.125*(1 + gt[0]);
	Gsr[0][6] =  0.125*(1 + gt[0]);
	Gsr[0][7] = -0.125*(1 + gt[0]);
	
	Gtr[0][0] =  0.125*(1 - gs[0]);
	Gtr[0][1] = -0.125*(1 - gs[0]);
	Gtr[0][2] = -0.125*(1 + gs[0]);
	Gtr[0][3] =  0.125*(1 + gs[0]);
	Gtr[0][4] = -0.125*(1 - gs[0]);
	Gtr[0][5] =  0.125*(1 - gs[0]);
	Gtr[0][6] =  0.125*(1 + gs[0]);
	Gtr[0][7] = -0.125*(1 + gs[0]);
	
	Grs[0][0] =  0.125*(1 - gt[0]);
	Grs[0][1] = -0.125*(1 - gt[0]);
	Grs[0][2] =  0.125*(1 - gt[0]);
	Grs[0][3] = -0.125*(1 - gt[0]);
	Grs[0][4] =  0.125*(1 + gt[0]);
	Grs[0][5] = -0.125*(1 + gt[0]);
	Grs[0][6] =  0.125*(1 + gt[0]);
	Grs[0][7] = -0.125*(1 + gt[0]);
	
	Gss[0][0] = 0.0;
	Gss[0][1] = 0.0;
	Gss[0][2] = 0.0;
	Gss[0][3] = 0.0;
	Gss[0][4] = 0.0;
	Gss[0][5] = 0.0;
	Gss[0][6] = 0.0;
	Gss[0][7] = 0.0;
	
	Gts[0][0] =  0.125*(1 - gr[0]);
	Gts[0][1] =  0.125*(1 + gr[0]);
	Gts[0][2] = -0.125*(1 + gr[0]);
	Gts[0][3] = -0.125*(1 - gr[0]);
	Gts[0][4] = -0.125*(1 - gr[0]);
	Gts[0][5] = -0.125*(1 + gr[0]);
	Gts[0][6] =  0.125*(1 + gr[0]);
	Gts[0][7] =  0.125*(1 - gr[0]);
	
	Grt[0][0] =  0.125*(1 - gs[0]);
	Grt[0][1] = -0.125*(1 - gs[0]);
	Grt[0][2] = -0.125*(1 + gs[0]);
	Grt[0][3] =  0.125*(1 + gs[0]);
	Grt[0][4] = -0.125*(1 - gs[0]);
	Grt[0][5] =  0.125*(1 - gs[0]);
	Grt[0][6] =  0.125*(1 + gs[0]);
	Grt[0][7] = -0.125*(1 + gs[0]);
	
	Gst[0][0] =  0.125*(1 - gr[0]);
	Gst[0][1] =  0.125*(1 + gr[0]);
	Gst[0][2] = -0.125*(1 + gr[0]);
	Gst[0][3] = -0.125*(1 - gr[0]);
	Gst[0][4] = -0.125*(1 - gr[0]);
	Gst[0][5] = -0.125*(1 + gr[0]);
	Gst[0][6] =  0.125*(1 + gr[0]);
	Gst[0][7] =  0.125*(1 - gr[0]);
	
	Gtt[0][0] = 0.0;
	Gtt[0][1] = 0.0;
	Gtt[0][2] = 0.0;
	Gtt[0][3] = 0.0;
	Gtt[0][4] = 0.0;
	Gtt[0][5] = 0.0;
	Gtt[0][6] = 0.0;
	Gtt[0][7] = 0.0;
}

//*****************************************************************************
//                          F E T E T E L E M E N T
//*****************************************************************************

void FETetElementTraits::init()
{
	int n;
	
	// gaussian integration for tetrahedral elements
	const double a = 0.58541020;
	const double b = 0.13819660;
	const double w = 1.0 / 24.0;
	
	gr[0] = b; gs[0] = b; gt[0] = b; gw[0] = w;
	gr[1] = a; gs[1] = b; gt[1] = b; gw[1] = w;
	gr[2] = b; gs[2] = a; gt[2] = b; gw[2] = w;
	gr[3] = b; gs[3] = b; gt[3] = a; gw[3] = w;
	
	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 1 - gr[n] - gs[n] - gt[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
		H[n][3] = gt[n];
	}

	Hi = H.inverse();
	
	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -1;
		Gr[n][1] =  1;
		Gr[n][2] =  0;
		Gr[n][3] =  0;
		
		Gs[n][0] = -1;
		Gs[n][1] =  0;
		Gs[n][2] =  1;
		Gs[n][3] =  0;
		
		Gt[n][0] = -1;
		Gt[n][1] =  0;
		Gt[n][2] =  0;
		Gt[n][3] =  1;
	}
	
	// calculate local second derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Grr[n][0] =  0.0;
		Grr[n][1] =  0.0;
		Grr[n][2] =  0.0;
		Grr[n][3] =  0.0;
		
		Gsr[n][0] =  0.0;
		Gsr[n][1] =  0.0;
		Gsr[n][2] =  0.0;
		Gsr[n][3] =  0.0;
		
		Gtr[n][0] =  0.0;
		Gtr[n][1] =  0.0;
		Gtr[n][2] =  0.0;
		Gtr[n][3] =  0.0;
		
		Grs[n][0] =  0.0;
		Grs[n][1] =  0.0;
		Grs[n][2] =  0.0;
		Grs[n][3] =  0.0;
		
		Gss[n][0] =  0.0;
		Gss[n][1] =  0.0;
		Gss[n][2] =  0.0;
		Gss[n][3] =  0.0;
		
		Gts[n][0] =  0.0;
		Gts[n][1] =  0.0;
		Gts[n][2] =  0.0;
		Gts[n][3] =  0.0;
		
		Grt[n][0] =  0.0;
		Grt[n][1] =  0.0;
		Grt[n][2] =  0.0;
		Grt[n][3] =  0.0;
		
		Gst[n][0] =  0.0;
		Gst[n][1] =  0.0;
		Gst[n][2] =  0.0;
		Gst[n][3] =  0.0;
		
		Gtt[n][0] =  0.0;
		Gtt[n][1] =  0.0;
		Gtt[n][2] =  0.0;
		Gtt[n][3] =  0.0;
	}
	
}


//*****************************************************************************
//                          F E G 1 T E T E L E M E N T
//*****************************************************************************

void FEG1TetElementTraits::init()
{
	// gaussian integration for tetrahedral elements
	const double a = 0.25;
	const double w = 1.0 / 6.0;
	
	gr[0] = a; gs[0] = a; gt[0] = a; gw[0] = w;
	
	// calculate shape function values at gauss points
	H[0][0] = 1 - gr[0] - gs[0] - gt[0];
	H[0][1] = gr[0];
	H[0][2] = gs[0];
	H[0][3] = gt[0];

//	Hi = H.inverse();
	
	// calculate local derivatives of shape functions at gauss points
	Gr[0][0] = -1;
	Gr[0][1] =  1;
	Gr[0][2] =  0;
	Gr[0][3] =  0;
	
	Gs[0][0] = -1;
	Gs[0][1] =  0;
	Gs[0][2] =  1;
	Gs[0][3] =  0;
	
	Gt[0][0] = -1;
	Gt[0][1] =  0;
	Gt[0][2] =  0;
	Gt[0][3] =  1;
	
	// calculate local second derivatives of shape functions at gauss points
	Grr[0][0] =  0.0;
	Grr[0][1] =  0.0;
	Grr[0][2] =  0.0;
	Grr[0][3] =  0.0;
	
	Gsr[0][0] =  0.0;
	Gsr[0][1] =  0.0;
	Gsr[0][2] =  0.0;
	Gsr[0][3] =  0.0;
	
	Gtr[0][0] =  0.0;
	Gtr[0][1] =  0.0;
	Gtr[0][2] =  0.0;
	Gtr[0][3] =  0.0;
	
	Grs[0][0] =  0.0;
	Grs[0][1] =  0.0;
	Grs[0][2] =  0.0;
	Grs[0][3] =  0.0;
	
	Gss[0][0] =  0.0;
	Gss[0][1] =  0.0;
	Gss[0][2] =  0.0;
	Gss[0][3] =  0.0;
	
	Gts[0][0] =  0.0;
	Gts[0][1] =  0.0;
	Gts[0][2] =  0.0;
	Gts[0][3] =  0.0;
	
	Grt[0][0] =  0.0;
	Grt[0][1] =  0.0;
	Grt[0][2] =  0.0;
	Grt[0][3] =  0.0;
	
	Gst[0][0] =  0.0;
	Gst[0][1] =  0.0;
	Gst[0][2] =  0.0;
	Gst[0][3] =  0.0;
	
	Gtt[0][0] =  0.0;
	Gtt[0][1] =  0.0;
	Gtt[0][2] =  0.0;
	Gtt[0][3] =  0.0;
	
}


//*****************************************************************************
//                          F E T E T 1 0 E L E M E N T
//*****************************************************************************
// I think this assumes that the tetrahedron in natural space is essentially
// a constant metric tet with the edge nodes at the center of the edges.
void FETet10ElementTraits::init()
{
	int n;
	
	// integration point coordinates
	const double a = 0.58541020;
	const double b = 0.13819660;
	const double w = 0.25 / 6.0;
	gr[ 0] = a; gs[ 0] = b; gt[ 0] = b; gw[ 0] = w;
	gr[ 1] = b; gs[ 1] = a; gt[ 1] = b; gw[ 1] = w;
	gr[ 2] = b; gs[ 2] = b; gt[ 2] = a; gw[ 2] = w;
	gr[ 3] = b; gs[ 3] = b; gt[ 3] = b; gw[ 3] = w;
	
	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		double r1 = 1.0 - gr[n] - gs[n] - gt[n];
		double r2 = gr[n];
		double r3 = gs[n];
		double r4 = gt[n];

		H[n][0] = r1*(2.0*r1 - 1.0);
		H[n][1] = r2*(2.0*r2 - 1.0);
		H[n][2] = r3*(2.0*r3 - 1.0);
		H[n][3] = r4*(2.0*r4 - 1.0);
		H[n][4] = 4.0*r1*r2;
		H[n][5] = 4.0*r2*r3;
		H[n][6] = 4.0*r3*r1;
		H[n][7] = 4.0*r1*r4;
		H[n][8] = 4.0*r2*r4;
		H[n][9] = 4.0*r3*r4;
	}

//	Hi = H.inverse();
	
	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -3.0 + 4.0*gr[n] + 4.0*(gs[n] + gt[n]);
		Gr[n][1] =  4.0*gr[n] - 1.0;
		Gr[n][2] =  0.0;
		Gr[n][3] =  0.0;
		Gr[n][4] =  4.0 - 8.0*gr[n] - 4.0*(gs[n] + gt[n]);
		Gr[n][5] =  4.0*gs[n];
		Gr[n][6] = -4.0*gs[n];
		Gr[n][7] = -4.0*gt[n];
		Gr[n][8] =  4.0*gt[n];
		Gr[n][9] =  0.0;

		Gs[n][0] = -3.0 + 4.0*gs[n] + 4.0*(gr[n] + gt[n]);
		Gs[n][1] =  0.0;
		Gs[n][2] =  4.0*gs[n] - 1.0;
		Gs[n][3] =  0.0;
		Gs[n][4] = -4.0*gr[n];
		Gs[n][5] =  4.0*gr[n];
		Gs[n][6] =  4.0 - 8.0*gs[n] - 4.0*(gr[n] + gt[n]);
		Gs[n][7] = -4.0*gt[n];
		Gs[n][8] =  0.0;
		Gs[n][9] =  4.0*gt[n];

		Gt[n][0] = -3.0 + 4.0*gt[n] + 4.0*(gr[n] + gs[n]);
		Gt[n][1] =  0.0;
		Gt[n][2] =  0.0;
		Gt[n][3] =  4.0*gt[n] - 1.0;
		Gt[n][4] = -4.0*gr[n];
		Gt[n][5] =  0.0;
		Gt[n][6] = -4.0*gs[n];
		Gt[n][7] =  4.0 - 8.0*gt[n] - 4.0*(gr[n] + gs[n]);
		Gt[n][8] =  4.0*gr[n];
		Gt[n][9] =  4.0*gs[n];
	}
	
	// TODO: calculate local second derivatives of shape functions at gauss points (need for biphasic problems)
}

//*****************************************************************************
//                          F E P E N T A E L E M E N T
//*****************************************************************************

void FEPentaElementTraits::init()
{
	int n;
	
	//gauss intergration points
	const double a = 1.0/6.0;
	const double b = 2.0/3.0;
	const double c = 1.0 / sqrt(3.0);
	const double w = 1.0 / 6.0;
	
	gr[0] = a; gs[0] = a; gt[0] = -c; gw[0] = w;
	gr[1] = b; gs[1] = a; gt[1] = -c; gw[1] = w;
	gr[2] = a; gs[2] = b; gt[2] = -c; gw[2] = w;
	gr[3] = a; gs[3] = a; gt[3] =  c; gw[3] = w;
	gr[4] = b; gs[4] = a; gt[4] =  c; gw[4] = w;
	gr[5] = a; gs[5] = b; gt[5] =  c; gw[5] = w;
	
	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.5*(1 - gt[n])*(1 - gr[n] - gs[n]);
		H[n][1] = 0.5*(1 - gt[n])*gr[n];
		H[n][2] = 0.5*(1 - gt[n])*gs[n];
		H[n][3] = 0.5*(1 + gt[n])*(1 - gr[n] - gs[n]);
		H[n][4] = 0.5*(1 + gt[n])*gr[n];
		H[n][5] = 0.5*(1 + gt[n])*gs[n];
	}

	Hi = H.inverse();
	
	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.5*(1 - gt[n]);
		Gr[n][1] =  0.5*(1 - gt[n]);
		Gr[n][2] =  0.0;
		Gr[n][3] = -0.5*(1 + gt[n]);
		Gr[n][4] =  0.5*(1 + gt[n]);
		Gr[n][5] =  0.0;
		
		Gs[n][0] = -0.5*(1 - gt[n]);
		Gs[n][1] =  0.0;
		Gs[n][2] =  0.5*(1 - gt[n]);
		Gs[n][3] = -0.5*(1 + gt[n]);
		Gs[n][4] =  0.0;
		Gs[n][5] =  0.5*(1 + gt[n]);
		
		Gt[n][0] = -0.5*(1 - gr[n] - gs[n]);
		Gt[n][1] = -0.5*gr[n];
		Gt[n][2] = -0.5*gs[n];
		Gt[n][3] =  0.5*(1 - gr[n] - gs[n]);
		Gt[n][4] =  0.5*gr[n];
		Gt[n][5] =  0.5*gs[n];
	}
	
	// calculate local second derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Grr[n][0] =  0.0;
		Grr[n][1] =  0.0;
		Grr[n][2] =  0.0;
		Grr[n][3] =  0.0;
		Grr[n][4] =  0.0;
		Grr[n][5] =  0.0;
		
		Gsr[n][0] =  0.0;
		Gsr[n][1] =  0.0;
		Gsr[n][2] =  0.0;
		Gsr[n][3] =  0.0;
		Gsr[n][4] =  0.0;
		Gsr[n][5] =  0.0;
		
		Gtr[n][0] =  0.5;
		Gtr[n][1] = -0.5;
		Gtr[n][2] =  0.0;
		Gtr[n][3] = -0.5;
		Gtr[n][4] =  0.5;
		Gtr[n][5] =  0.0;
		
		Grs[n][0] =  0.0;
		Grs[n][1] =  0.0;
		Grs[n][2] =  0.0;
		Grs[n][3] =  0.0;
		Grs[n][4] =  0.0;
		Grs[n][5] =  0.0;
		
		Gss[n][0] =  0.0;
		Gss[n][1] =  0.0;
		Gss[n][2] =  0.0;
		Gss[n][3] =  0.0;
		Gss[n][4] =  0.0;
		Gss[n][5] =  0.0;
		
		Gts[n][0] =  0.5;
		Gts[n][1] =  0.0;
		Gts[n][2] = -0.5;
		Gts[n][3] = -0.5;
		Gts[n][4] =  0.0;
		Gts[n][5] =  0.5;
		
		Grt[n][0] =  0.5;
		Grt[n][1] = -0.5;
		Grt[n][2] =  0.0;
		Grt[n][3] = -0.5;
		Grt[n][4] =  0.5;
		Grt[n][5] =  0.0;
		
		Gst[n][0] =  0.5;
		Gst[n][1] =  0.0;
		Gst[n][2] = -0.5;
		Gst[n][3] = -0.5;
		Gst[n][4] =  0.0;
		Gst[n][5] =  0.5;
		
		Gtt[n][0] =  0.0;
		Gtt[n][1] =  0.0;
		Gtt[n][2] =  0.0;
		Gtt[n][3] =  0.0;
		Gtt[n][4] =  0.0;
		Gtt[n][5] =  0.0;
	}
	
}

//*****************************************************************************
//                          F E Q U A D E L E M E N T
//*****************************************************************************

void FEQuadElementTraits::init()
{
	int n;
	
	const double a = 1.0 / sqrt(3.0);
	
	gr[0] = -a; gs[0] = -a; gw[0] = 1;
	gr[1] =  a; gs[1] = -a; gw[1] = 1;
	gr[2] =  a; gs[2] =  a; gw[2] = 1;
	gr[3] = -a; gs[3] =  a; gw[3] = 1;
	
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.25*(1-gr[n])*(1-gs[n]);
		H[n][1] = 0.25*(1+gr[n])*(1-gs[n]);
		H[n][2] = 0.25*(1+gr[n])*(1+gs[n]);
		H[n][3] = 0.25*(1-gr[n])*(1+gs[n]);
	}
	
	Hi = H.inverse();
	
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.25*(1-gs[n]);
		Gr[n][1] =  0.25*(1-gs[n]);
		Gr[n][2] =  0.25*(1+gs[n]);
		Gr[n][3] = -0.25*(1+gs[n]);
		
		Gs[n][0] = -0.25*(1-gr[n]);
		Gs[n][1] = -0.25*(1+gr[n]);
		Gs[n][2] =  0.25*(1+gr[n]);
		Gs[n][3] =  0.25*(1-gr[n]);
	}
}

void FEQuadElementTraits::shape(double* H, double r, double s)
{
	H[0] = 0.25*(1-r)*(1-s);
	H[1] = 0.25*(1+r)*(1-s);
	H[2] = 0.25*(1+r)*(1+s);
	H[3] = 0.25*(1-r)*(1+s);
}

void FEQuadElementTraits::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -0.25*(1-s); Hs[0] = -0.25*(1-r);
	Hr[1] =  0.25*(1-s); Hs[1] = -0.25*(1+r);
	Hr[2] =  0.25*(1+s); Hs[2] =  0.25*(1+r);
	Hr[3] = -0.25*(1+s); Hs[3] =  0.25*(1-r);
}

void FEQuadElementTraits::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] = 0; Hrs[0] =  0.25; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = -0.25; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] =  0.25; Hss[2] = 0;
	Hrr[3] = 0; Hrs[3] = -0.25; Hss[3] = 0;
}

void FEQuadElementTraits::project_to_nodes(double* ai, double* ao)
{
	int ni = NINT;
	int ne = NELN;
	assert(ni == ne);
	for (int i=0; i<ne; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<ni; ++j) ao[i] += Hi[i][j]*ai[j];
	}
}

//*****************************************************************************
//                          F E N I Q U A D E L E M E N T
//*****************************************************************************

void FENIQuadElementTraits::init()
{
	int n;
	
	gr[0] = -1; gs[0] = -1; gw[0] = 1;
	gr[1] =  1; gs[1] = -1; gw[1] = 1;
	gr[2] =  1; gs[2] =  1; gw[2] = 1;
	gr[3] = -1; gs[3] =  1; gw[3] = 1;
	
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.25*(1-gr[n])*(1-gs[n]);
		H[n][1] = 0.25*(1+gr[n])*(1-gs[n]);
		H[n][2] = 0.25*(1+gr[n])*(1+gs[n]);
		H[n][3] = 0.25*(1-gr[n])*(1+gs[n]);
	}
	
	Hi = H.inverse();
	
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -0.25*(1-gs[n]);
		Gr[n][1] =  0.25*(1-gs[n]);
		Gr[n][2] =  0.25*(1+gs[n]);
		Gr[n][3] = -0.25*(1+gs[n]);
		
		Gs[n][0] = -0.25*(1-gr[n]);
		Gs[n][1] = -0.25*(1+gr[n]);
		Gs[n][2] =  0.25*(1+gr[n]);
		Gs[n][3] =  0.25*(1-gr[n]);
	}
}

void FENIQuadElementTraits::shape(double* H, double r, double s)
{
	H[0] = 0.25*(1-r)*(1-s);
	H[1] = 0.25*(1+r)*(1-s);
	H[2] = 0.25*(1+r)*(1+s);
	H[3] = 0.25*(1-r)*(1+s);
}

void FENIQuadElementTraits::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -0.25*(1-s); Hs[0] = -0.25*(1-r);
	Hr[1] =  0.25*(1-s); Hs[1] = -0.25*(1+r);
	Hr[2] =  0.25*(1+s); Hs[2] =  0.25*(1+r);
	Hr[3] = -0.25*(1+s); Hs[3] =  0.25*(1-r);
}

void FENIQuadElementTraits::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] = 0; Hrs[0] =  0.25; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = -0.25; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] =  0.25; Hss[2] = 0;
	Hrr[3] = 0; Hrs[3] = -0.25; Hss[3] = 0;
}

void FENIQuadElementTraits::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[1];
	ao[2] = ai[2];
	ao[3] = ai[3];
}

//*****************************************************************************
//                          F E T R I E L E M E N T
//*****************************************************************************

void FETriElementTraits::init()
{
	int n;
	
	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;
	
	gr[0] = a; gs[0] = a; gw[0] = a;
	gr[1] = b; gs[1] = a; gw[1] = a;
	gr[2] = a; gs[2] = b; gw[2] = a;
	
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 1-gr[n]-gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}
	
	Hi = H.inverse();
	
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -1;
		Gr[n][1] =  1;
		Gr[n][2] =  0;
		
		Gs[n][0] = -1;
		Gs[n][1] =  0;
		Gs[n][2] =  1;
	}
}

void FETriElementTraits::shape(double* H, double r, double s)
{
	H[0] = 1.0 - r - s;
	H[1] = r;
	H[2] = s;
}

void FETriElementTraits::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -1; Hs[0] = -1;
	Hr[1] =  1; Hs[1] =  0;
	Hr[2] =  0; Hs[2] =  1;
}

void FETriElementTraits::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] = 0; Hrs[0] = 0; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = 0; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] = 0; Hss[2] = 0;
}

void FETriElementTraits::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[0];
	ao[2] = ai[0];
}

//*****************************************************************************
//                          F E N I T R I E L E M E N T
//*****************************************************************************

void FENITriElementTraits::init()
{
	int n;
	
	const double a = 1.0 / 6.0;
	
	gr[0] = 0; gs[0] = 0; gw[0] = a;
	gr[1] = 1; gs[1] = 0; gw[1] = a;
	gr[2] = 0; gs[2] = 1; gw[2] = a;
	
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 1-gr[n]-gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}
	
	Hi = H.inverse();
	
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -1;
		Gr[n][1] =  1;
		Gr[n][2] =  0;
		
		Gs[n][0] = -1;
		Gs[n][1] =  0;
		Gs[n][2] =  1;
	}
}

void FENITriElementTraits::shape(double* H, double r, double s)
{
	H[0] = 1.0 - r - s;
	H[1] = r;
	H[2] = s;
}

void FENITriElementTraits::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -1; Hs[0] = -1;
	Hr[1] =  1; Hs[1] =  0;
	Hr[2] =  0; Hs[2] =  1;
}

void FENITriElementTraits::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] = 0; Hrs[0] = 0; Hss[0] = 0;
	Hrr[1] = 0; Hrs[1] = 0; Hss[1] = 0;
	Hrr[2] = 0; Hrs[2] = 0; Hss[2] = 0;
}

void FENITriElementTraits::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[1];
	ao[2] = ai[2];
}

//*****************************************************************************
//                          F E T R I 6 E L E M E N T
//*****************************************************************************

void FETri6ElementTraits::init()
{
	int n;
	
	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;
	
	gr[0] = a; gs[0] = a; gw[0] = a;
	gr[1] = b; gs[1] = a; gw[1] = a;
	gr[2] = a; gs[2] = b; gw[2] = a;
	
	// calculate shape function values at gauss points
	for (n=0; n<NINT; ++n)
	{
		double r1 = 1.0 - gr[n] - gs[n];
		double r2 = gr[n];
		double r3 = gs[n];

		H[n][0] = r1*(2.0*r1 - 1.0);
		H[n][1] = r2*(2.0*r2 - 1.0);
		H[n][2] = r3*(2.0*r3 - 1.0);
		H[n][3] = 4.0*r1*r2;
		H[n][4] = 4.0*r2*r3;
		H[n][5] = 4.0*r3*r1;
	}

//	Hi = H.inverse();

	// calculate local derivatives of shape functions at gauss points
	for (n=0; n<NINT; ++n)
	{
		Gr[n][0] = -3.0 + 4.0*gr[n] + 4.0*gs[n];
		Gr[n][1] =  4.0*gr[n] - 1.0;
		Gr[n][2] =  0.0;
		Gr[n][3] =  4.0 - 8.0*gr[n] - 4.0*gs[n];
		Gr[n][4] =  4.0*gs[n];
		Gr[n][5] = -4.0*gs[n];

		Gs[n][0] = -3.0 + 4.0*gs[n] + 4.0*gr[n];
		Gs[n][1] =  0.0;
		Gs[n][2] =  4.0*gs[n] - 1.0;
		Gs[n][3] = -4.0*gr[n];
		Gs[n][4] =  4.0*gr[n];
		Gs[n][5] =  4.0 - 8.0*gs[n] - 4.0*gr[n];
	}
}

void FETri6ElementTraits::shape(double* H, double r, double s)
{
	double r1 = 1.0 - r - s;
	double r2 = r;
	double r3 = s;

	H[0] = r1*(2.0*r1 - 1.0);
	H[1] = r2*(2.0*r2 - 1.0);
	H[2] = r3*(2.0*r3 - 1.0);
	H[3] = 4.0*r1*r2;
	H[4] = 4.0*r2*r3;
	H[5] = 4.0*r3*r1;
}

void FETri6ElementTraits::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -3.0 + 4.0*r + 4.0*s;
	Hr[1] =  4.0*r - 1.0;
	Hr[2] =  0.0;
	Hr[3] =  4.0 - 8.0*r - 4.0*s;
	Hr[4] =  4.0*s;
	Hr[5] = -4.0*s;

	Hs[0] = -3.0 + 4.0*s + 4.0*r;
	Hs[1] =  0.0;
	Hs[2] =  4.0*s - 1.0;
	Hs[3] = -4.0*r;
	Hs[4] =  4.0*r;
	Hs[5] =  4.0 - 8.0*s - 4.0*r;
}

void FETri6ElementTraits::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] =  4.0; Hrs[0] =  4.0; Hss[0] =  4.0;
	Hrr[1] =  4.0; Hrs[1] =  0.0; Hss[1] =  0.0;
	Hrr[2] =  0.0; Hrs[2] =  0.0; Hss[2] =  4.0;
	Hrr[3] = -8.0; Hrs[3] = -4.0; Hss[3] =  0.0;
	Hrr[4] =  0.0; Hrs[4] =  4.0; Hss[4] =  0.0;
	Hrr[5] =  0.0; Hrs[5] = -4.0; Hss[5] = -8.0;
}

void FETri6ElementTraits::project_to_nodes(double* ai, double* ao)
{
	matrix H(3, 3);
	for (int n=0; n<3; ++n)
	{
		H[n][0] = 1.0 - gr[n] - gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}
	H.inverse();

	for (int i=0; i<3; ++i)
	{
		ao[i] = 0;
		for (int j=0; j<3; ++j) ao[i] += H[i][j]*ai[j];
	}

	ao[3] = 0.5*(ao[0] + ao[1]);
	ao[4] = 0.5*(ao[1] + ao[2]);
	ao[5] = 0.5*(ao[2] + ao[0]);
}

//*****************************************************************************
//                          F E N I T R I 6 E L E M E N T
//*****************************************************************************

void FENITri6ElementTraits::init()
{
	const double a = 0.0;
	const double b = 1.0/6.0;
	gr[0] = 0.0; gs[0] = 0.0; gw[0] = a;
	gr[1] = 1.0; gs[1] = 0.0; gw[1] = a;
	gr[2] = 0.0; gs[2] = 1.0; gw[2] = a;
	gr[3] = 0.5; gs[3] = 0.0; gw[3] = b;
	gr[4] = 0.5; gs[4] = 0.5; gw[4] = b;
	gr[5] = 0.0; gs[5] = 0.5; gw[5] = b;
	
	// calculate shape function values at gauss points
	for (int n=0; n<NINT; ++n)
	{
		double r1 = 1.0 - gr[n] - gs[n];
		double r2 = gr[n];
		double r3 = gs[n];

		H[n][0] = r1*(2.0*r1 - 1.0);
		H[n][1] = r2*(2.0*r2 - 1.0);
		H[n][2] = r3*(2.0*r3 - 1.0);
		H[n][3] = 4.0*r1*r2;
		H[n][4] = 4.0*r2*r3;
		H[n][5] = 4.0*r3*r1;
	}

//	Hi = H.inverse();

	// calculate local derivatives of shape functions at gauss points
	for (int n=0; n<NINT; ++n)
	{
		Gr[n][0] = -3.0 + 4.0*gr[n] + 4.0*gs[n];
		Gr[n][1] =  4.0*gr[n] - 1.0;
		Gr[n][2] =  0.0;
		Gr[n][3] =  4.0 - 8.0*gr[n] - 4.0*gs[n];
		Gr[n][4] =  4.0*gs[n];
		Gr[n][5] = -4.0*gs[n];

		Gs[n][0] = -3.0 + 4.0*gs[n] + 4.0*gr[n];
		Gs[n][1] =  0.0;
		Gs[n][2] =  4.0*gs[n] - 1.0;
		Gs[n][3] = -4.0*gr[n];
		Gs[n][4] =  4.0*gr[n];
		Gs[n][5] =  4.0 - 8.0*gs[n] - 4.0*gr[n];
	}
}

void FENITri6ElementTraits::shape(double* H, double r, double s)
{
	double r1 = 1.0 - r - s;
	double r2 = r;
	double r3 = s;

	H[0] = r1*(2.0*r1 - 1.0);
	H[1] = r2*(2.0*r2 - 1.0);
	H[2] = r3*(2.0*r3 - 1.0);
	H[3] = 4.0*r1*r2;
	H[4] = 4.0*r2*r3;
	H[5] = 4.0*r3*r1;
}

void FENITri6ElementTraits::shape_deriv(double* Hr, double* Hs, double r, double s)
{
	Hr[0] = -3.0 + 4.0*r + 4.0*s;
	Hr[1] =  4.0*r - 1.0;
	Hr[2] =  0.0;
	Hr[3] =  4.0 - 8.0*r - 4.0*s;
	Hr[4] =  4.0*s;
	Hr[5] = -4.0*s;

	Hs[0] = -3.0 + 4.0*s + 4.0*r;
	Hs[1] =  0.0;
	Hs[2] =  4.0*s - 1.0;
	Hs[3] = -4.0*r;
	Hs[4] =  4.0*r;
	Hs[5] =  4.0 - 8.0*s - 4.0*r;
}

void FENITri6ElementTraits::shape_deriv2(double* Hrr, double* Hrs, double* Hss, double r, double s)
{
	Hrr[0] =  4.0; Hrs[0] =  4.0; Hss[0] =  4.0;
	Hrr[1] =  4.0; Hrs[1] =  0.0; Hss[1] =  0.0;
	Hrr[2] =  0.0; Hrs[2] =  0.0; Hss[2] =  4.0;
	Hrr[3] = -8.0; Hrs[3] = -4.0; Hss[3] =  0.0;
	Hrr[4] =  0.0; Hrs[4] =  4.0; Hss[4] =  0.0;
	Hrr[5] =  0.0; Hrs[5] = -4.0; Hss[5] = -8.0;
}

void FENITri6ElementTraits::project_to_nodes(double* ai, double* ao)
{
	ao[0] = ai[0];
	ao[1] = ai[1];
	ao[2] = ai[2];
	ao[3] = ai[3];
	ao[4] = ai[4];
	ao[5] = ai[5];
}

//*****************************************************************************
//                          F E S H E L L Q U A D E L E M E N T
//*****************************************************************************

void FEShellQuadElementTraits::init()
{
	int n;
	
	const double a = 1.0 / sqrt(3.0);
	const double b = sqrt(3.0/5.0);
	const double w = 5.0 / 9.0;
	
	gr[ 0] = -a; gs[ 0] = -a; gt[ 0] = -b; gw[ 0] = w;
	gr[ 1] =  a; gs[ 1] = -a; gt[ 1] = -b; gw[ 1] = w;
	gr[ 2] =  a; gs[ 2] =  a; gt[ 2] = -b; gw[ 2] = w;
	gr[ 3] = -a; gs[ 3] =  a; gt[ 3] = -b; gw[ 3] = w;
	
	gr[ 4] = -a; gs[ 4] = -a; gt[ 4] =  0; gw[ 4] = 8.0/9.0;
	gr[ 5] =  a; gs[ 5] = -a; gt[ 5] =  0; gw[ 5] = 8.0/9.0;
	gr[ 6] =  a; gs[ 6] =  a; gt[ 6] =  0; gw[ 6] = 8.0/9.0;
	gr[ 7] = -a; gs[ 7] =  a; gt[ 7] =  0; gw[ 7] = 8.0/9.0;
	
	gr[ 8] = -a; gs[ 8] = -a; gt[ 8] =  b; gw[ 8] = w;
	gr[ 9] =  a; gs[ 9] = -a; gt[ 9] =  b; gw[ 9] = w;
	gr[10] =  a; gs[10] =  a; gt[10] =  b; gw[10] = w;
	gr[11] = -a; gs[11] =  a; gt[11] =  b; gw[11] = w;
	
	
	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 0.25*(1-gr[n])*(1-gs[n]);
		H[n][1] = 0.25*(1+gr[n])*(1-gs[n]);
		H[n][2] = 0.25*(1+gr[n])*(1+gs[n]);
		H[n][3] = 0.25*(1-gr[n])*(1+gs[n]);
	}

//	Hi = H.inverse();

	for (n=0; n<NINT; ++n)
	{
		Hr[n][0] = -0.25*(1-gs[n]);
		Hr[n][1] =  0.25*(1-gs[n]);
		Hr[n][2] =  0.25*(1+gs[n]);
		Hr[n][3] = -0.25*(1+gs[n]);
		
		Hs[n][0] = -0.25*(1-gr[n]);
		Hs[n][1] = -0.25*(1+gr[n]);
		Hs[n][2] =  0.25*(1+gr[n]);
		Hs[n][3] =  0.25*(1-gr[n]);
	}
}

//*****************************************************************************
//                          F E S H E L L T R I E L E M E N T
//*****************************************************************************

void FEShellTriElementTraits::init()
{
	int n;

	const double a = 1.0 / 6.0;
	const double b = 2.0 / 3.0;
	const double c = sqrt(3.0/5.0);
	const double w = 5.0 / 9.0;

	gr[0] = a; gs[0] = a; gt[0] = -b; gw[0] = a*w;
	gr[1] = b; gs[1] = a; gt[1] = -b; gw[1] = a*w;
	gr[2] = a; gs[2] = b; gt[2] = -b; gw[2] = a*w;

	gr[3] = a; gs[3] = a; gt[3] =  0; gw[3] = a*w;
	gr[4] = b; gs[4] = a; gt[4] =  0; gw[4] = a*w;
	gr[5] = a; gs[5] = b; gt[5] =  0; gw[5] = a*w;

	gr[6] = a; gs[6] = a; gt[6] =  b; gw[6] = a*w;
	gr[7] = b; gs[7] = a; gt[7] =  b; gw[7] = a*w;
	gr[8] = a; gs[8] = b; gt[8] =  b; gw[8] = a*w;

	for (n=0; n<NINT; ++n)
	{
		H[n][0] = 1-gr[n]-gs[n];
		H[n][1] = gr[n];
		H[n][2] = gs[n];
	}

//	Hi = H.inverse();

	for (n=0; n<NINT; ++n)
	{
		Hr[n][0] = -1;
		Hr[n][1] =  1;
		Hr[n][2] =  0;

		Hs[n][0] = -1;
		Hs[n][1] =  0;
		Hs[n][2] =  1;
	}
}

//*****************************************************************************
//                          F E T R U S S E L E M E N T
//*****************************************************************************

void FETrussElementTraits::init()
{

}
