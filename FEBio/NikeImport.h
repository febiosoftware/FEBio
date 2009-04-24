#pragma once

#include "FileImport.h"

//-----------------------------------------------------------------------------
//! Implements a class to import NIKE files

class FENIKEImport : public FEFileImport
{
protected:
	struct SLIDING_INTERFACE
	{
		int		nss;		// nr of slave elements
		int		nms;		// nr of master elements
		int		itype;		// type of sliding interface
		double	sfac;		// penalty scale factor
		double	mus;		// static friction coeff
		int		iaug;		// nr of augmentations
		double	altoln;		// normal augmentation tolerance
		double	altolt;		// tangent augmentation tolerance
		double	tkmult;		// tangent stiffness multiplier
	};

public:
	bool Load(FEM& fem, const char* szfile);

private:
	char* read_line(FILE* fp, char* szline, int n);

	bool ReadControlDeck (FEM& fem);
	bool ReadMaterialDeck(FEM& fem);
	bool ReadGeometry    (FEM& fem);
	bool ReadCurveDeck   (FEM& fem);
	bool ReadBCDecks     (FEM& fem);

protected:
	int	m_nrn;		//!< nr of rigid nodes cards
	int	m_nmat;		//!< nr of materials
	int	m_nlc;		//!< nr of loadcurves
	int m_nn;		//!< nr of nodes
	int m_nbel;		//!< nr of brick elements
	int	m_nsel;		//!< nr of shell elements
	int m_numsi;	//!< nr of sliding interfaces
	int m_ndis;		//!< nr of displacement cards
	int m_ncnf;		//!< nr of nodal force cards
	int m_npr;		//!< nr of face pressure cards
};
