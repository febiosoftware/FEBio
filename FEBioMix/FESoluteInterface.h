#pragma once
#include <FECore/fecore_api.h>

class FESolute;

//------------------------------------------------------------------------
// This class should be used by all materials that support solutes
// TODO: This is a work in progress. The goal is to reduce the dynamic_casts to materials
//       that support solutes, and instead provide a single consistent interface to features
//       that need access to solute data (e.g. plot variables).
class FECORE_API FESoluteInterface
{
public:
	FESoluteInterface(){}
	virtual ~FESoluteInterface(){}

// derived classes need to implement the following functions
public:
	// return the number of solutes in the material
	virtual int Solutes() = 0;

	// return a solute material
	virtual FESolute* GetSolute(int i) = 0;

// additional member functions
public:
	// return the local index of a global solute ID (or -1, if the solute is not in this material)
	int FindLocalSoluteID(int soluteID);
};
