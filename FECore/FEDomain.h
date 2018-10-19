#pragma once
#include "FEMeshPartition.h"

// forward declaration of material class
class FEMaterial;

// Base class for solid and shell parts. Domains can also have materials assigned.
class FEDomain : public FEMeshPartition
{
public:
	FEDomain(int nclass, FEModel* fem);

	//! get the material of this domain
	virtual FEMaterial* GetMaterial() { return 0; }

	// assign a material to this domain
	virtual void SetMaterial(FEMaterial* pm);

	//! set the material ID of all elements
	void SetMatID(int mid);

	//! Allocate material point data for the elements
	//! This is called after elements get read in from the input file.
	//! And must be called before material point data can be accessed.
	//! \todo Perhaps I can make this part of the "creation" routine
	void CreateMaterialPointData();

	//! Initialize material points in the domain (optional)
	virtual void InitMaterialPoints() {}

	// serialization
	void Serialize(DumpStream& ar) override;
};
