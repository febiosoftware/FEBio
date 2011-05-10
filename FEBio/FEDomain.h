#pragma once

#include "FECore/FEElement.h"
#include "FECore/vector.h"
#include "FECore/DumpFile.h"
#include "FECore/FE_enum.h"

class FENode;
class FEMesh;
class FESolidSolver;
class FEHeatSolver;
class FEMaterial;

//-----------------------------------------------------------------------------
//! This class describes a physical domain that will be divided into elements
//! of a specific type. All elements in the domain have to have the same type
//!
class FEDomain
{
public:
	FEDomain(int ntype, FEMesh* pm, FEMaterial* pmat) { m_pMesh = pm; m_ntype = ntype; m_pMat = pmat; }
	virtual ~FEDomain() {}

	int Type() { return m_ntype; }

	void SetMesh(FEMesh* pm) { m_pMesh = pm; }
	FEMesh* GetMesh() { return m_pMesh; }

	void SetMaterial(FEMaterial* pmat) { m_pMat = pmat; }
	FEMaterial* GetMaterial() { return m_pMat; }

	virtual void create(int n) = 0;

	virtual int Elements() = 0;

	virtual FEDomain* Clone() { assert(false); return 0; }

	virtual void Reset() {}

	virtual void Serialize(DumpFile& ar) {}

	virtual bool Initialize(FEM& fem) { return true; }

	// TODO: this is temporary and will be moved to a different class
	virtual void UpdateStresses(FEM& fem) {}

	virtual void InitElements() {}

	virtual void StiffnessMatrix(FESolidSolver* psolver) {}

	virtual void Residual(FESolidSolver* psolver, vector<double>& R) {}

	//!< Initialize material point data for the elements
	void InitMaterialPointData();

	// TODO: this is not the preferred interface but I've added it for now
	virtual FEElement& ElementRef(int i) = 0;

	FEElement* FindElementFromID(int nid);

	virtual void UnpackElement(FEElement& el, unsigned int nflags = FE_UNPACK_ALL) = 0;

	void SetMatID(int mid)
	{
		for (int i=0; i<Elements(); ++i) ElementRef(i).SetMatID(mid);
	}

	virtual int Nodes() = 0;
	virtual FENode& Node(int i) = 0;

protected:
	FEMesh*		m_pMesh;	//!< the mesh that this domain is a part of
	FEMaterial*	m_pMat;		//!< the material for this domain

	int	m_ntype;
};
