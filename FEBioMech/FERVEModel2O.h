#pragma once
#include <FECore/FEModel.h>
#include <FECore/mat3d.h>
#include <FECore/tens3d.h>
#include <FECore/tens4d.h>
#include <FECore/tens5d.h>
#include <FECore/tens6d.h>

//-----------------------------------------------------------------------------
// Class describing the second-order RVE model.
// This is used by the second-order homogenization code.
class FERVEModel2O : public FEModel
{
public:
	enum RVE_TYPE
	{
		DISPLACEMENT,		// prescribed displacement
		PERIODIC_AL,		// periodic, augmented Lagrangian
		PERIODIC_LC			// periodic, linear constraints
	};

public:
	FERVEModel2O();
	~FERVEModel2O();

	//! one time initialization
	bool InitRVE(int rveType, const char* szbc);

	// scale the geometry
	void ScaleGeometry(double scale);

public:
	//! Return the initial volume (calculated in Init)
	double InitialVolume() const { return m_V0; }

	//! periodicity flag
	int RVEType() const { return m_rveType; }

	//! get boundary nodes
	const vector<int>& BoundaryList() const { return m_BN; }

protected:
	//! Calculate the initial volume
	void EvalInitialVolume();

	//! find the list of boundary nodes
	void FindBoundaryNodes(vector<int>& BN);

	//! Center the RVE
	void CenterRVE();

	bool PrepDisplacementBC();
	bool PrepPeriodicBC(const char* szbc);
	bool PrepPeriodicLC(const char* szbc);

private:
	double			m_V0;			//!< initial volume
	int				m_rveType;		//!< type of RVE
	FEBoundingBox	m_bb;			//!< bounding box of mesh
	vector<int>		m_BN;			//!< boundary node flags
};

//-----------------------------------------------------------------------------
//! class representing the RVE models at the material points.
class FEMicroModel2O : public FEModel
{
public:
	FEMicroModel2O();
	~FEMicroModel2O();

public:
	//! Initialize the micro model from the "master" RVE model
	bool Init(FERVEModel2O& rve);

	//! Solve the RVE model
	bool Solve(const mat3d& F, const tens3drs& G);

	//! calculate average PK1 stress
	mat3d AveragedStressPK1(FEMaterialPoint &mp);

	void AveragedStress2O   (mat3d&  Pa, tens3drs&  Qa);
//	void AveragedStress2OPK1(FEMaterialPoint &mp, mat3d& PK1a, tens3drs& QK1a);
//	void AveragedStress2OPK2(FEMaterialPoint &mp, mat3ds&  Sa, tens3ds&    Ta);
	
	void AveragedStiffness(FEMaterialPoint &mp, tens4d& C, tens5d& L, tens5d& H, tens6d& J);

protected:
	//! Update the boundary conditions
	void UpdateBC(const mat3d& F, const tens3drs& G);

	//! see if node is boundary node
	bool IsBoundaryNode(int i) const { return ((*m_BN)[i]==1); }

protected:
	int		m_rveType;		//!< RVE type flag
	double	m_V0;			//!< initial RVE volume
	const vector<int>*	m_BN;
};
