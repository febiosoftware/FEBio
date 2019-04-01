#pragma once
#include <FECore/FEPrescribedBC.h>
#include <FECore/tens3d.h>

//-----------------------------------------------------------------------------
class FEBCPrescribedDeformation : public FEPrescribedBC
{
public:
	FEBCPrescribedDeformation(FEModel* pfem);

public:
	void AddNode(int n);
	void AddNodes(const FENodeSet& set) override;

	void Activate() override;

	void Deactivate() override;

	void PrepStep(std::vector<double>& ui, bool brel) override;

	void Update() override;

	int Items() const { return (int) m_node.size(); }

	int NodeID(int index) const { return m_node[index]; }

	void SetDeformationGradient(const mat3d& F);

	void CopyFrom(FEPrescribedBC* pbc) override;

protected:
	vec3d NodeValue(const vec3d& X);

protected:
	double	m_scale;
	mat3d	m_F;
	vector<int>	m_node;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEBCPrescribedDeformation2O : public FEPrescribedBC
{
public:
	FEBCPrescribedDeformation2O(FEModel* pfem);

public:
	void SetReferenceNode(int n);
	void AddNode(int n);
	void AddNodes(const FENodeSet& set) override;

	bool Init() override;

	void Activate() override;

	void Deactivate() override;

	void PrepStep(std::vector<double>& ui, bool brel) override;

	void Update() override;

	int Items() const { return (int)m_node.size(); }

	int NodeID(int index) const { return m_node[index]; }

	void SetDeformationGradient(const mat3d& F);
	void SetDeformationHessian(const tens3drs& G);

	void CopyFrom(FEPrescribedBC* pbc) override;

public:
	void SetScale(double s, int lc = -1);

protected:
	vec3d NodeValue(const vec3d& X1, const vec3d& X);

protected:
	double	m_scale;
	mat3d	m_F;
	tens3drs m_G;
	vector<int>	m_node;
	int	m_refNode;

	DECLARE_FECORE_CLASS();
};
