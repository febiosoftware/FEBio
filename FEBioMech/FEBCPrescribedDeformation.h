#include <FECore/BC.h>
#include <FECore/tens3d.h>

//-----------------------------------------------------------------------------
class FEBCPrescribedDeformation : public FEPrescribedBC
{
public:
	FEBCPrescribedDeformation(FEModel* pfem);

public:
	void AddNode(int n);
	void AddNodes(const FENodeSet& set);

	void Activate();

	void Deactivate();

	void PrepStep(std::vector<double>& ui, bool brel);

	void Update();

	int Items() const { return (int) m_node.size(); }

	int NodeID(int index) const { return m_node[index]; }

	void SetDeformationGradient(const mat3d& F);
	void SetDeformationHessian(const tens3drs& G);

	void CopyFrom(FEPrescribedBC* pbc);

protected:
	vec3d NodeValue(const vec3d& X);

protected:
	double	m_scale;
	mat3d	m_F;
	tens3drs m_G;
	vector<int>	m_node;

	DECLARE_PARAMETER_LIST();
};
