#include <FECore/BC.h>
#include <FECore/tens3d.h>

//-----------------------------------------------------------------------------
class FEBCPrescribedDeformation : public FEPrescribedBC
{
public:
	FEBCPrescribedDeformation(FEModel* pfem);

public:
	void AddNodes(const FENodeSet& set);

	void Activate();

	void Deactivate();

	void PrepStep(std::vector<double>& ui, bool brel);

	void Update();

protected:
	vec3d NodeValue(const vec3d& X);

protected:
	double	m_scale;
	mat3d	m_F;
	tens3drs m_G;
	vector<int>	m_node;

	DECLARE_PARAMETER_LIST();
};
