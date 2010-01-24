#include "stdafx.h"
#include "FEPlotData.h"
#include "fem.h"

//-----------------------------------------------------------------------------
void FEPlotNodeDisplacement::Save(FEM& fem, Archive& ar)
{
	float xf[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_rt.x;
		xf[1] = (float) node.m_rt.y;
		xf[2] = (float) node.m_rt.z;

		ar.write(xf, sizeof(float), 3);
	}
}

//-----------------------------------------------------------------------------
void FEPlotNodeVelocity::Save(FEM& fem, Archive& ar)
{
	float xf[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_vt.x;
		xf[1] = (float) node.m_vt.y;
		xf[2] = (float) node.m_vt.z;

		ar.write(xf, sizeof(float), 3);
	}
}

//-----------------------------------------------------------------------------
void FEPlotNodeAcceleration::Save(FEM& fem, Archive& ar)
{
	float xf[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_at.x;
		xf[1] = (float) node.m_at.y;
		xf[2] = (float) node.m_at.z;

		ar.write(xf, sizeof(float), 3);
	}
}

//-----------------------------------------------------------------------------
void FEPlotElementStress::Save(FEM& fem, Archive& ar)
{
	int i, j;

	FEMesh& mesh = fem.m_mesh;

	// write solid element data
	float s[6] = {0};
	double f;
	int nint;
	for (i=0; i<mesh.SolidElements(); ++i)
	{
		FESolidElement& el = mesh.SolidElement(i);

		for (j=0; j<6; ++j) s[j] = 0;

		nint = el.GaussPoints();

		f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint* ppt = (el.m_State[j]->ExtractData<FEElasticMaterialPoint>());

			if (ppt)
			{
				FEElasticMaterialPoint& pt = *ppt;
				s[0] += (float) (f*pt.s.xx());
				s[1] += (float) (f*pt.s.yy());
				s[2] += (float) (f*pt.s.zz());
				s[3] += (float) (f*pt.s.xy());
				s[4] += (float) (f*pt.s.yz());
				s[5] += (float) (f*pt.s.xz());
			}
		}

		ar.write(s, sizeof(float), 6);
	}
}
