#include "stdafx.h"
#include "FEDataExport.h"
#include "vec3d.h"

//-----------------------------------------------------------------------------
void FEDataExport::Serialize(vector<float>& d)
{
	if ((m_type == PLT_VEC3F)&&(m_fmt == FMT_NODE))
	{
		vector<vec3d>& v = *(static_cast<vector<vec3d>*>(m_pd));

		int n = (int) v.size();
		for (int i=0; i<n; ++i)
		{
			d.push_back((float) v[i].x);
			d.push_back((float) v[i].y);
			d.push_back((float) v[i].z);
		}
	}
}
