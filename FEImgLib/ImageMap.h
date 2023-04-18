#pragma once
#include "Image.h"
#include <FECore/vec3d.h>
#include <FECore/mat3d.h>
#include "feimglib_api.h"

class FEIMGLIB_API ImageMap
{
public:
	struct POINT
	{
		int		i, j, k;
		double	h[8];
	};

public:
	ImageMap(Image& img);
	~ImageMap(void);

	void SetRange(vec3d r0, vec3d r1);

	// map a vector to the image domain
	POINT map(const vec3d& p);

	// evaluate image
	double value(const POINT& p);
	double value(const vec3d& r) { return value(map(r)); }

	// image gradient
	vec3d gradient(const vec3d& r);

	// image hessian
	mat3ds hessian(const vec3d& r);

	// pixel dimensions
	double dx() { return (m_r1.x - m_r0.x)/(double) (m_img.width () - 1); }
	double dy() { return (m_r1.y - m_r0.y)/(double) (m_img.height() - 1); }
	double dz() { int nz = m_img.depth(); if (nz == 1) return 1.0; else return (m_r1.z - m_r0.z)/(double) (m_img.depth () - 1); }

protected:
	double grad_x(int i, int j, int k);
	double grad_y(int i, int j, int k);
	double grad_z(int i, int j, int k);

	double hessian_xx(int i, int j, int k);
	double hessian_yy(int i, int j, int k);
	double hessian_zz(int i, int j, int k);
	double hessian_xy(int i, int j, int k);
	double hessian_yz(int i, int j, int k);
	double hessian_xz(int i, int j, int k);

protected:
	Image&	m_img;
	vec3d	m_r0;
	vec3d	m_r1;	
};
