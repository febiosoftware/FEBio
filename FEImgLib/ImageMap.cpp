//#include "stdafx.h"
#include "ImageMap.h"
#include "math.h"

ImageMap::ImageMap(Image& img) : m_img(img)
{
	m_r0 = vec3d(0,0,0);
	m_r1 = vec3d(0,0,0);
}

ImageMap::~ImageMap(void)
{
}

void ImageMap::SetRange(vec3d r0, vec3d r1)
{
	m_r0 = r0;
	m_r1 = r1;
}

ImageMap::POINT ImageMap::map(const vec3d& p)
{
	int nx = m_img.width ();
	int ny = m_img.height();
	int nz = m_img.depth ();

	double x = (p.x - m_r0.x)/(m_r1.x - m_r0.x);
	double y = (p.y - m_r0.y)/(m_r1.y - m_r0.y);
	double z = (nz==1?0:(p.z - m_r0.z)/(m_r1.z - m_r0.z));

	double dx = 1.0 / (double) (nx - 1);
	double dy = 1.0 / (double) (ny - 1);
	double dz = (nz>1?1.0 / (double) (nz - 1):1.0);

	int i = (int) floor(x*(nx - 1));
	int j = (int) floor(y*(ny - 1));
	int k = (int) floor(z*(nz - 1));

	if (i == nx-1) i--;
	if (j == ny-1) j--;
	if ((k == nz-1)&&(nz>1)) k--;

	double r = (x - i*dx)/dx;
	double s = (y - j*dy)/dy;
	double t = (z - k*dz)/dz;

	POINT pt;
	pt.i = i;
	pt.j = j;
	pt.k = k;

	if (nz == 1)
	{
		pt.h[0] = (1-r)*(1-s);
		pt.h[1] = r*(1-s);
		pt.h[2] = r*s;
		pt.h[3] = s*(1-r);
		pt.h[4] = pt.h[5] = pt.h[6] = pt.h[7] = 0.0;
	}
	else
	{
		pt.h[0] = (1-r)*(1-s)*(1-t);
		pt.h[1] = r*(1-s)*(1-t);
		pt.h[2] = r*s*(1-t);
		pt.h[3] = s*(1-r)*(1-t);
		pt.h[4] = (1-r)*(1-s)*t;
		pt.h[5] = r*(1-s)*t;
		pt.h[6] = r*s*t;
		pt.h[7] = s*(1-r)*t;
	}

	return pt;
}

double ImageMap::value(const POINT& p)
{
	int nx = m_img.width ();
	int ny = m_img.height();
	int nz = m_img.depth ();

	if (nz == 1)
	{
		if ((p.i<0) || (p.i >= nx-1)) return 0.0;
		if ((p.j<0) || (p.j >= ny-1)) return 0.0;

		double v[4];
		v[0] = m_img.value(p.i  , p.j  , 0);
		v[1] = m_img.value(p.i+1, p.j  , 0);
		v[2] = m_img.value(p.i+1, p.j+1, 0);
		v[3] = m_img.value(p.i  , p.j+1, 0);

		return (p.h[0]*v[0] + p.h[1]*v[1] + p.h[2]*v[2] + p.h[3]*v[3]);
	}
	else
	{
		if ((p.i<0) || (p.i >= nx-1)) return 0.0;
		if ((p.j<0) || (p.j >= ny-1)) return 0.0;
		if ((p.k<0) || (p.k >= nz-1)) return 0.0;

		double v[8];
		v[0] = m_img.value(p.i  , p.j  , p.k  );
		v[1] = m_img.value(p.i+1, p.j  , p.k  );
		v[2] = m_img.value(p.i+1, p.j+1, p.k  );
		v[3] = m_img.value(p.i  , p.j+1, p.k  );
		v[4] = m_img.value(p.i  , p.j  , p.k+1);
		v[5] = m_img.value(p.i+1, p.j  , p.k+1);
		v[6] = m_img.value(p.i+1, p.j+1, p.k+1);
		v[7] = m_img.value(p.i  , p.j+1, p.k+1);

		return (p.h[0]*v[0] + p.h[1]*v[1] + p.h[2]*v[2] + p.h[3]*v[3] + p.h[4]*v[4] + p.h[5]*v[5] + p.h[6]*v[6] + p.h[7]*v[7]);
	}
}

vec3d ImageMap::gradient(const vec3d& r)
{
	POINT p = map(r);

	if (m_img.depth() == 1)
	{
		// get the x-gradient values
		double gx[4];
		gx[0] = grad_x(p.i  , p.j  , 0);
		gx[1] = grad_x(p.i+1, p.j  , 0);
		gx[2] = grad_x(p.i+1, p.j+1, 0);
		gx[3] = grad_x(p.i  , p.j+1, 0);

		// get the y-gradient values
		double gy[4];
		gy[0] = grad_y(p.i  , p.j  , 0);
		gy[1] = grad_y(p.i+1, p.j  , 0);
		gy[2] = grad_y(p.i+1, p.j+1, 0);
		gy[3] = grad_y(p.i  , p.j+1, 0);

		// get the gradient at point r
		double hx = (p.h[0]*gx[0]+p.h[1]*gx[1]+p.h[2]*gx[2]+p.h[3]*gx[3]);
		double hy = (p.h[0]*gy[0]+p.h[1]*gy[1]+p.h[2]*gy[2]+p.h[3]*gy[3]);

		double Dx = dx();
		double Dy = dy();

		return vec3d(hx/Dx, hy/Dy, 0.0);
	}
	else
	{
		// get the x-gradient values
		double gx[8];
		gx[0] = grad_x(p.i  , p.j  , p.k  );
		gx[1] = grad_x(p.i+1, p.j  , p.k  );
		gx[2] = grad_x(p.i+1, p.j+1, p.k  );
		gx[3] = grad_x(p.i  , p.j+1, p.k  );
		gx[4] = grad_x(p.i  , p.j  , p.k+1);
		gx[5] = grad_x(p.i+1, p.j  , p.k+1);
		gx[6] = grad_x(p.i+1, p.j+1, p.k+1);
		gx[7] = grad_x(p.i  , p.j+1, p.k+1);

		// get the y-gradient values
		double gy[8];
		gy[0] = grad_y(p.i  , p.j  , p.k  );
		gy[1] = grad_y(p.i+1, p.j  , p.k  );
		gy[2] = grad_y(p.i+1, p.j+1, p.k  );
		gy[3] = grad_y(p.i  , p.j+1, p.k  );
		gy[4] = grad_y(p.i  , p.j  , p.k+1);
		gy[5] = grad_y(p.i+1, p.j  , p.k+1);
		gy[6] = grad_y(p.i+1, p.j+1, p.k+1);
		gy[7] = grad_y(p.i  , p.j+1, p.k+1);

		// get the z-gradient values
		double gz[8];
		gz[0] = grad_z(p.i  , p.j  , p.k  );
		gz[1] = grad_z(p.i+1, p.j  , p.k  );
		gz[2] = grad_z(p.i+1, p.j+1, p.k  );
		gz[3] = grad_z(p.i  , p.j+1, p.k  );
		gz[4] = grad_z(p.i  , p.j  , p.k+1);
		gz[5] = grad_z(p.i+1, p.j  , p.k+1);
		gz[6] = grad_z(p.i+1, p.j+1, p.k+1);
		gz[7] = grad_z(p.i  , p.j+1, p.k+1);

		// get the gradient at point r
		double hx = (p.h[0]*gx[0]+p.h[1]*gx[1]+p.h[2]*gx[2]+p.h[3]*gx[3]+p.h[4]*gx[4]+p.h[5]*gx[5]+p.h[6]*gx[6]+p.h[7]*gx[7]);
		double hy = (p.h[0]*gy[0]+p.h[1]*gy[1]+p.h[2]*gy[2]+p.h[3]*gy[3]+p.h[4]*gy[4]+p.h[5]*gy[5]+p.h[6]*gy[6]+p.h[7]*gy[7]);
		double hz = (p.h[0]*gz[0]+p.h[1]*gz[1]+p.h[2]*gz[2]+p.h[3]*gz[3]+p.h[4]*gz[4]+p.h[5]*gz[5]+p.h[6]*gz[6]+p.h[7]*gz[7]);

		double Dx = dx();
		double Dy = dy();
		double Dz = dz();

		return vec3d(hx/Dx, hy/Dy, hz/Dz);
	}
}

double ImageMap::grad_x(int i, int j, int k)
{
	int nx = m_img.width();
	if (i == 0   ) return m_img.value(i+1, j, k) - m_img.value(i  , j, k);
	if (i == nx-1) return m_img.value(i  , j, k) - m_img.value(i-1, j, k);
	return 0.5*(m_img.value(i+1, j, k) - m_img.value(i-1, j, k));
}

double ImageMap::grad_y(int i, int j, int k)
{
	int ny = m_img.height();
	if (j == 0   ) return m_img.value(i, j+1, k) - m_img.value(i, j  , k);
	if (j == ny-1) return m_img.value(i, j  , k) - m_img.value(i, j-1, k);
	return 0.5*(m_img.value(i, j+1, k) - m_img.value(i, j-1, k));
}

double ImageMap::grad_z(int i, int j, int k)
{
	int nz = m_img.depth();
	if (nz == 1) return 0.0;
	if (k == 0   ) return m_img.value(i, j, k+1) - m_img.value(i, j, k  );
	if (k == nz-1) return m_img.value(i, j, k  ) - m_img.value(i, j, k-1);
	return 0.5*(m_img.value(i, j, k+1) - m_img.value(i, j, k-1));
}

mat3ds ImageMap::hessian(const vec3d& r)
{
	POINT p = map(r);

	if (m_img.depth() == 1)
	{
		// get the xx-hessian values
		double hxx[4];
		hxx[0] = hessian_xx(p.i  , p.j  , 0);
		hxx[1] = hessian_xx(p.i+1, p.j  , 0);
		hxx[2] = hessian_xx(p.i+1, p.j+1, 0);
		hxx[3] = hessian_xx(p.i  , p.j+1, 0);

		// get the yy-hessian values
		double hyy[4];
		hyy[0] = hessian_yy(p.i  , p.j  , 0);
		hyy[1] = hessian_yy(p.i+1, p.j  , 0);
		hyy[2] = hessian_yy(p.i+1, p.j+1, 0);
		hyy[3] = hessian_yy(p.i  , p.j+1, 0);

		// get the xy-hessian values
		double hxy[4];
		hxy[0] = hessian_xy(p.i  , p.j  , 0);
		hxy[1] = hessian_xy(p.i+1, p.j  , 0);
		hxy[2] = hessian_xy(p.i+1, p.j+1, 0);
		hxy[3] = hessian_xy(p.i  , p.j+1, 0);

		// get the hessian at point r
		double Hxx = (p.h[0]*hxx[0]+p.h[1]*hxx[1]+p.h[2]*hxx[2]+p.h[3]*hxx[3]);
		double Hyy = (p.h[0]*hyy[0]+p.h[1]*hyy[1]+p.h[2]*hyy[2]+p.h[3]*hyy[3]);
		double Hxy = (p.h[0]*hxy[0]+p.h[1]*hxy[1]+p.h[2]*hxy[2]+p.h[3]*hxy[3]);

		double Dxx = dx()*dx();
		double Dyy = dy()*dy();
		double Dxy = dx()*dy();

		return mat3ds(Hxx/Dxx, Hyy/Dyy, 0, Hxy/Dxy, 0.0, 0.0);
	}
	else
	{
		// get the xx-hessian values
		double hxx[8];
		hxx[0] = hessian_xx(p.i  , p.j  , p.k);
		hxx[1] = hessian_xx(p.i+1, p.j  , p.k);
		hxx[2] = hessian_xx(p.i+1, p.j+1, p.k);
		hxx[3] = hessian_xx(p.i  , p.j+1, p.k);
		hxx[4] = hessian_xx(p.i  , p.j  , p.k+1);
		hxx[5] = hessian_xx(p.i+1, p.j  , p.k+1);
		hxx[6] = hessian_xx(p.i+1, p.j+1, p.k+1);
		hxx[7] = hessian_xx(p.i  , p.j+1, p.k+1);

		// get the yy-hessian values
		double hyy[8];
		hyy[0] = hessian_yy(p.i  , p.j  , p.k);
		hyy[1] = hessian_yy(p.i+1, p.j  , p.k);
		hyy[2] = hessian_yy(p.i+1, p.j+1, p.k);
		hyy[3] = hessian_yy(p.i  , p.j+1, p.k);
		hyy[4] = hessian_yy(p.i  , p.j  , p.k+1);
		hyy[5] = hessian_yy(p.i+1, p.j  , p.k+1);
		hyy[6] = hessian_yy(p.i+1, p.j+1, p.k+1);
		hyy[7] = hessian_yy(p.i  , p.j+1, p.k+1);

		// get the zz-hessian values
		double hzz[8];
		hzz[0] = hessian_zz(p.i  , p.j  , p.k);
		hzz[1] = hessian_zz(p.i+1, p.j  , p.k);
		hzz[2] = hessian_zz(p.i+1, p.j+1, p.k);
		hzz[3] = hessian_zz(p.i  , p.j+1, p.k);
		hzz[4] = hessian_zz(p.i  , p.j  , p.k+1);
		hzz[5] = hessian_zz(p.i+1, p.j  , p.k+1);
		hzz[6] = hessian_zz(p.i+1, p.j+1, p.k+1);
		hzz[7] = hessian_zz(p.i  , p.j+1, p.k+1);

		// get the xy-hessian values
		double hxy[8];
		hxy[0] = hessian_xy(p.i  , p.j  , p.k);
		hxy[1] = hessian_xy(p.i+1, p.j  , p.k);
		hxy[2] = hessian_xy(p.i+1, p.j+1, p.k);
		hxy[3] = hessian_xy(p.i  , p.j+1, p.k);
		hxy[4] = hessian_xy(p.i  , p.j  , p.k+1);
		hxy[5] = hessian_xy(p.i+1, p.j  , p.k+1);
		hxy[6] = hessian_xy(p.i+1, p.j+1, p.k+1);
		hxy[7] = hessian_xy(p.i  , p.j+1, p.k+1);

		// get the yz-hessian values
		double hyz[8];
		hyz[0] = hessian_yz(p.i  , p.j  , p.k);
		hyz[1] = hessian_yz(p.i+1, p.j  , p.k);
		hyz[2] = hessian_yz(p.i+1, p.j+1, p.k);
		hyz[3] = hessian_yz(p.i  , p.j+1, p.k);
		hyz[4] = hessian_yz(p.i  , p.j  , p.k+1);
		hyz[5] = hessian_yz(p.i+1, p.j  , p.k+1);
		hyz[6] = hessian_yz(p.i+1, p.j+1, p.k+1);
		hyz[7] = hessian_yz(p.i  , p.j+1, p.k+1);

		// get the xz-hessian values
		double hxz[8];
		hxz[0] = hessian_xz(p.i  , p.j  , p.k);
		hxz[1] = hessian_xz(p.i+1, p.j  , p.k);
		hxz[2] = hessian_xz(p.i+1, p.j+1, p.k);
		hxz[3] = hessian_xz(p.i  , p.j+1, p.k);
		hxz[4] = hessian_xz(p.i  , p.j  , p.k+1);
		hxz[5] = hessian_xz(p.i+1, p.j  , p.k+1);
		hxz[6] = hessian_xz(p.i+1, p.j+1, p.k+1);
		hxz[7] = hessian_xz(p.i  , p.j+1, p.k+1);

		// get the hessian at point r
		double Hxx = (p.h[0]*hxx[0]+p.h[1]*hxx[1]+p.h[2]*hxx[2]+p.h[3]*hxx[3]+p.h[4]*hxx[4]+p.h[5]*hxx[5]+p.h[6]*hxx[6]+p.h[7]*hxx[7]);
		double Hyy = (p.h[0]*hyy[0]+p.h[1]*hyy[1]+p.h[2]*hyy[2]+p.h[3]*hyy[3]+p.h[4]*hyy[4]+p.h[5]*hyy[5]+p.h[6]*hyy[6]+p.h[7]*hyy[7]);
		double Hzz = (p.h[0]*hzz[0]+p.h[1]*hzz[1]+p.h[2]*hzz[2]+p.h[3]*hzz[3]+p.h[4]*hzz[4]+p.h[5]*hzz[5]+p.h[6]*hzz[6]+p.h[7]*hzz[7]);
		double Hxy = (p.h[0]*hxy[0]+p.h[1]*hxy[1]+p.h[2]*hxy[2]+p.h[3]*hxy[3]+p.h[4]*hxy[4]+p.h[5]*hxy[5]+p.h[6]*hxy[6]+p.h[7]*hxy[7]);
		double Hyz = (p.h[0]*hyz[0]+p.h[1]*hyz[1]+p.h[2]*hyz[2]+p.h[3]*hyz[3]+p.h[4]*hyz[4]+p.h[5]*hyz[5]+p.h[6]*hyz[6]+p.h[7]*hyz[7]);
		double Hxz = (p.h[0]*hxz[0]+p.h[1]*hxz[1]+p.h[2]*hxz[2]+p.h[3]*hxz[3]+p.h[4]*hxz[4]+p.h[5]*hxz[5]+p.h[6]*hxz[6]+p.h[7]*hxz[7]);

		double Dxx = dx()*dx();
		double Dyy = dy()*dy();
		double Dzz = dz()*dz();
		double Dxy = dx()*dy();
		double Dyz = dy()*dz();
		double Dxz = dx()*dz();

		return mat3ds(Hxx/Dxx, Hyy/Dyy, Hzz/Dzz, Hxy/Dxy, Hyz/Dyz, Hxz/Dxz);
	}
}

double ImageMap::hessian_xx(int i, int j, int k)
{
	int nx = m_img.width();
	if (nx <= 2) return 0.0;
	if (i==0   ) return (grad_x(i+1,j,k) - grad_x(i  ,j,k));
	if (i==nx-1) return (grad_x(i  ,j,k) - grad_x(i-1,j,k));
	return 0.5*(grad_x(i+1,j,k) - grad_x(i-1,j,k));
}

double ImageMap::hessian_yy(int i, int j, int k)
{
	int ny = m_img.height();
	if (ny <= 2) return 0.0;
	if (j==0   ) return (grad_y(i,j+1,k) - grad_y(i  ,j,k));
	if (j==ny-1) return (grad_y(i  ,j,k) - grad_y(i,j-1,k));
	return 0.5*(grad_y(i,j+1,k) - grad_y(i,j-1,k));
}

double ImageMap::hessian_zz(int i, int j, int k)
{
	int nz = m_img.depth();
	if (nz <= 2) return 0.0;
	if (k==0   ) return (grad_z(i,j,k+1) - grad_z(i  ,j,k));
	if (k==nz-1) return (grad_z(i  ,j,k) - grad_z(i,j,k-1));
	return 0.5*(grad_z(i,j,k+1) - grad_z(i,j,k-1));
}

double ImageMap::hessian_xy(int i, int j, int k)
{
	int nx = m_img.width();
	if (nx <= 2) return 0.0;
	if (i==0   ) return (grad_y(i+1,j,k) - grad_y(i  ,j,k));
	if (i==nx-1) return (grad_y(i  ,j,k) - grad_y(i-1,j,k));
	return 0.5*(grad_y(i+1,j,k) - grad_y(i-1,j,k));
}

double ImageMap::hessian_yz(int i, int j, int k)
{
	int ny = m_img.height();
	if (ny <= 2) return 0.0;
	if (j==0   ) return (grad_z(i,j+1,k) - grad_z(i  ,j,k));
	if (j==ny-1) return (grad_z(i  ,j,k) - grad_z(i,j-1,k));
	return 0.5*(grad_z(i,j+1,k) - grad_z(i,j-1,k));
}

double ImageMap::hessian_xz(int i, int j, int k)
{
	int nx = m_img.width();
	if (nx <= 2) return 0.0;
	if (i==0   ) return (grad_z(i+1,j,k) - grad_z(i  ,j,k));
	if (i==nx-1) return (grad_z(i  ,j,k) - grad_z(i-1,j,k));
	return 0.5*(grad_z(i+1,j,k) - grad_z(i-1,j,k));
}
