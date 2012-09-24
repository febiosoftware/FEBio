#pragma once

//-----------------------------------------------------------------------------
// This class implements a simple single precision grayscale 3D-image
class Image
{
public:
	// constructor
	Image(void);

	// destructor
	~Image(void);

	// allocate storage for image data
	void Create(int nx, int ny, int nz);

	// load raw data from file
	bool Load(const char* szfile);

	// return size attributes
	int width () { return m_nx; }
	int height() { return m_ny; }
	int depth () { return m_nz; }

	// get a particular data value
	float& value(int x, int y, int z) { return m_pf[(z*m_ny + (m_ny-y-1))*m_nx+x]; }

	// zero image data
	void zero();

protected:
	float*	m_pf;				// image data
	int		m_nx, m_ny, m_nz;	// image dimensions
};

//-----------------------------------------------------------------------------
// helper functions for calculating image derivatives
void image_derive_x(Image& s, Image& d);
void image_derive_y(Image& s, Image& d);
void image_derive_z(Image& s, Image& d);
