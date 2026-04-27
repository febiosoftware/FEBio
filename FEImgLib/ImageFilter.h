#pragma once
#include "Image.h"
#include <FECore/FECoreClass.h>
#include <FECore/mat3d.h>
#include "feimglib_api.h"
#include "image_tools.h"

#ifdef HAVE_FFTW
#include <fftw3.h>
#endif // HAVE_FFTW

class FEIMGLIB_API ImageFilter : public FECoreClass
{
	FECORE_BASE_CLASS(ImageFilter)

public:
	//! default constructor
	ImageFilter(FEModel* fem);

	//! initialize filter
	virtual bool Init() = 0;

	//! evaluate the filter at the current position
	virtual void Update(Image& trg, Image& src) = 0;

	//! evaluate the filter at the current position
	virtual float Apply(Image& img, int m_pos[3], int m_range[3], int m_dir) = 0;

	//! return blur value
	virtual double GetBlur() = 0;

protected:
};

class FEIMGLIB_API IterativeBlur : public ImageFilter
{
public:
	IterativeBlur(FEModel* fem);

	bool Init() override;

	void Update(Image& trg, Image& src) override;

	float Apply(Image& img, int m_pos[3], int m_range[3], int m_dir) override;

	double GetBlur() override { return m_blur; };

protected:
	//! flag to normalize data so that blur units coincide with physical dimensions rather than img dimensions
	bool m_norm_flag;
	double m_blur;
	DECLARE_FECORE_CLASS();
};

class FEIMGLIB_API BoxBlur : public ImageFilter
{
public:
	BoxBlur(FEModel* fem);

	bool Init() override;

	void Update(Image& trg, Image& src) override;

	float Apply(Image& img, int m_pos[3], int m_range[3], int m_dir) override;

	//! return the current blur value
	double GetBlur() override { return m_blur; };

public:
	bool m_norm_flag = false;
	double m_blur; // sigma (units of length)
	double m_rp; // previous modulus
	double m_tp; // previous time
	int m_K; // iterations of box blur to apply
	double m_res[3]; // voxel resolution (L / px)
	int m_ri[3]; // effective radii along each direction
	DECLARE_FECORE_CLASS();
};

#ifdef HAVE_FFTW
class FEIMGLIB_API FFTWBlur : public ImageFilter
{
public:
	FFTWBlur(FEModel* fem);

	bool Init() override;

	void Update(Image& trg, Image& src) override;

	float Apply(Image& img, int m_pos[3], int m_range[3], int m_dir) override;

	double GetBlur() override { return m_blur; };

public:
	bool m_norm_flag = false;
	double m_blur; // sigma (units of length)
	double m_res[3]; // voxel resolution (L / px)
	double m_sigma[3];
	DECLARE_FECORE_CLASS();
};
#endif // HAVE_FFTW