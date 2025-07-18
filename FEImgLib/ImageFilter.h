#pragma once
#include "Image.h"
#include <FECore/FECoreClass.h>
#include <FECore/mat3d.h>
#include "feimglib_api.h"

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
	virtual double Apply(Image& img, int m_pos[3], int m_range[3], int m_dir) = 0;

protected:
	//Image& m_img;
	//ImageFilter& m_filt;
};

class FEIMGLIB_API IterativeBlur1D : public ImageFilter
{
public:
	IterativeBlur1D(FEModel* fem);

	bool Init() override;

	void Update(Image& trg, Image& src) override;

	double Apply(Image& img, int m_pos[3], int m_range[3], int m_dir) override;

protected:
	//! flag to normalize data so that blur units coincide with physical dimensions rather than img dimensions
	bool m_norm_flag;
	double m_blur;
	DECLARE_FECORE_CLASS();
};

class FEIMGLIB_API IterativeBlur3D : public ImageFilter
{
public:
	IterativeBlur3D(FEModel* fem);

	bool Init() override;

	void Update(Image& trg, Image& src) override;

	double Apply(Image& img, int m_pos[3], int m_range[3], int m_dir) override;

protected:
	//! flag to normalize data so that blur units coincide with physical dimensions rather than img dimensions
	bool m_norm_flag;
	double m_blur;
	DECLARE_FECORE_CLASS();
};

class FEIMGLIB_API BoxBlur1D : public ImageFilter
{
public:
	BoxBlur1D(FEModel* fem);

	bool Init() override;

	void Update(Image& trg, Image& src) override;

	double Apply(Image& img, int m_pos[3], int m_range[3], int m_dir) override;

public:
	bool m_norm_flag = false;
	double m_blur;
	DECLARE_FECORE_CLASS();
};

class FEIMGLIB_API BoxBlur3D : public ImageFilter
{
public:
	BoxBlur3D(FEModel* fem);

	bool Init() override;

	void Update(Image& trg, Image& src) override;

	double Apply(Image& img, int m_pos[3], int m_range[3], int m_dir) override;

public:
	bool m_norm_flag = false;
	double m_blur;
	DECLARE_FECORE_CLASS();
};