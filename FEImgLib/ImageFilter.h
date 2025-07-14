#pragma once
#include "Image.h"
#include <FECore/mat3d.h>
#include "feimglib_api.h"

class FEIMGLIB_API ImageFilter
{
public:
	//! default constructor
	ImageFilter();

	//!destructor
	~ImageFilter(void) {};

	//! evaluate the filter at the current position
	virtual void eval1D(Image& trg, Image& src, float d) = 0;

	//! evaluate the filter at the current position
	virtual void eval3D(Image& trg, Image& src, float d) = 0;

	//! evaluate the filter at the current position
	virtual double apply1D(Image& img, int m_pos[3], int m_range[3], int m_dir) = 0;

	//! evaluate the filter at the current position
	virtual double apply3D(Image& img, int m_pos[3], int m_range[3]) = 0;

protected:
	//Image& m_img;
	//ImageFilter& m_filt;
};

class IterativeBlur : public ImageFilter
{
public:
	IterativeBlur();

	void eval1D(Image& trg, Image& src, float d) override;

	void eval3D(Image& trg, Image& src, float d) override;

	double apply1D(Image& img, int m_pos[3], int m_range[3], int m_dir) override;

	double apply3D(Image& img, int m_pos[3], int m_range[3]) override;

public:
	//! flag to normalize data so that blur units coincide with physical dimensions rather than img dimensions
	bool m_norm_flag = false;
};

class BoxBlur : public ImageFilter
{
public:
	BoxBlur();

	void eval1D(Image& trg, Image& src, float d) override;

	void eval3D(Image& trg, Image& src, float d) override;

	double apply1D(Image& img, int m_pos[3], int m_range[3], int m_dir) override;

	double apply3D(Image& img, int m_pos[3], int m_range[3]) override;

public:
	bool m_norm_flag = false;
};