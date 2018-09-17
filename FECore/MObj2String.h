#pragma once
#include "MathObject.h"
#include <string>
using namespace std;

// this class converts a MathObject to a string
class MObj2String
{
public:
	string Convert(MathObject& o);
	string Convert(MObjList* po);

protected:
	string Convert(MItem* pi);

	string Constant (MConstant* pc);
	string Fraction (MFraction* pc);
	string NamedCt  (MNamedCt*  pc);
	string Variable (MVarRef*   pv);
	string OpNeg    (MNeg*      po);
	string OpAdd    (MAdd*      po);
	string OpSub    (MSub*      po);
	string OpMul    (MMul*      po);
	string OpDiv    (MDiv*      po);
	string OpPow    (MPow*      po);
	string OpEqual  (MEquation* po);
	string OpFnc1D  (MFunc1D*   po);
	string OpFnc2D  (MFunc2D*   po);
	string OpFncND  (MFuncND*   po);
	string OpSFnc   (MSFuncND*  po);
};
