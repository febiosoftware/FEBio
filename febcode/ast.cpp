#include "ast.h"
#include <iostream>

std::ostream& operator << (std::ostream& o, const febcode::Value& v)
{
	if      (isVoid  (v)) return o << "null";
	else if (isInt   (v)) return o << getInt(v);
	else if (isDouble(v)) return o << getDouble(v);
	else if (isBool  (v)) return o << (getBool(v) ? "true" : "false");
	else if (isString(v)) return o << getString(v);
	else if (isArray(v))
	{
		const febcode::ArrayValue& arr = getArray(v);
		o << "[";
		for (size_t i = 0; i < arr.size(); ++i)
		{
			o << arr.elements[i];
			if (i != arr.size() - 1) o << ", ";
		}
		return o << "]";
	}
	else if (isStruct(v))
	{
		const febcode::StructValue& s = getStruct(v);
		o << "{";
		for (size_t i = 0; i < s.fields.size(); ++i)
		{
			o << s.fields[i];
			if (i != s.fields.size() - 1) o << ", ";
		}
		return o << "}";
	}
	else if (isVec2(v))
	{
		const febcode::vec2& vec = getVec2(v);
		return o << "vec2(" << vec.x << ", " << vec.y << ")";
	}
	else if (isVec3(v))
	{
		const febcode::vec3& vec = getVec3(v);
		return o << "vec3(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
	}
	else if (isRef(v))
	{
		const febcode::Ref& ref = v.ref;
		return o << "ref";
	}
	else
		return o << "<unknown value>";
}

std::string ValueToString(const febcode::Value& v)
{
	std::string s;
	if      (isVoid  (v)) s = "null";
	else if (isBool  (v)) s = getBool(v) ? "true" : "false";
	else if (isInt   (v)) s = std::to_string(getInt (v));
	else if (isDouble(v)) s = std::to_string(getDouble(v));
	else if (isVec2  (v)) s = "vec2";
	else if (isVec3  (v)) s = "vec3";
	else if (isString(v)) s = "\"" + getString(v) + "\"";
	else if (isArray (v)) { auto& p = getArrayPtr (v); s = "[array:"  + std::to_string(p.use_count()) + "]"; }
	else if (isStruct(v)) { 
		const febcode::StructValue& o = febcode::getStruct(v);
		const std::string& name = o.type->name;
		auto& p = getStructPtr(v);
		s = "[" + name;
		s = "{";
		for (size_t i = 0; i < o.fields.size(); ++i)
		{
			s += ValueToString(o.fields[i]);
			if (i != o.fields.size() - 1) s += ", ";
		}
		s += "}";
		s += ":" + std::to_string(p.use_count()) + "]";
	}
	else if (isRef(v))
	{
		const febcode::Ref& ref = getRef(v);
		s = "ref->";
		switch (ref.type)
		{
		case febcode::TypeKind::Value:
		{
			febcode::Value* v = static_cast<febcode::Value*>(ref.ptr);
			if (v)
			{
				s += ValueTypeToString(*v);
			}
			else s += "null"; 
			break;
		}
		case febcode::TypeKind::Double: s += "double"; break;
		}
		
	}
	else s = "unknown";
	return s;
}

std::string ValueTypeToString(const febcode::Value& v)
{
	if      (isVoid  (v)) return "void";
	else if (isBool  (v)) return "bool";
	else if (isInt   (v)) return "int";
	else if (isDouble(v)) return "double";
	else if (isString(v)) return "string";
	else if (isArray (v)) return "array";
	else if (isStruct(v)) return "struct";
	else if (isVec2  (v)) return "vec2";
	else if (isVec3  (v)) return "vec3";
	else if (isRef    (v)) return "ref";
	else return "<unknown type>";
}
