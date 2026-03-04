#include "ast.h"
#include <iostream>

std::ostream& operator << (std::ostream& o, const febcode::Value& v)
{
	if      (isVoid  (v)) return o << "null";
	else if (isInt   (v)) return o << std::get<int>(v);
	else if (isDouble(v)) return o << std::get<double>(v);
	else if (isBool  (v)) return o << (std::get<bool>(v) ? "true" : "false");
	else if (isString(v)) return o << std::get<std::string>(v);
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
	else if (isString(v)) s = "\"" + getString(v) + "\"";
	else if (isArray (v)) { auto& p = getArrayPtr (v); s = "[array:"  + std::to_string(p.use_count()) + "]"; }
	else if (isStruct(v)) { 
		const febcode::StructValue& o = febcode::getStruct(v);
		const std::string& name = o.type->name;
		auto& p = getStructPtr(v);
		s = "[" + name + ":" + std::to_string(p.use_count()) + "]";
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
	else return "unknown type";
}
