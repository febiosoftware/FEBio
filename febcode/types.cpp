#include "types.h"
using namespace febcode;


ArrayValuePtr febcode::copyArray(const ArrayValue& src)
{
	auto arr = std::make_shared<ArrayValue>();
	arr->type = src.type;
	for (int i = 0; i < src.size(); ++i)
	{
		if (isStruct(src.elements[i]))
		{
			arr->elements.emplace_back(copyStruct(getStruct(src.elements[i])));
		}
		else if (isArray(src.elements[i]))
		{
			arr->elements.emplace_back(copyArray(getArray(src.elements[i])));
		}
		else
			arr->elements.emplace_back(src.elements[i]);
	}
	return arr;
}

StructValuePtr febcode::copyStruct(const StructValue& src)
{
	auto obj = std::make_shared<StructValue>();
	obj->type = src.type;
	for (size_t i = 0; i < src.fields.size(); ++i)
	{
		if (isStruct(src.fields[i]))
		{
			obj->fields.emplace_back(copyStruct(getStruct(src.fields[i])));
		}
		else if (isArray(src.fields[i]))
		{
			obj->fields.emplace_back(copyArray(getArray(src.fields[i])));
		}
		else
			obj->fields.emplace_back(src.fields[i]);
	}

	return obj;
}

TypeRegistry::TypeRegistry()
{
	m_void   = { TypeKind::Void  };
	m_bool   = { TypeKind::Bool  };
	m_int    = { TypeKind::Int   };
	m_double = { TypeKind::Double};
	m_vec2   = { TypeKind::Vec2  };
	m_vec3   = { TypeKind::Vec3  };
	m_mat2   = { TypeKind::Mat2  };
	m_mat3   = { TypeKind::Mat3  };
}

void TypeRegistry::clear()
{
//	m_arrayTypes.clear();
}

Type TypeRegistry::getVoidType() const { return &m_void; }
Type TypeRegistry::getBoolType() const { return &m_bool; }
Type TypeRegistry::getIntType() const { return &m_int; }
Type TypeRegistry::getDoubleType() const { return &m_double; }
Type TypeRegistry::getVec2Type() const { return &m_vec2; }
Type TypeRegistry::getVec3Type() const { return &m_vec3; }
Type TypeRegistry::getMat2Type() const { return &m_mat2; }
Type TypeRegistry::getMat3Type() const { return &m_mat3; }

Type TypeRegistry::getTypeFromKind(TypeKind kind) const
{
	switch (kind)
	{
	case TypeKind::Void  : return getVoidType();
	case TypeKind::Bool  : return getBoolType();
	case TypeKind::Int   : return getIntType();
	case TypeKind::Double: return getDoubleType();
	case TypeKind::Vec2  : return getVec2Type();
	case TypeKind::Vec3  : return getVec3Type();
	case TypeKind::Mat2  : return getMat2Type();
	case TypeKind::Mat3  : return getMat3Type();
	default:
		return nullptr;
	}
}

Type TypeRegistry::getArrayType(Type element, size_t size)
{
	if (size == 0)
		throw std::runtime_error("Array size must be greater than zero.");

	// Check if already exists
	for (auto& t : m_arrayTypes)
	{
		if ((t->elementType == element) && (t->arraySize == size))
			return t.get();
	}

	// Create new one
	auto newType = std::make_unique<TypeStruct>();
	newType->kind = TypeKind::Array;
	newType->typeIndex = (int)m_arrayTypes.size();
	newType->elementType = element;
	newType->arraySize = size;

	m_arrayTypes.push_back(std::move(newType));
	return m_arrayTypes.back().get();
}

Type TypeRegistry::getArrayType(Type element, const std::vector<size_t>& sizes) const
{
	Type type = element;
	for (auto it = sizes.rbegin(); it != sizes.rend(); ++it)
	{
		type = getArrayType(type, *it);
	}
	return type;
}

Type TypeRegistry::getArrayType(int index) const
{
	if (index < 0 || index >= (int)m_arrayTypes.size())
		throw std::runtime_error("Invalid array type index: " + std::to_string(index));
	return m_arrayTypes[index].get();
}

Type TypeRegistry::getArrayType(Type element, size_t size) const
{
	if (size == 0)
		throw std::runtime_error("Array size must be greater than zero.");

	for (auto& t : m_arrayTypes)
	{
		if ((t->elementType == element) && (t->arraySize == size))
			return t.get();
	}

	return nullptr;
}

Type TypeRegistry::getBuiltinType(const Value& v) const
{
	if (isVoid  (v)) return getVoidType();
	if (isBool  (v)) return getBoolType();
	if (isInt   (v)) return getIntType();
	if (isDouble(v)) return getDoubleType();
	if (isVec2  (v)) return getVec2Type();
	if (isVec3  (v)) return getVec3Type();
	if (isMat2  (v)) return getMat2Type();
	if (isMat3  (v)) return getMat3Type();
	if (isArray (v)) { const ArrayValue&  a = getArray (v); return a.type; }
	if (isStruct(v)) { const StructValue& s = getStruct(v); return s.type; }

	throw std::runtime_error("Unknown value type");
	return getVoidType();
}

Type TypeRegistry::getType(const std::string& name) const
{
	if (name == "void"  ) return getVoidType();
	if (name == "bool"  ) return getBoolType();
	if (name == "int"   ) return getIntType();
	if (name == "double") return getDoubleType();
	if (name == "vec2"  ) return getVec2Type();
	if (name == "vec3"  ) return getVec3Type();
	if (name == "mat2"  ) return getMat2Type();
	if (name == "mat3"  ) return getMat3Type();

	Type s = getStructType(name);
	if (s) return s;

	throw std::runtime_error("Unknown type name: " + name);
}

Type TypeRegistry::getStructType(const std::string& name) const
{
	auto it = std::find_if(m_structTypes.begin(), m_structTypes.end(),
		[&name](const std::unique_ptr<TypeStruct>& t) { return t->name == name; });
	if (it == m_structTypes.end()) return nullptr;
	return it->get();
}

int TypeRegistry::getStructTypeIndex(const std::string& name) const
{
	auto it = std::find_if(m_structTypes.begin(), m_structTypes.end(),
		[&name](const std::unique_ptr<TypeStruct>& t) { return t->name == name; });
	if (it == m_structTypes.end()) return -1;
	return (*it)->typeIndex;
}

Type TypeRegistry::getStructType(int index) const
{
	if (index < 0 || index >= (int)m_structTypes.size())
		throw std::runtime_error("Invalid struct type index: " + std::to_string(index));
	return m_structTypes[index].get();
}

Type TypeRegistry::defineStructType(
	const std::string& name,
	const std::vector<std::pair<Type, std::string>>& fields)
{
	Type s = getStructType(name);
	if (s)
		throw std::runtime_error("Struct already defined: " + name);

	auto t = std::make_unique<TypeStruct>();
	t->kind = TypeKind::Struct;
	t->typeIndex = (int)m_structTypes.size();
	t->name = name;
	t->fields = fields;

	Type result = t.get();
	m_structTypes.push_back(std::move(t));

	return result;
}

Type TypeRegistry::defineStructType(const std::string& name, const std::vector<std::pair<TypeKind, std::string>>& fields)
{
	std::vector<std::pair<Type, std::string>> resolvedFields;
	for (const auto& f : fields)
	{
		Type t = nullptr;
		switch (f.first)
		{
		case TypeKind::Bool  : t = getBoolType(); break;
		case TypeKind::Int   : t = getIntType(); break;
		case TypeKind::Double: t = getDoubleType(); break;
		case TypeKind::Vec2  : t = getVec2Type(); break;
		case TypeKind::Vec3  : t = getVec3Type(); break;
		case TypeKind::Mat2  : t = getMat2Type(); break;
		case TypeKind::Mat3  : t = getMat3Type(); break;
		default:
			throw std::runtime_error("Unsupported field type in RegisterStruct: " + std::to_string((int)f.first));
		}
		resolvedFields.push_back({ t, f.second });
	}
	return defineStructType(name, resolvedFields);
}

std::string febcode::TypeToString(Type type)
{
	if (type == nullptr) return "<null>";

	switch (type->kind)
	{
	case TypeKind::Void  : return "void";
	case TypeKind::Bool  : return "bool";
	case TypeKind::Int   : return "int";
	case TypeKind::Double: return "double";
	case TypeKind::Vec2  : return "vec2";
	case TypeKind::Vec3  : return "vec3";
	case TypeKind::Mat2  : return "mat2";
	case TypeKind::Mat3  : return "mat3";
	case TypeKind::Struct: return type->name;
	case TypeKind::Array:
	{
		std::string elemStr = TypeToString(type->elementType);
		return "(" + elemStr + ")[" + std::to_string(type->arraySize) + "]";
	}
	default:
		return "<unknown type>";
	}
}
