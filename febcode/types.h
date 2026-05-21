#pragma once
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <functional>
#include <cstring>
#include <cmath>
#include <assert.h>

#include <FECore/vec2d.h>
#include <FECore/vec3d.h>
#include <FECore/mat2d.h>
#include <FECore/mat3d.h>

namespace febcode
{
	enum class TypeKind : uint8_t {
		Void,
		Bool,
		Int,
		Double,
		Vec2,
		Vec3,
		Mat2,
		Mat3,
		Array,
		Struct,
	};

	using vec2 = vec2d;
	using vec3 = vec3d;
	using mat2 = mat2d;
	using mat3 = mat3d;

	struct TypeStruct;
	using Type = const TypeStruct*;

	using StructField = std::pair<Type, std::string>; // type and name of a struct field

	struct TypeStruct {
		TypeKind kind = TypeKind::Void;

		// index in type registry (for arrays and structs)
		int typeIndex = -1;

		// for arrays
		Type elementType = nullptr;
		size_t arraySize = 0;

		// For struct types
		std::string name;
		std::vector<StructField> fields;

		size_t size() const {
			switch (kind)
			{
			case TypeKind::Void:   return 1;
			case TypeKind::Bool:   return 1;
			case TypeKind::Int:    return 1;
			case TypeKind::Double: return 1;
			case TypeKind::Vec2:   return 2;
			case TypeKind::Vec3:   return 3;
			case TypeKind::Mat2:   return 4;
			case TypeKind::Mat3:   return 9;
			case TypeKind::Array:  return elementType->size()* arraySize;
			case TypeKind::Struct:
				{
					size_t totalSize = 0;
					for (const auto& field : fields)
						totalSize += field.first->size();
					return totalSize;
				}
			default:
				throw std::runtime_error("Unknown type kind.");
			}
		}
	};

	struct ArrayValue;
	struct StructValue;

	using ArrayValuePtr = std::shared_ptr<ArrayValue>;
	using StructValuePtr = std::shared_ptr<StructValue>;

	// make sure this matches the order of the Value variant below (and of TypeKind above)
	enum ValueIndex {
		VOID = 0,
		BOOL,
		INT,
		DOUBLE,
		VEC2,
		VEC3,
		MAT2,
		MAT3,
		ARRAY,
		STRUCT,
		REF
	};

	struct Void {};

	struct Value
	{
		ValueIndex index = ValueIndex::VOID;

		Value() {}
		Value(bool   a) : index(ValueIndex::BOOL), b(a) {}
		Value(int    a) : index(ValueIndex::INT), i(a) {}
		Value(double a) : index(ValueIndex::DOUBLE), d(a) {}
		Value(const vec2& a) : index(ValueIndex::VEC2), vec2Value(a) {}
		Value(const vec3& a) : index(ValueIndex::VEC3), vec3Value(a) {}
		Value(const mat2& a) : index(ValueIndex::MAT2), mat2Value(a) {}
		Value(const mat3& a) : index(ValueIndex::MAT3), mat3Value(a) {}
		Value(const ArrayValuePtr& a) : index(ValueIndex::ARRAY), arrayValue(a) {}
		Value(const StructValuePtr& a) : index(ValueIndex::STRUCT), structValue(a) {}

		Value(const Value& v)
		{
			copyFrom(v);
		}

		void operator = (const Value& other)
		{
			// for simple types we use the fast path
			if ((index < ValueIndex::ARRAY) && (other.index < ValueIndex::ARRAY))
			{
				index = other.index;
				switch (index)
				{
				case ValueIndex::VOID  : break;	// no data to copy
				case ValueIndex::BOOL  : b = other.b; break;
				case ValueIndex::INT   : i = other.i; break;
				case ValueIndex::DOUBLE: d = other.d; break;
				case ValueIndex::VEC2  : vec2Value = other.vec2Value; break;
				case ValueIndex::VEC3  : vec3Value = other.vec3Value; break;
				case ValueIndex::MAT2  : mat2Value = other.mat2Value; break;
				case ValueIndex::MAT3  : mat3Value = other.mat3Value; break;
				default: 
					assert(false);
					break; // should not happen
				}
			}
			else if (this != &other)
			{
				destroy();
				copyFrom(other);
			}
		}

		~Value()
		{
			if (index >= ValueIndex::ARRAY)
				destroy();
		}

		bool operator == (const Value& other) const
		{
			if (index != other.index)
				return false;
			switch (index)
			{
			case ValueIndex::VOID  : return true; // all void values are equal
			case ValueIndex::BOOL  : return b == other.b;
			case ValueIndex::INT   : return i == other.i;
			case ValueIndex::DOUBLE: return d == other.d;
			case ValueIndex::VEC2  : return vec2Value == other.vec2Value;
			case ValueIndex::VEC3  : return vec3Value == other.vec3Value;
			case ValueIndex::MAT2  : return mat2Value == other.mat2Value;
			case ValueIndex::MAT3  : return mat3Value == other.mat3Value;
			case ValueIndex::ARRAY : return arrayValue == other.arrayValue;
			case ValueIndex::STRUCT: return structValue == other.structValue;
			}
			return false;
		}

		bool operator != (const Value& other) const
		{
			return !(*this == other);
		}

		union {
			Void v;
			bool b;
			int i;
			double d;
			vec2 vec2Value;
			vec3 vec3Value;
			mat2 mat2Value;
			mat3 mat3Value;
			ArrayValuePtr arrayValue;
			StructValuePtr structValue;
		};

	private:
		void copyFrom(const Value& v) 
		{
			index = v.index;
			switch (index)
			{
			case ValueIndex::VOID: break;	// no data to copy
			case ValueIndex::BOOL: b = v.b; break;
			case ValueIndex::INT: i = v.i; break;
			case ValueIndex::DOUBLE: d = v.d; break;
			case ValueIndex::VEC2: vec2Value = v.vec2Value; break;
			case ValueIndex::VEC3: vec3Value = v.vec3Value; break;
			case ValueIndex::MAT2: mat2Value = v.mat2Value; break;
			case ValueIndex::MAT3: mat3Value = v.mat3Value; break;
			case ValueIndex::ARRAY: new (&arrayValue) ArrayValuePtr(v.arrayValue); break;
			case ValueIndex::STRUCT: new (&structValue) StructValuePtr(v.structValue); break;
			}
		}

		void destroy()
		{
			switch (index)
			{
				case ValueIndex::ARRAY:
					arrayValue.~shared_ptr();
					break;
				case ValueIndex::STRUCT:
					structValue.~shared_ptr();
					break;
				default:
					break; // no special handling needed for other types
			}
			index = ValueIndex::VOID; // reset to void after destruction
		}
	};

	struct ArrayValue {
		Type type;
		std::vector<Value> elements;
		size_t size() const { return elements.size(); }
	};

	struct StructValue {
		Type type;
		std::vector<Value> fields;
	};

	inline ArrayValuePtr createArray(Type elementType, size_t size)
	{
		auto arr = std::make_shared<ArrayValue>();
		arr->type = elementType;
		arr->elements.resize(size);
		return arr;
	}

	inline StructValuePtr createStruct(Type type)
	{
		auto obj = std::make_shared<StructValue>();
		obj->type = type;
		obj->fields.resize(type->fields.size());
		return obj;
	}

	ArrayValuePtr copyArray(const ArrayValue& src);

	StructValuePtr copyStruct(const StructValue& src);

	inline bool isVoid  (const Value& v) { return v.index == ValueIndex::VOID;}
	inline bool isBool  (const Value& v) { return v.index == ValueIndex::BOOL;}
	inline bool isInt   (const Value& v) { return v.index == ValueIndex::INT; }
	inline bool isDouble(const Value& v) { return v.index == ValueIndex::DOUBLE; }
	inline bool isVec2  (const Value& v) { return v.index == ValueIndex::VEC2; }
	inline bool isVec3  (const Value& v) { return v.index == ValueIndex::VEC3; }
	inline bool isMat2  (const Value& v) { return v.index == ValueIndex::MAT2; }
	inline bool isMat3  (const Value& v) { return v.index == ValueIndex::MAT3; }
	inline bool isArray (const Value& v) { return v.index == ValueIndex::ARRAY; }
	inline bool isStruct(const Value& v) { return v.index == ValueIndex::STRUCT; }
	inline bool isRef   (const Value& v) { return v.index == ValueIndex::REF; }

	inline bool isVoidType  (Type type) { return type->kind == TypeKind::Void; }
	inline bool isBoolType  (Type type) { return type->kind == TypeKind::Bool; }
	inline bool isIntType   (Type type) { return type->kind == TypeKind::Int; }
	inline bool isDoubleType(Type type) { return type->kind == TypeKind::Double; }
	inline bool isVec2Type  (Type type) { return type->kind == TypeKind::Vec2; }
	inline bool isVec3Type  (Type type) { return type->kind == TypeKind::Vec3; }
	inline bool isMat2Type  (Type type) { return type->kind == TypeKind::Mat2; }
	inline bool isMat3Type  (Type type) { return type->kind == TypeKind::Mat3; }
	inline bool isArrayType (Type type) { return type->kind == TypeKind::Array; }
	inline bool isStructType(Type type) { return type->kind == TypeKind::Struct; }

	inline bool isLogicalType(Type type) { return isBoolType(type) || isIntType(type) || isDoubleType(type); }

	inline bool isScalarType(Type type) { return isIntType(type) || isDoubleType(type); }

	inline bool isNumericType(Type type) { return isIntType(type) || isDoubleType(type) || isVec2Type(type) || isVec3Type(type) || isMat2Type(type) || isMat3Type(type); }

	inline bool   getBool(const Value& v) { assert(v.index == ValueIndex::BOOL); return v.b; }
	inline int    getInt   (const Value& v) { assert(v.index == ValueIndex::INT); return v.i; }
	inline double getDouble(const Value& v) { assert(v.index == ValueIndex::DOUBLE); return v.d; }
	inline const vec2& getVec2(const Value& v) { assert(v.index == ValueIndex::VEC2); return v.vec2Value; }
	inline vec2& getVec2(Value& v) { assert(v.index == ValueIndex::VEC2); return v.vec2Value; }
	inline const vec3& getVec3(const Value& v) { assert(v.index == ValueIndex::VEC3); return v.vec3Value; }
	inline vec3& getVec3(Value& v) { assert(v.index == ValueIndex::VEC3); return v.vec3Value; }
	inline const mat2& getMat2(const Value& v) { assert(v.index == ValueIndex::MAT2); return v.mat2Value; }
	inline mat2& getMat2(Value& v) { assert(v.index == ValueIndex::MAT2); return v.mat2Value; }
	inline const mat3& getMat3(const Value& v) { assert(v.index == ValueIndex::MAT3); return v.mat3Value; }
	inline mat3& getMat3(Value& v) { assert(v.index == ValueIndex::MAT3); return v.mat3Value; }
	inline const ArrayValue&  getArray (const Value& v) { assert(v.index == ValueIndex::ARRAY); return *v.arrayValue; }
	inline const StructValue& getStruct(const Value& v) { assert(v.index == ValueIndex::STRUCT); return *v.structValue; }

	inline const ArrayValuePtr& getArrayPtr(const Value& v) { assert(v.index == ValueIndex::ARRAY); return v.arrayValue; }
	inline const StructValuePtr& getStructPtr(const Value& v) { assert(v.index == ValueIndex::STRUCT); return v.structValue; }

	inline ArrayValue&  getArray (Value& v) { assert(v.index == ValueIndex::ARRAY); return *v.arrayValue; }
	inline StructValue& getStruct(Value& v) { assert(v.index == ValueIndex::STRUCT); return *v.structValue; }

	Type coerce(Type from, Type to);

	inline Type commonType(Type l, Type r)
	{
		if (l == r) return l;
		if (isDoubleType(l) && isScalarType(r)) return l;
		if (isDoubleType(r) && isScalarType(l)) return r;
		return nullptr;
	}

	inline bool isZero(const Value& v)
	{
		if (isInt(v) && (getInt(v) == 0)) return true;
		if (isDouble(v) && (getDouble(v) == 0.0)) return true;
		if (isVec2(v) && (getVec2(v) == vec2(0, 0))) return true;
		if (isVec3(v) && (getVec3(v) == vec3(0, 0, 0))) return true;
		if (isMat2(v) && (getMat2(v) == mat2(0.0))) return true;
		if (isMat3(v) && (getMat3(v) == mat3(0.0))) return true;
		return false;
	}

	inline bool isOne(const Value& v)
	{
		if (isInt(v) && (getInt(v) == 1)) return true;
		if (isDouble(v) && (getDouble(v) == 1.0)) return true;
		return false;
	}

	inline bool isIdentity(const mat2& m)
	{
		return (m == mat2(1.0));
	}

	inline bool isIdentity(const mat3& m)
	{
		return (m == mat3(1.0));
	}

	inline bool isSymmetric(const mat2& m)
	{
		return (m(0,1) == m(1,0));
	}

	inline bool isSymmetric(const mat3& m)
	{
		return (m(0,1) == m(1,0)) && (m(0,2) == m(2,0)) && (m(1,2) == m(2,1));
	}

	inline bool isIntNumber(const Value& v)
	{
		if (isInt(v)) return true;
		if (isDouble(v))
		{
			double d = getDouble(v);
			return (d == std::floor(d));
		}
		return false;
	}

	inline int toIntNumber(const Value& v)
	{
		assert(isIntNumber(v));
		if (isInt(v)) return getInt(v);
		if (isDouble(v)) return (int)getDouble(v);
		return 0;
	}

	class TypeRegistry {
	public:
		TypeRegistry();
		void clear();

		Type Void() const;
		Type Int() const;
		Type Bool() const;
		Type Double() const;
		Type Vec2() const;
		Type Vec3() const;
		Type Mat2() const;
		Type Mat3() const;
		Type getArrayType(TypeKind type, size_t size);
		Type getArrayType(Type element, size_t size);
		Type getArrayType(Type element, size_t size) const;
		Type getArrayType(Type element, const std::vector<size_t>& size) const;
		Type getArrayType(int index) const;

		Type getTypeFromKind(TypeKind kind) const;

		Type getBuiltinType(const Value& v) const;

		Type getType(const std::string& name) const;

		Type getStructType(const std::string& name) const;
		int getStructTypeIndex(const std::string& name) const;
		Type getStructType(int index) const;

		Type defineStructType(const std::string& name, const std::vector<std::pair<Type, std::string>>& fields);
		Type defineStructType(const std::string& name, const std::vector<std::pair<TypeKind, std::string>>& fields);

	private:
		// Primitive canonical types
		TypeStruct m_void;
		TypeStruct m_bool;
		TypeStruct m_int;
		TypeStruct m_double;
		TypeStruct m_vec2;
		TypeStruct m_vec3;
		TypeStruct m_mat2;
		TypeStruct m_mat3;

		// Interned array types
		std::vector<std::unique_ptr<TypeStruct>> m_arrayTypes;

		// User-defined struct types
		std::vector<std::unique_ptr<TypeStruct>> m_structTypes;
	};

	inline TypeKind ValueType(const Value& v)
	{
		return (TypeKind)v.index;
	}

	const char* TypeKindToString(TypeKind kind);
	std::string TypeToString(Type type);

	struct BinaryOpSignature
	{
		Type leftType;
		Type rightType;

		Type resultType;
	};
}

namespace febcode
{
	struct FuncArgs {
		double* stack = nullptr;
		size_t count = 0; // number of arguments passed to the function
		size_t index = 0;

		bool   getBool()   { return (bool  )(stack[index++]); }
		int    getInt()    { return (int   )(stack[index++]); }
		double getDouble() { return stack[index++]; }
		vec2 getVec2() {
			double x = stack[index++];
			double y = stack[index++];
			return vec2(x, y);
		}
		vec3 getVec3() {
			double x = stack[index++];
			double y = stack[index++];
			double z = stack[index++];
			return vec3(x, y, z);
		}
		mat2 getMat2() {
			double a00 = stack[index++];
			double a01 = stack[index++];
			double a10 = stack[index++];
			double a11 = stack[index++];
			return mat2(a00, a01, a10, a11);
		}
		mat3 getMat3() {
			double a00 = stack[index++]; double a01 = stack[index++]; double a02 = stack[index++];
			double a10 = stack[index++]; double a11 = stack[index++]; double a12 = stack[index++];
			double a20 = stack[index++]; double a21 = stack[index++]; double a22 = stack[index++];
			return mat3(a00, a01, a02,
						a10, a11, a12,
						a20, a21, a22);
		}
	};

	using NativeFnc = std::function<Value(FuncArgs args)>;

	struct SourceLocation {
		int line = 0;
		int column = 0;
	};

	class FatalError : public std::runtime_error {
	public:
		FatalError(const std::string& message, const SourceLocation& loc)
			: std::runtime_error(message), location(loc) {}
		SourceLocation location;
	};
}
