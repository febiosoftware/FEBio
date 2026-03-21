#pragma once
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>
#include <functional>
#include <assert.h>

namespace febcode
{
	enum class TypeKind : uint8_t {
		Void,
		Bool,
		Int,
		Double,
		Vec2,
		Vec3,
		Array,
		Struct,
	};

	enum class RefType : uint8_t {
		Invalid,
		Double,
		Value
	};

	struct Ref {
		void* ptr;
		RefType type;
		Ref() : ptr(nullptr), type(RefType::Invalid) {}
	};

	struct vec2
	{
		double x, y;

		vec2() : x(0), y(0) {}
		vec2(double x, double y) : x(x), y(y) {}

		vec2 operator+(const vec2& other) const { return vec2(x + other.x, y + other.y); }
		vec2 operator-(const vec2& other) const { return vec2(x - other.x, y - other.y); }
		vec2 operator*(double scalar) const { return vec2(x * scalar, y * scalar); }

		// dot product
		double operator*(const vec2& other) const { return x * other.x + y * other.y; }

		bool operator==(const vec2& other) const { return (x == other.x) && (y == other.y); }
		bool operator!=(const vec2& other) const { return !(*this == other); }
	};

	struct vec3
	{
		double x, y, z;

		vec3() : x(0), y(0), z(0) {}
		vec3(double x, double y, double z) : x(x), y(y), z(z) {}

		vec3 operator+(const vec3& other) const { return vec3(x + other.x, y + other.y, z + other.z); }
		vec3 operator-(const vec3& other) const { return vec3(x - other.x, y - other.y, z - other.z); }
		vec3 operator*(double scalar) const { return vec3(x * scalar, y * scalar, z * scalar); }

		// dot product
		double operator*(const vec3& other) const { return x * other.x + y * other.y + z * other.z; }

		bool operator==(const vec3& other) const { return (x == other.x) && (y == other.y) && (z == other.z); }
		bool operator!=(const vec3& other) const { return !(*this == other); }

		vec3 cross(const vec3& other) const {
			return vec3(
				y * other.z - z * other.y,
				z * other.x - x * other.z,
				x * other.y - y * other.x
			);
		}
	};

	struct TypeStruct;
	using Type = const TypeStruct*;

	using StructField = std::pair<Type, std::string>; // type and name of a struct field

	struct TypeStruct {
		TypeKind kind = TypeKind::Void;

		// for arrays
		Type elementType = nullptr;
		size_t arraySize = 0;

		// For struct types
		int typeIndex = -1;
		std::string name;
		std::vector<StructField> fields;
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
			ArrayValuePtr arrayValue;
			StructValuePtr structValue;
			Ref ref;
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
			case ValueIndex::ARRAY: new (&arrayValue) ArrayValuePtr(v.arrayValue); break;
			case ValueIndex::STRUCT: new (&structValue) StructValuePtr(v.structValue); break;
			case ValueIndex::REF   : ref = v.ref; break;
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
	inline bool isArray (const Value& v) { return v.index == ValueIndex::ARRAY; }
	inline bool isStruct(const Value& v) { return v.index == ValueIndex::STRUCT; }
	inline bool isRef   (const Value& v) { return v.index == ValueIndex::REF; }

	inline bool isVoidType  (Type type) { return type->kind == TypeKind::Void; }
	inline bool isBoolType  (Type type) { return type->kind == TypeKind::Bool; }
	inline bool isIntType   (Type type) { return type->kind == TypeKind::Int; }
	inline bool isDoubleType(Type type) { return type->kind == TypeKind::Double; }
	inline bool isVec2Type  (Type type) { return type->kind == TypeKind::Vec2; }
	inline bool isVec3Type  (Type type) { return type->kind == TypeKind::Vec3; }
	inline bool isArrayType (Type type) { return type->kind == TypeKind::Array; }
	inline bool isStructType(Type type) { return type->kind == TypeKind::Struct; }

	inline bool isNumericType(Type type) { return isIntType(type) || isDoubleType(type); }

	inline bool   getBool(const Value& v) { assert(v.index == ValueIndex::BOOL); return v.b; }
	inline int    getInt   (const Value& v) { assert(v.index == ValueIndex::INT); return v.i; }
	inline double getDouble(const Value& v) { assert(v.index == ValueIndex::DOUBLE); return v.d; }
	inline const vec2& getVec2(const Value& v) { assert(v.index == ValueIndex::VEC2); return v.vec2Value; }
	inline vec2& getVec2(Value& v) { assert(v.index == ValueIndex::VEC2); return v.vec2Value; }
	inline const vec3& getVec3(const Value& v) { assert(v.index == ValueIndex::VEC3); return v.vec3Value; }
	inline vec3& getVec3(Value& v) { assert(v.index == ValueIndex::VEC3); return v.vec3Value; }
	inline const ArrayValue&  getArray (const Value& v) { assert(v.index == ValueIndex::ARRAY); return *v.arrayValue; }
	inline const StructValue& getStruct(const Value& v) { assert(v.index == ValueIndex::STRUCT); return *v.structValue; }
	inline const Ref& getRef(const Value& v) { assert(v.index == ValueIndex::REF); return v.ref; }

	inline const ArrayValuePtr& getArrayPtr(const Value& v) { assert(v.index == ValueIndex::ARRAY); return v.arrayValue; }
	inline const StructValuePtr& getStructPtr(const Value& v) { assert(v.index == ValueIndex::STRUCT); return v.structValue; }

	inline ArrayValue&  getArray (Value& v) { assert(v.index == ValueIndex::ARRAY); return *v.arrayValue; }
	inline StructValue& getStruct(Value& v) { assert(v.index == ValueIndex::STRUCT); return *v.structValue; }

	inline bool isZero(const Value& v)
	{
		if (isInt(v) && (getInt(v) == 0)) return true;
		if (isDouble(v) && (getDouble(v) == 0.0)) return true;
		return false;
	}

	inline bool isOne(const Value& v)
	{
		if (isInt(v) && (getInt(v) == 1)) return true;
		if (isDouble(v) && (getDouble(v) == 1.0)) return true;
		return false;
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

	using NativeFnc = std::function<Value(const Value* arg, int argc)>;

	using BinaryFnc = std::function<Value(const Value&, const Value&)>;

	class TypeRegistry {
	public:
		TypeRegistry();
		void clear();

		Type getVoidType() const;
		Type getIntType() const;
		Type getBoolType() const;
		Type getDoubleType() const;
		Type getVec2Type() const;
		Type getVec3Type() const;
		Type getArrayType(Type element, size_t size);
		Type getArrayType(Type element, size_t size) const;
		Type getArrayType(Type element, const std::vector<size_t>& size) const;

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

		// Interned array types
		std::vector<std::unique_ptr<TypeStruct>> m_arrayTypes;

		// User-defined struct types
		std::vector<std::unique_ptr<TypeStruct>> m_structTypes;
	};

	inline TypeKind ValueType(const Value& v)
	{
		return (TypeKind)v.index;
	}

	std::string TypeToString(Type type);

}

inline febcode::Value operator + (const febcode::Value& a, const febcode::Value& b)
{
	if (isInt(a) && isInt(b))
		return getInt(a) + getInt(b);
	if (isDouble(a) && isDouble(b))
		return getDouble(a) + getDouble(b);
	throw std::runtime_error("Unsupported operand types for +");
}

inline febcode::Value operator - (const febcode::Value& a, const febcode::Value& b)
{
	if (isInt(a) && isInt(b))
		return getInt(a) - getInt(b);
	if (isDouble(a) && isDouble(b))
		return getDouble(a) - getDouble(b);
	throw std::runtime_error("Unsupported operand types for +");
}

inline febcode::Value operator * (const febcode::Value& a, const febcode::Value& b)
{
	if (isInt(a) && isInt(b))
		return getInt(a) * getInt(b);
	if (isDouble(a) && isDouble(b))
		return getDouble(a) * getDouble(b);
	throw std::runtime_error("Unsupported operand types for +");
}

inline febcode::Value operator / (const febcode::Value& a, const febcode::Value& b)
{
	if (isDouble(a) && isDouble(b))
		return getDouble(a) / getDouble(b);
	throw std::runtime_error("Unsupported operand types for +");
}
