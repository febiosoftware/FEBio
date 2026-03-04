#pragma once
#include <memory>
#include <variant>
#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <functional>

namespace febcode
{
	enum class TypeKind : uint8_t {
		Void,
		Bool,
		Int,
		Double,
		String,
		Array,
		Struct
	};

	struct TypeStruct;
	using Type = const TypeStruct*;

	struct TypeStruct {
		TypeKind kind = TypeKind::Void;

		// for arrays
		Type elementType = nullptr;
		size_t arraySize = 0;

		// For struct types
		int typeIndex = -1;
		std::string name;
		std::vector<std::pair<Type, std::string>> fields;
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
		STRING,
		ARRAY,
		STRUCT
	};

	// variant for storing either a number or a boolean
	using Value = std::variant <
		std::monostate,
		bool,
		int,
		double,
		std::string,
		ArrayValuePtr,
		StructValuePtr
	>;

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

	inline bool isVoid  (const Value& v) { return v.index() == ValueIndex::VOID;}
	inline bool isBool  (const Value& v) { return v.index() == ValueIndex::BOOL;}
	inline bool isInt   (const Value& v) { return v.index() == ValueIndex::INT; }
	inline bool isDouble(const Value& v) { return v.index() == ValueIndex::DOUBLE; }
	inline bool isString(const Value& v) { return v.index() == ValueIndex::STRING; }
	inline bool isArray (const Value& v) { return v.index() == ValueIndex::ARRAY; }
	inline bool isStruct(const Value& v) { return v.index() == ValueIndex::STRUCT; }

	inline bool isVoidType  (Type type) { return type->kind == TypeKind::Void; }
	inline bool isBoolType  (Type type) { return type->kind == TypeKind::Bool; }
	inline bool isIntType   (Type type) { return type->kind == TypeKind::Int; }
	inline bool isDoubleType(Type type) { return type->kind == TypeKind::Double; }
	inline bool isStringType(Type type) { return type->kind == TypeKind::String; }
	inline bool isStructType(Type type) { return type->kind == TypeKind::Array; }
	inline bool isArrayType (Type type) { return type->kind == TypeKind::Struct; }

	inline bool isNumericType(Type type) { return isIntType(type) || isDoubleType(type); }

	inline bool   getBool  (const Value& v) { return std::get<bool  >(v); }
	inline int    getInt   (const Value& v) { return std::get<int   >(v); }
	inline double getDouble(const Value& v) { return std::get<double>(v); }
	inline const std::string& getString(const Value& v) { return std::get<std::string>(v); }
	inline const ArrayValue&  getArray (const Value& v) { return *std::get<std::shared_ptr<ArrayValue>>(v); }
	inline const StructValue& getStruct(const Value& v) { return *std::get<std::shared_ptr<StructValue>>(v); }

	inline const ArrayValuePtr& getArrayPtr(const Value& v) { return std::get<std::shared_ptr<ArrayValue>>(v); }
	inline const StructValuePtr& getStructPtr(const Value& v) { return std::get<std::shared_ptr<StructValue>>(v); }

	inline std::string& getString(Value& v) { return std::get<std::string>(v); }
	inline ArrayValue&  getArray (Value& v) { return *std::get<std::shared_ptr<ArrayValue>>(v); }
	inline StructValue& getStruct(Value& v) { return *std::get<std::shared_ptr<StructValue>>(v); }

	using NativeFnc = std::function<Value(const std::vector<Value>&)>;

	using BinaryFnc = std::function<Value(const Value&, const Value&)>;

	class TypeRegistry {
	public:
		TypeRegistry();
		void clear();

		Type getVoidType() const;
		Type getIntType() const;
		Type getBoolType() const;
		Type getDoubleType() const;
		Type getStringType() const;
		Type getArrayType(Type element, size_t size);
		Type getArrayType(Type element, size_t size) const;

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
		TypeStruct m_string;

		// Interned array types
		std::vector<std::unique_ptr<TypeStruct>> m_arrayTypes;

		// User-defined struct types
		std::vector<std::unique_ptr<TypeStruct>> m_structTypes;
	};

	inline TypeKind ValueType(const Value& v)
	{
		return (TypeKind)v.index();
	}

	std::string TypeToString(Type type);

}

inline febcode::Value operator + (const febcode::Value& a, const febcode::Value& b)
{
	if (isInt(a) && isInt(b))
		return std::get<int>(a) + std::get<int>(b);
	if (isDouble(a) && isDouble(b))
		return std::get<double>(a) + std::get<double>(b);
	if (isString(a) && isString(b))
		return std::get<std::string>(a) + std::get<std::string>(b);
	throw std::runtime_error("Unsupported operand types for +");
}
