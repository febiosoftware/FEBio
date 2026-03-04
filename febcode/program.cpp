#include "program.h"
using namespace febcode;

Program::Program()
{
	// Add "main" function (must be the first function, so it gets index 0)
	functions.push_back(FunctionInfo{
		"main",
		nullptr,
		std::vector<Type>{},
		0,
		0,
		false
		});
}

int Program::addGlobal(const std::string& name, Type type, const Value& initializer, bool immutable)
{
	if (globals.find(name) != globals.end())
		throw std::runtime_error("Global variable '" + name + "' is already defined.");

	int slot = (int)globals.size();
	globals[name] = { type, slot, !isVoid(initializer), immutable };
	globalInitializers[slot] = initializer;
	return slot;
}

int Program::addGlobal(const std::string& name, Type type, std::initializer_list<Value> values, bool immutable)
{
	if (globals.find(name) != globals.end())
		throw std::runtime_error("Global variable '" + name + "' is already defined.");

	int slot = (int)globals.size();
	globals[name] = { type, slot, true, immutable };

	if (type->kind == TypeKind::Struct)
	{
		if (values.size() != type->fields.size())
			throw std::runtime_error("Struct initializer has incorrect number of fields.");

		auto obj = createStruct(type);
		std::copy(values.begin(), values.end(), obj->fields.begin());

		globalInitializers[slot] = obj;
	}
	else if (type->kind == TypeKind::Array)
	{
		if (values.size() != type->arraySize)
			throw std::runtime_error("Array initializer has incorrect number of elements.");

		auto arr = createArray(type->elementType, type->arraySize);
		std::copy(values.begin(), values.end(), arr->elements.begin());
		globalInitializers[slot] = arr;
	}
	else
	{
		throw std::runtime_error("Initializer list can only be used for struct or array types.");
	}

	return slot;
}

int Program::addGlobal(const std::string& name, double v, bool immutable)
{
	return addGlobal(name, types.getDoubleType(), v, immutable);
}

Type Program::RegisterStruct(const std::string& name, const std::vector<std::pair<Type, std::string>>& fields)
{
	return types.defineStructType(name, fields);
}

Type Program::RegisterStruct(const std::string& name, const std::vector<std::pair<TypeKind, std::string>>& fields)
{
	return types.defineStructType(name, fields);
}


void Program::registerNative(const std::string& name, Type returnType, std::vector<Type> argTypes, NativeFnc fn)
{
	size_t slot = functions.size();
	functions.push_back(FunctionInfo{
		name,
		returnType,
		argTypes,
		slot,
		0,
		true,
		fn
		});
}

void Program::registerNative(const std::string& name, double (*f)(double))
{
	registerNative(
		name,
		types.getDoubleType(),
		{ types.getDoubleType() },
		[f](const std::vector<Value>& args) -> Value {
			double arg = std::get<double>(args[0]);
			return f(arg);
		});
}

void Program::RegisterBinaryOperator(BinaryOp op, Type retType, Type type_l, Type type_r, BinaryFnc f)
{
	size_t index = operators.size();
	operators.push_back(BinaryOperatorInfo{
		op,
		retType,
		type_l, type_r,
		index,
		f
		});
}

BinaryOperatorInfo* Program::findBinaryOperatorOverload(BinaryOp op, Type type_l, Type type_r)
{
	for (auto& overload : operators)
	{
		if (overload.op == op &&
			overload.leftType == type_l &&
			overload.rightType == type_r)
		{
			return &overload;
		}
	}
	return nullptr;
}
