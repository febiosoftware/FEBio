#include "program.h"
using namespace febcode;

Program::Program()
{
	ast = std::make_unique<AST>();

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

int Program::addGlobal(const std::string& name, Type type)
{
	// make sure the global variable name is unique
	auto it = globalIndices.find(name);
	if (it != globalIndices.end())
		throw std::runtime_error("Global variable '" + name + "' is already defined.");

	int slot = (int)globals.size();
	globals.push_back({ type, (int)globalStackSize, false, false });
	globalIndices[name] = slot;

	globalStackSize += type->size(); // reserve stack slots for this global variable

	return slot;
}

int Program::injectGlobal(const std::string& name, Type type)
{
	// make sure the global variable name is unique
	auto it = globalIndices.find(name);
	if (it != globalIndices.end())
		throw std::runtime_error("Global variable '" + name + "' is already defined.");

	int slot = (int)globals.size();
	globals.push_back({ type, (int)globalStackSize, true, true });
	globalIndices[name] = slot;

	globalStackSize += type->size(); // reserve stack slots for this global variable

	return slot;
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
	// calculate required stack size for arguments
	int argSize = 0;
	for (const auto& argType : argTypes)
	{
		argSize += (int)argType->size();
	}

	size_t slot = functions.size();
	functions.push_back(FunctionInfo{
		name,
		returnType,
		argTypes,
		slot,
		argSize,
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
		[f](FuncArgs args) -> Value {
			double arg = args.getDouble();
			return f(arg);
		});
}
