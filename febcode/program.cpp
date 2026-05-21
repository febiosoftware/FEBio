#include "program.h"
#include "module_math.h"
#include "module_vec2.h"
#include "module_vec3.h"
#include "module_mat2.h"
#include "module_mat3.h"

using namespace febcode;

Program::Program()
{
	// create (empty) AST
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

	// register built-in binary operators for basic types
	binaryOps[BinaryOp::Plus] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Int   () },
		{ types.Double(), types.Double(), types.Double() },
		{ types.Double(), types.Int   (), types.Double() },
		{ types.Int   (), types.Double(), types.Double() },
	};
	binaryOps[BinaryOp::Minus] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Int   () },
		{ types.Double(), types.Double(), types.Double() },
		{ types.Double(), types.Int   (), types.Double() },
		{ types.Int   (), types.Double(), types.Double() },
	};
	binaryOps[BinaryOp::Multiply] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Int   () },
		{ types.Double(), types.Double(), types.Double() },
		{ types.Double(), types.Int   (), types.Double() },
		{ types.Int   (), types.Double(), types.Double() },
	};
	binaryOps[BinaryOp::Divide] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Int   () },
		{ types.Double(), types.Double(), types.Double() },
		{ types.Double(), types.Int   (), types.Double() },
		{ types.Int   (), types.Double(), types.Double() },
	};
	binaryOps[BinaryOp::Exponent] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Int   () },
		{ types.Double(), types.Double(), types.Double() },
		{ types.Double(), types.Int   (), types.Double() },
		{ types.Int   (), types.Double(), types.Double() },
	};
	binaryOps[BinaryOp::Greater] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Bool() },
		{ types.Double(), types.Double(), types.Bool() },
		{ types.Double(), types.Int   (), types.Bool() },
		{ types.Int   (), types.Double(), types.Bool() },
	};
	binaryOps[BinaryOp::GreaterEqual] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Bool() },
		{ types.Double(), types.Double(), types.Bool() },
		{ types.Double(), types.Int   (), types.Bool() },
		{ types.Int   (), types.Double(), types.Bool() },
	};
	binaryOps[BinaryOp::Less] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Bool() },
		{ types.Double(), types.Double(), types.Bool() },
		{ types.Double(), types.Int   (), types.Bool() },
		{ types.Int   (), types.Double(), types.Bool() },
	};
	binaryOps[BinaryOp::LessEqual] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Bool() },
		{ types.Double(), types.Double(), types.Bool() },
		{ types.Double(), types.Int   (), types.Bool() },
		{ types.Int   (), types.Double(), types.Bool() },
	};
	binaryOps[BinaryOp::EqualEqual] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Bool() },
		{ types.Double(), types.Double(), types.Bool() },
		{ types.Double(), types.Int   (), types.Bool() },
		{ types.Int   (), types.Double(), types.Bool() },
	};
	binaryOps[BinaryOp::NotEqual] = {
		//   LHS           RHS             Result
		{ types.Int   (), types.Int   (), types.Bool() },
		{ types.Double(), types.Double(), types.Bool() },
		{ types.Double(), types.Int   (), types.Bool() },
		{ types.Int   (), types.Double(), types.Bool() },
	};
	binaryOps[BinaryOp::AndAnd] = {
		//   LHS           RHS             Result
		{ types.Bool(), types.Bool(), types.Bool() },
	};
	binaryOps[BinaryOp::OrOr] = {
		//   LHS           RHS             Result
		{ types.Bool(), types.Bool(), types.Bool() },
	};

	// register modules
	MathModule mathModule;
	mathModule.Register(*this);

	Vec2Module vec2Module;
	vec2Module.Register(*this);

	Vec3Module vec3Module;
	vec3Module.Register(*this);

	Mat2Module mat2Module;
	mat2Module.Register(*this);

	Mat3Module mat3Module;
	mat3Module.Register(*this);
}

uint8_t Program::addConstant(const Value& value)
{
	// see if this value already exists in the constant pool
	for (size_t i = 0; i < constants.size(); ++i)
	{
		if (constants[i] == value)
			return (uint8_t)i;
	}

	// add new constant to the pool
	constants.push_back(value);
	return (uint8_t)(constants.size() - 1);
}

int Program::addGlobal(const std::string& name, Type type)
{
	// make sure the global variable name is unique
	auto it = globalIndices.find(name);
	if (it != globalIndices.end())
		throw std::runtime_error("Global variable '" + name + "' is already defined.");

	int slot = (int)globals.size();
	globals.push_back({ type, (int)globalStackSize, false, 0 });
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
	globals.push_back({ type, (int)globalStackSize, true, 0 });
	globalIndices[name] = slot;

	globalStackSize += type->size(); // reserve stack slots for this global variable

	return slot;
}

int Program::addInput(const std::string& name, Type type)
{
	int slot = injectGlobal(name, type);
	inputIndices[name] = inputs.size();
	inputs.push_back({ type, name, slot });
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
		types.Double(),
		{ types.Double() },
		[f](FuncArgs args) -> Value {
			double arg = args.getDouble();
			return f(arg);
		});
}

Type Program::globalType(const std::string& name) const
{
	auto it = globalIndices.find(name);
	if (it == globalIndices.end())
		throw std::runtime_error("Global variable '" + name + "' is not defined.");
	size_t slot = it->second;
	return globals[slot].type;
}

int Program::resolveFunction(const std::string& name, std::vector<Type> args)
{
	for (int i = 0; i < (int)functions.size(); ++i)
	{
		if ((functions[i].name == name) && (functions[i].args == args))
		{
			return i;
		}
	}
	return -1;
}

BinaryOpSignature Program::resolveBinaryOp(BinaryOp op, Type left, Type right)
{
	auto it = binaryOps.find(op);
	if (it != binaryOps.end())
	{
		for (const auto& sig : it->second)
		{
			if ((sig.leftType == left) && (sig.rightType == right))
			{
				return sig;
			}
			if (sig.leftType == left)
			{
				// try implicit conversion of right operand
				Type convertedRight = coerce(right, sig.rightType);
				if (convertedRight != nullptr)
					return sig;
			}
			if (sig.rightType == right)
			{
				// try implicit conversion of left operand
				Type convertedLeft = coerce(left, sig.leftType);
				if (convertedLeft != nullptr)
					return sig;
			}
		}
	}

	return BinaryOpSignature({ nullptr, nullptr, nullptr });
}
