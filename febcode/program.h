#pragma once
#include "types.h"
#include "ast.h"

namespace febcode {

	struct FunctionInfo
	{
		std::string name;
		Type returnType = nullptr;
		std::vector<Type> args;
		size_t entry = 0;
		int argSize = 0; // stack size needed for arguments.
		bool isNative = false;
		NativeFnc fnc;

		size_t maxStackSize = 0;
	};

	class Program
	{
	public:
		struct Global
		{
			Type type = nullptr;
			int slot = -1;	// start index in stack
			bool immutable = false;
			int refcount = 0;
		};

		struct Input
		{
			Type type = nullptr;
			std::string name;
			int slot = -1;	// index in global's list
		};

	public:
		Program();

		// add a constant value to the program's constant pool and return its index
		uint8_t addConstant(const Value& value);

		// adds a global variable (from compiled code)
		int addGlobal(const std::string& name, Type type);

		// adds a global variable (from native code)
		int injectGlobal(const std::string& name, Type type);

		int addInput(const std::string& name, Type type);

		Type RegisterStruct(const std::string& name, const std::vector<std::pair<Type, std::string>>& fields);
		Type RegisterStruct(const std::string& name, const std::vector<std::pair<TypeKind, std::string>>& fields);

		void registerNative(const std::string& name, Type returnType, std::vector<Type> argTypes, NativeFnc f);

		void registerNative(const std::string& name, double (*f)(double));

		Type globalType(const std::string& name) const;

	public:
		int resolveFunction(const std::string& name, std::vector<Type> args);

		BinaryOpSignature resolveBinaryOp(BinaryOp op, Type left, Type right);

	public:
		std::unique_ptr<AST> ast;
		TypeRegistry types;

		Type returnType = nullptr; // expected return type of the program

		std::vector<uint8_t> code;
		std::vector<Value> constants;
		std::vector<FunctionInfo> functions;

		size_t globalStackSize = 0; // next available slot index for global variables
		std::vector<Global> globals;
		std::unordered_map<std::string, size_t> globalIndices; // maps global variable names to their slot index

		std::vector<Input> inputs;
		std::unordered_map<std::string, size_t> inputIndices; // maps global variable names to their slot index

		size_t maxStackSize = 0;

		// binary operator signatures for operator overloading
		std::unordered_map<BinaryOp, std::vector<BinaryOpSignature>> binaryOps;
	};
}
