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
		int localCount = 0;
		bool isNative = false;
		NativeFnc fnc;

		size_t maxStackSize = 0;
	};

	struct BinaryOperatorInfo
	{
		BinaryOp op;
		Type returnType = nullptr;
		Type leftType = nullptr;
		Type rightType = nullptr;
		size_t index = 0;
		BinaryFnc fnc;
	};

	class Program
	{
	public:
		struct Global
		{
			Type type = nullptr;
			int slot = -1;
			bool isInitialized = false;
			bool immutable = false;
			int refcount = 0;
		};

	public:
		Program();

		int addGlobal(const std::string& name, Type type, const Value& initializer = Value{}, bool immutable = false);
		int addGlobal(const std::string& name, Type type, std::initializer_list<Value> values, bool immutable = false);
		int addGlobal(const std::string& name, double v, bool immutable = false);


		Type RegisterStruct(const std::string& name, const std::vector<std::pair<Type, std::string>>& fields);
		Type RegisterStruct(const std::string& name, const std::vector<std::pair<TypeKind, std::string>>& fields);

		void registerNative(const std::string& name, Type returnType, std::vector<Type> argTypes, NativeFnc f);

		void registerNative(const std::string& name, double (*f)(double));

		void RegisterBinaryOperator(BinaryOp op, Type retType, Type type_l, Type type_r, BinaryFnc f);

		BinaryOperatorInfo* findBinaryOperatorOverload(BinaryOp op, Type type_l, Type type_r);

	public:
		std::unique_ptr<AST> ast;
		TypeRegistry types;

		std::vector<uint8_t> code;
		std::vector<Value> constants;
		std::vector<FunctionInfo> functions;
		std::vector<BinaryOperatorInfo> operators;

		std::unordered_map<std::string, Global> globals;
		std::unordered_map<int, Value> globalInitializers;
		std::unordered_map<std::string, NativeFnc> m_specialFns;

		size_t maxStackSize = 0;
	};
}
