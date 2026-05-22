#include "compiler.h"
#include "scanner.h"
#include "parser.h"
#include "resolver.h"

using namespace febcode;

void febcode::CompileSource(Program& prg, const std::string& source)
{
	// 1. Tokenize
	Scanner scanner(source);
	std::vector<Token> tokens = scanner.scanTokens();

	// 2. Parse -> AST
	Parser parser(prg);
	parser.parse(tokens);

	// Optional safety check
	if (prg.ast->empty())
		throw std::runtime_error("Parser produced no statements.");

	// 3. semantic analysis
	Resolver resolve(prg);
	resolve.resolve();

	// 4. Compile -> bytecode
	Compiler compiler(prg);
	compiler.compile();
}

//
// ================= IMPLEMENTATION =================
//

Compiler::Compiler(Program& prg) : prg(prg)
{
}

//
// ===== Scope =====
//

void Compiler::beginScope()
{
	m_scopeDepth++;
}

void Compiler::endScope()
{
	while (!m_locals.empty() && m_locals.back().depth == m_scopeDepth)
	{
		Local& local = m_locals.back();

		switch (local.type->kind)
		{
		case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
		case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
		case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
		case TypeKind::Vec2  : emit(OpCode::POP_VEC2  ); break;
		case TypeKind::Vec3  : emit(OpCode::POP_VEC3  ); break;
		case TypeKind::Mat2  : emit(OpCode::POP_MAT2  ); break;
		case TypeKind::Mat3  : emit(OpCode::POP_MAT3  ); break;
		case TypeKind::Array : 
			if (local.type->size() > 255)
				throw std::runtime_error("Array size exceeds maximum of 255.");
			emit(OpCode::POP_ARRAY, (int)local.type->size());
			emitUint8((uint8_t)local.type->size());
			break;
		case TypeKind::Struct:
			emit(OpCode::POP_STRUCT, (int)local.type->size());
			emitUint8((uint8_t)local.type->typeIndex);
		break;		
		default:
			throw std::runtime_error("Unknown type kind in endScope.");
			break;
		}

		m_locals.pop_back();
		localStackSize -= (int)local.type->size();
	}
	m_scopeDepth--;
}

int Compiler::resolveLocal(const std::string& name)
{
	for (int i = (int)m_locals.size() - 1; i >= 0; --i)
		if (m_locals[i].name == name)
			return i;
	return -1;
}

int Compiler::resolveGlobal(const std::string& name)
{
	auto it = prg.globalIndices.find(name);
	if (it != prg.globalIndices.end())
		return (int)it->second;

	return -1;
}

Type Compiler::resolveVariableType(const std::string& name)
{
	int localIndex = resolveLocal(name);
	if (localIndex != -1)
		return m_locals[localIndex].type;
	int globalIndex = resolveGlobal(name);
	if (globalIndex != -1)
		return prg.globals[globalIndex].type;
	throw std::runtime_error("Undefined variable: " + name);
}

//
// ===== Bytecode Helpers =====
//

void Compiler::emit(OpCode op, int arg)
{
	prg.code.push_back((uint8_t)op);

	stackDepth += stackEffect(op, arg); assert(stackDepth >= 0);
	if (stackDepth > maxStackDepth)
		maxStackDepth = stackDepth;
}

void Compiler::emitUint8(uint8_t v)
{
	prg.code.push_back(v);
}

void Compiler::emitUint16(uint16_t v)
{
	prg.code.push_back((v >> 8) & 0xff);
	prg.code.push_back(v & 0xff);
}

int Compiler::stackEffect(OpCode op, int arg)
{
	switch (op)
	{
	case OpCode::PUSH_VOID: 
	case OpCode::PUSH_BOOL: 
	case OpCode::PUSH_INT: 
	case OpCode::PUSH_DOUBLE: 
		return +1;

	case OpCode::PUSH_LOCAL: return +arg;

	case OpCode::PUSH_VEC2: return +2;
	case OpCode::PUSH_VEC3: return +3;
	case OpCode::PUSH_MAT2: return +4;
	case OpCode::PUSH_MAT3: return +9;

	case OpCode::GET_GLOBAL_BOOL  : return +1;
	case OpCode::GET_GLOBAL_INT   : return +1;
	case OpCode::GET_GLOBAL_DOUBLE: return +1;
	case OpCode::GET_GLOBAL_VEC2  : return +2;
	case OpCode::GET_GLOBAL_VEC3  : return +3;
	case OpCode::GET_GLOBAL_MAT2  : return +4;
	case OpCode::GET_GLOBAL_MAT3  : return +9;
	case OpCode::GET_GLOBAL_ARRAY : return +arg;
	case OpCode::GET_GLOBAL_STRUCT: return +arg;

	case OpCode::SET_GLOBAL_BOOL: 
	case OpCode::SET_GLOBAL_INT: 
	case OpCode::SET_GLOBAL_DOUBLE: 
	case OpCode::SET_GLOBAL_VEC2: 
	case OpCode::SET_GLOBAL_VEC3: 
	case OpCode::SET_GLOBAL_MAT2: 
	case OpCode::SET_GLOBAL_MAT3: 
	case OpCode::SET_GLOBAL_ARRAY: 
	case OpCode::SET_GLOBAL_STRUCT: 
		return 0;

	case OpCode::GET_GLOBAL_REF       : return +1;

	case OpCode::GET_LOCAL_BOOL  : return +1;
	case OpCode::GET_LOCAL_INT   : return +1;
	case OpCode::GET_LOCAL_DOUBLE: return +1;
	case OpCode::GET_LOCAL_VEC2  : return +2;
	case OpCode::GET_LOCAL_VEC3  : return +3;
	case OpCode::GET_LOCAL_MAT2  : return +4;
	case OpCode::GET_LOCAL_MAT3  : return +9;
	case OpCode::GET_LOCAL_ARRAY : return +arg;
	case OpCode::GET_LOCAL_STRUCT: return +arg;

	case OpCode::GET_LOCAL_REF       : return +1;

	case OpCode::GET_PROPERTY_BOOL  : return +1;
	case OpCode::GET_PROPERTY_INT   : return +1;
	case OpCode::GET_PROPERTY_DOUBLE: return +1;
	case OpCode::GET_PROPERTY_VEC2  : return +2;
	case OpCode::GET_PROPERTY_VEC3  : return +3;
	case OpCode::GET_PROPERTY_MAT2  : return +4;
	case OpCode::GET_PROPERTY_MAT3  : return +9;
	case OpCode::GET_PROPERTY_ARRAY : return +1;
	case OpCode::GET_PROPERTY_STRUCT: return +1;

	case OpCode::GET_MEMBER_REF: return 0;

	case OpCode::GET_INDEX_BOOL  : return +arg;
	case OpCode::GET_INDEX_INT   : return +arg;
	case OpCode::GET_INDEX_DOUBLE: return +arg;
	case OpCode::GET_INDEX_VEC2  : return +arg;
	case OpCode::GET_INDEX_VEC3  : return +arg;
	case OpCode::GET_INDEX_MAT2  : return +arg;
	case OpCode::GET_INDEX_MAT3  : return +arg;
	case OpCode::GET_INDEX_ARRAY : return +arg;
	case OpCode::GET_INDEX_STRUCT: return +arg;

	case OpCode::GET_GLOBAL_INDEX_DOUBLE: return +1;

	case OpCode::GET_INDEX_REF:
	case OpCode::GET_INDEX_REF_BOOL:
	case OpCode::GET_INDEX_REF_INT:
	case OpCode::GET_INDEX_REF_DOUBLE:
	case OpCode::GET_INDEX_REF_VEC2:
	case OpCode::GET_INDEX_REF_VEC3:
	case OpCode::GET_INDEX_REF_MAT2:
	case OpCode::GET_INDEX_REF_MAT3:
		return 0;

	case OpCode::GET_VEC2_X: return -1;
	case OpCode::GET_VEC2_Y: return -1;
	case OpCode::GET_VEC2_X_REF: return 0;
	case OpCode::GET_VEC2_Y_REF: return 0;
	case OpCode::GET_VEC2_SWIZZLE: return +1; // in case the swizzle returns a vec3
	case OpCode::GET_VEC2_INDEX: return -2; // pops the index and vec2, and pushes the component

	case OpCode::GET_VEC3_X: return 0;
	case OpCode::GET_VEC3_Y: return 0;
	case OpCode::GET_VEC3_Z: return 0;
	case OpCode::GET_VEC3_X_REF: return 0;
	case OpCode::GET_VEC3_Y_REF: return 0;
	case OpCode::GET_VEC3_Z_REF: return 0;
	case OpCode::GET_VEC3_SWIZZLE: return +1;
	case OpCode::GET_VEC3_INDEX: return -3; // pops the index and vec3, and pushes the component

	case OpCode::NEG_INT: return 0;
	case OpCode::ADD_INT: return -1;
	case OpCode::SUB_INT: return -1;
	case OpCode::MUL_INT: return -1;
	case OpCode::DIV_INT: return -1;
	case OpCode::EXP_INT: return -1;
	case OpCode::GT_INT: return -1;
	case OpCode::LT_INT: return -1;
	case OpCode::GE_INT: return -1;
	case OpCode::LE_INT: return -1;
	case OpCode::NEG_DOUBLE: return 0;
	case OpCode::SQR_DOUBLE: return 0;
	case OpCode::SQRT_DOUBLE: return 0;
	case OpCode::ADD_DOUBLE: return -1;
	case OpCode::SUB_DOUBLE: return -1;
	case OpCode::MUL_DOUBLE: return -1;
	case OpCode::DIV_DOUBLE: return -1;
	case OpCode::EXP_DOUBLE: return -1;
	case OpCode::GT_DOUBLE: return -1;
	case OpCode::LT_DOUBLE: return -1;
	case OpCode::GE_DOUBLE: return -1;
	case OpCode::LE_DOUBLE: return -1;

	case OpCode::CREATE_VEC2_1ARG: return +1; // pops the double, pushes the vec2
	case OpCode::NEG_VEC2: return 0;
	case OpCode::ADD_VEC2: return -2;
	case OpCode::SUB_VEC2: return -2;
	case OpCode::DOT_VEC2: return -3;
	case OpCode::MUL_VEC2_DOUBLE: return -1;
	case OpCode::MUL_DOUBLE_VEC2: return -1;
	case OpCode::DIV_VEC2_DOUBLE: return -1;

	case OpCode::CREATE_VEC3_1ARG: return +2; // pops the double, pushes the vec3
	case OpCode::NEG_VEC3: return 0;
	case OpCode::ADD_VEC3: return -3;
	case OpCode::SUB_VEC3: return -3;
	case OpCode::DOT_VEC3: return -5;
	case OpCode::MUL_VEC3_DOUBLE: return -1;
	case OpCode::MUL_DOUBLE_VEC3: return -1;
	case OpCode::DIV_VEC3_DOUBLE: return -1;

	case OpCode::NEG_MAT2: return 0;
	case OpCode::ADD_MAT2: return -4;
	case OpCode::SUB_MAT2: return -4;
	case OpCode::MUL_MAT2: return -4;
	case OpCode::MUL_MAT2_DOUBLE: return -1;
	case OpCode::MUL_DOUBLE_MAT2: return -1;
	case OpCode::DIV_MAT2_DOUBLE: return -1;
	case OpCode::MUL_MAT2_VEC2 : return -4;
	case OpCode::GET_MAT2_INDEX: return -3;
	case OpCode::CREATE_MAT2_DIAG: return +3;

	case OpCode::NEG_MAT3: return 0;
	case OpCode::ADD_MAT3: return -9;
	case OpCode::SUB_MAT3: return -9;
	case OpCode::MUL_MAT3: return -9;
	case OpCode::MUL_MAT3_DOUBLE: return -1;
	case OpCode::DIV_MAT3_DOUBLE: return -1;
	case OpCode::MUL_DOUBLE_MAT3: return -1;
	case OpCode::MUL_MAT3_VEC3  : return -9;
	case OpCode::GET_MAT3_INDEX : return -7;
	case OpCode::ADD_GLOBAL_MAT3: return +9;
	case OpCode::SUB_GLOBAL_MAT3: return +9;
	case OpCode::MUL_GLOBAL_MAT3: return +9;
	case OpCode::CREATE_MAT3_DIAG: return +8;
	case OpCode::CREATE_MAT3_VEC3: return +0;

	case OpCode::NOT: return 0;

	case OpCode::EQUAL_BOOL: 
	case OpCode::EQUAL_INT: 
	case OpCode::EQUAL_DOUBLE: 
		return -1;

	case OpCode::NEQ_BOOL: 
	case OpCode::NEQ_INT: 
	case OpCode::NEQ_DOUBLE: 
		return -1;

	case OpCode::JUMP: return 0;
	case OpCode::JUMP_IF_FALSE: return +1; // technically +0, but each jump adds two pops (one for each branch)
	case OpCode::JUMP_IF_TRUE: return +1; // technically +0, but each jump adds two pops (one for each branch)
	case OpCode::LOOP: return 0;

	case OpCode::STORE_BOOL:
	case OpCode::STORE_INT:
	case OpCode::STORE_DOUBLE:
	case OpCode::STORE_VEC2:
	case OpCode::STORE_VEC3:
	case OpCode::STORE_MAT2:
	case OpCode::STORE_MAT3:
	case OpCode::STORE_ARRAY:
	case OpCode::STORE_STRUCT:
		return -1;

	case OpCode::POP_VOID  : return -1;
	case OpCode::POP_BOOL  : return -1;
	case OpCode::POP_INT   : return -1;
	case OpCode::POP_DOUBLE: return -1;
	case OpCode::POP_VEC2  : return -2;
	case OpCode::POP_VEC3  : return -3;
	case OpCode::POP_MAT2  : return -4;
	case OpCode::POP_MAT3  : return -9;
	case OpCode::POP_ARRAY : return -arg;
	case OpCode::POP_STRUCT: return -arg;

	case OpCode::CALL: return arg; // this is the net stack effect calculated when compiling.

	case OpCode::RETURN_VOID:
	case OpCode::RETURN_BOOL: 
	case OpCode::RETURN_INT: 
	case OpCode::RETURN_DOUBLE: 
	case OpCode::RETURN_VEC2: 
	case OpCode::RETURN_VEC3: 
	case OpCode::RETURN_MAT2: 
	case OpCode::RETURN_MAT3: 
	case OpCode::RETURN_ARRAY: 
	case OpCode::RETURN_STRUCT: 
		return 0;
	default:
		assert(false);
		return 0;
	}
}

int Compiler::emitJump(OpCode op)
{
	emit(op);
	emitUint16(0xffff);
	return (int)prg.code.size() - 2;
}

void Compiler::patchJump(int offset)
{
	int jump = (int)prg.code.size() - offset - 2;
	prg.code[offset] = (jump >> 8) & 0xff;
	prg.code[offset + 1] = jump & 0xff;
}

void Compiler::emitLoop(int loopStart)
{
	emit(OpCode::LOOP);
	int offset = (int)prg.code.size() - loopStart + 2;
	emitUint16(offset);
}

//
// ===== Program =====
//

void Compiler::compile()
{
	prg.functions[0].entry = 0;
	hasReturn = false;

	expectedReturnType = prg.returnType;

	for (auto& stmt : prg.ast->root.statements)
		compileStatement(stmt.get());

	// only add return if no return was encountered.
	if (!hasReturn)
		emit(OpCode::RETURN_VOID);

	prg.maxStackSize = maxStackDepth;
}

//
// ===== Statements =====
//

void Compiler::compileStatement(Statement* stmt)
{
	if      (auto b = dynamic_cast<ExpressionStmt*>(stmt)) compileExprStmt(b);
	else if (auto b = dynamic_cast<BlockStmt*     >(stmt)) compileBlock(b);
	else if (auto v = dynamic_cast<VarDeclStmt*   >(stmt)) compileVarDecl(v);
	else if (auto f = dynamic_cast<FunctionStmt*  >(stmt)) compileFunction(f);
	else if (auto i = dynamic_cast<IfStmt*        >(stmt)) compileIf(i);
	else if (auto w = dynamic_cast<WhileStmt*     >(stmt)) compileWhile(w);
	else if (auto l = dynamic_cast<ForStmt*       >(stmt)) compileFor(l);
	else if (auto r = dynamic_cast<ReturnStmt*    >(stmt)) compileReturn(r);
	else if (auto s = dynamic_cast<StructStmt*    >(stmt)) compileStruct(s);
	else
		error(stmt, "Unsupported statement type");
}

void Compiler::compileExprStmt(ExpressionStmt* stmt)
{
	Type type = compileExpression(stmt->expr.get());
	pop(type);
}

void Compiler::pop(Type type)
{
	switch (type->kind)
	{
	case TypeKind::Void  : emit(OpCode::POP_VOID  ); break;
	case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
	case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
	case TypeKind::Vec2  : emit(OpCode::POP_VEC2  ); break;
	case TypeKind::Vec3  : emit(OpCode::POP_VEC3  ); break;
	case TypeKind::Mat2  : emit(OpCode::POP_MAT2  ); break;
	case TypeKind::Mat3  : emit(OpCode::POP_MAT3  ); break;
	case TypeKind::Array : 
		if (type->size() > 255)
			throw std::runtime_error("Array size exceeds maximum of 255.");
		emit(OpCode::POP_ARRAY, (int)type->size());
		emitUint8((uint8_t)type->size());
		break;
	case TypeKind::Struct: 
		emit(OpCode::POP_STRUCT, (int)type->size());
		emitUint8((uint8_t)type->typeIndex);
		break;
	default:
		throw std::runtime_error("Unsupported expression type in expression statement");
	}
}

void Compiler::compileBlock(BlockStmt* stmt)
{
	beginScope();
	for (auto& s : stmt->statements)
		compileStatement(s.get());
	endScope();
}

Type Compiler::compileInitializer(InitExpr* init)
{
	if (init->elements.empty())
		error(init, "Initializer cannot be empty.");

	// deduce the array type from the first element
	Type elemType = compileExpression(init->elements[0].get());
	Type arrayType = prg.types.getArrayType(elemType, { init->elements.size() });

	// compile the rest of the elements and check that they match the deduced type
	for (size_t i = 1; i < init->elements.size(); ++i)
	{
		Type type = coerce(compileExpression(init->elements[i].get()), elemType);
		if (type != elemType)
			error(init->elements[i].get(), "Initializer element type mismatch.");
	}

	return arrayType;
}

Type Compiler::compileConstructor(ConstructorExpr* construct)
{
	// check the number of arguments for the constructor
	int nargs = (int)construct->args.size();
	switch (construct->valType->kind)
	{
	case TypeKind::Vec2:
		if ((nargs != 1) && (nargs != 2))
			throw std::runtime_error("Vec2 constructor must have 1 or 2 arguments.");
		break;
	case TypeKind::Vec3:
		if ((nargs != 1) && (nargs != 3))
			throw std::runtime_error("Vec3 constructor must have 1 or 3 arguments.");
		break;
	case TypeKind::Mat2:
		if (nargs != 1 && nargs != 4)
			throw std::runtime_error("Mat2 constructor must have either 1 or 4 arguments.");
		break;
	case TypeKind::Mat3:
		if (nargs != 1 && nargs != 3 && nargs != 9)
			throw std::runtime_error("Mat3 constructor must have either 1, 3, or 9 arguments.");
		break;
	case TypeKind::Struct:
		if (nargs != construct->valType->fields.size())
			throw std::runtime_error("Struct constructor must have exactly " + std::to_string(construct->valType->fields.size()) + " arguments.");
		break;
	default:
		break;
	}

	// process each argument and check their types
	std::vector<Type> argTypes;
	for (int i = 0; i < construct->args.size(); ++i)
	{
		Type argType = compileExpression(construct->args[i].get());
		argTypes.push_back(argType);

		// check the argument type against the expected type for the constructor
		switch (construct->valType->kind)
		{
		case TypeKind::Vec2:
		case TypeKind::Vec3:
		case TypeKind::Mat2:
			argType = coerce(argType, prg.types.Double());
			if (argType != prg.types.Double())
				throw std::runtime_error("Vec2 constructor arguments must be of type double.");
			break;
		case TypeKind::Mat3:
			if (!isScalarType(argType) && argType != prg.types.Vec3())
				throw std::runtime_error("Mat3 constructor arguments must be of type int, double or vec3.");
			break;
		case TypeKind::Struct:
			argType = coerce(argType, construct->valType->fields[i].first);
			if (argType != construct->valType->fields[i].first)
				throw std::runtime_error("Struct constructor argument type mismatch.");
			break;
		default:
			throw std::runtime_error("Unsupported constructor type");
		}
	}

	// emit special codes for matrix constructors
	if (nargs == 1)
	{
		if (construct->valType->kind == TypeKind::Vec2)
			emit(CREATE_VEC2_1ARG);
		else if (construct->valType->kind == TypeKind::Vec3)
			emit(CREATE_VEC3_1ARG);
		else if (construct->valType->kind == TypeKind::Mat2)
			emit(CREATE_MAT2_DIAG);
		else if (construct->valType->kind == TypeKind::Mat3)
			emit(CREATE_MAT3_DIAG);
	}

	if (nargs == 3)
	{
		if (construct->valType->kind == TypeKind::Mat3)
			emit(CREATE_MAT3_VEC3);
	}

	return construct->valType;
}

Type Compiler::expressionType(Expression* expr)
{
	if (auto l = dynamic_cast<LiteralExpr*>(expr))
	{
		return prg.types.getBuiltinType(l->value);
	}
	if (auto v = dynamic_cast<VariableExpr*>(expr))
	{
		return resolveVariableType(v->name);
	}
	if (auto u = dynamic_cast<UnaryExpr*>(expr))
	{
		return expressionType(u->right.get());
	}
	if (auto b = dynamic_cast<BinaryExpr*>(expr))
	{
		Type leftType = expressionType(b->left.get());
		Type rightType = expressionType(b->right.get());

		// this is only supported for multiplications!
		if (b->op == BinaryOp::Multiply)
		{
			// scalar * scalar --> scalar
			if (isScalarType(leftType) && isScalarType(rightType))
				return commonType(leftType, rightType);

			// scalar * type --> type
			if (isScalarType(leftType) && isVec2Type(rightType))
				return rightType;
			if (isScalarType(leftType) && isVec3Type(rightType))
				return rightType;
			if (isScalarType(leftType) && isMat2Type(rightType))
				return rightType;
			if (isScalarType(leftType) && isMat3Type(rightType))
				return rightType;

			// vec2 * type
			if (isVec2Type(leftType) && isScalarType(rightType))
				return leftType;
			if (isVec2Type(leftType) && isVec2Type(rightType)) // dot product
				return prg.types.Double();

			// vec3 * type
			if (isVec3Type(leftType) && isScalarType(rightType))
				return leftType;
			if (isVec3Type(leftType) && isVec3Type(rightType)) // dot product
				return prg.types.Double();

			// mat2 * type
			if (isMat2Type(leftType) && isScalarType(rightType))
				return leftType;
			if (isMat2Type(leftType) && isVec2Type(rightType))
				return rightType;
			if (isMat2Type(leftType) && isMat2Type(rightType))
				return leftType;

			// mat3 * type
			if (isMat3Type(leftType) && isScalarType(rightType))
				return leftType;
			if (isMat3Type(leftType) && isVec3Type(rightType))
				return rightType;
			if (isMat3Type(leftType) && isMat3Type(rightType))
				return leftType;
		}
	}
	if (auto c = dynamic_cast<CallExpr*>(expr))
	{
		std::vector<Type> args;
		for (auto& arg : c->arguments)
			args.push_back(expressionType(arg.get()));

		int index = prg.resolveFunction(c->name, args);
		if (index >=0 )
		{
			FunctionInfo& func = prg.functions[index];
			return func.returnType;
		}
	}
	if (auto i = dynamic_cast<ConstructorExpr*>(expr))
	{
		return i->valType;
	}
	throw std::runtime_error("Unsupported expression type in expressionType");
}

void Compiler::compileVarDecl(VarDeclStmt* decl)
{
	Type baseType = decl->type;
	for (auto& var : decl->vars)
	{
		Type type = baseType;
		if (var->arraySizes.size() > 0)
		{
			type = prg.types.getArrayType(baseType, var->arraySizes);
		}

		if (!decl->input && var->initializer)
		{
			Type initType = coerce(compileExpression(var->initializer.get()), type);
		}

		if (m_scopeDepth == 0)
		{
			if (decl->input)
				prg.addInput(var->name, type);
			else
			{
				prg.addGlobal(var->name, type);

				if (var->initializer)
				{
					auto it = prg.globalIndices.find(var->name);

					Program::Global& global = prg.globals[it->second];

					switch (type->kind)
					{
					case TypeKind::Bool  : emit(OpCode::SET_GLOBAL_BOOL  ); break;
					case TypeKind::Int   : emit(OpCode::SET_GLOBAL_INT   ); break;
					case TypeKind::Double: emit(OpCode::SET_GLOBAL_DOUBLE); break;
					case TypeKind::Vec2  : emit(OpCode::SET_GLOBAL_VEC2  ); break;
					case TypeKind::Vec3  : emit(OpCode::SET_GLOBAL_VEC3  ); break;
					case TypeKind::Mat2  : emit(OpCode::SET_GLOBAL_MAT2  ); break;
					case TypeKind::Mat3  : emit(OpCode::SET_GLOBAL_MAT3  ); break;
					case TypeKind::Array :
						emit(OpCode::SET_GLOBAL_ARRAY);
						emitUint8((uint8_t)type->size());
						break;
					case TypeKind::Struct:
						emit(OpCode::SET_GLOBAL_STRUCT);
						emitUint8((uint8_t)type->size());
						break;
					default:
						throw std::runtime_error("Unsupported global variable type");
					}

					emitUint8((uint8_t)global.slot);

					switch (type->kind)
					{
					case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
					case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
					case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
					case TypeKind::Vec2  : emit(OpCode::POP_VEC2  ); break;
					case TypeKind::Vec3  : emit(OpCode::POP_VEC3  ); break;
					case TypeKind::Mat2  : emit(OpCode::POP_MAT2  ); break;
					case TypeKind::Mat3  : emit(OpCode::POP_MAT3  ); break;
					case TypeKind::Array :
						if (type->size() > 255)
							throw std::runtime_error("Array size exceeds maximum of 255.");
						emit(OpCode::POP_ARRAY, type->size());
						emitUint8((uint8_t)type->size());
						break;
					case TypeKind::Struct:
						emit(OpCode::POP_STRUCT, type->size());
						emitUint8((uint8_t)type->typeIndex);
						break;
					default:
						throw std::runtime_error("Unsupported global variable type");
					}
				}
			}
		}
		else
		{
			// make sure the variable isn't already declared in this scope
			for (int i = (int)m_locals.size() - 1; i >= 0; --i)
			{
				if (m_locals[i].depth < m_scopeDepth)
					break;
				if (m_locals[i].name == var->name)
					throw std::runtime_error("Variable '" + var->name + "' is already declared in this scope.");
			}

			m_locals.push_back({ var->name, type, m_scopeDepth, localStackSize });

			localStackSize += (int)type->size();
			if (var->initializer == nullptr)
			{
				emit(OpCode::PUSH_LOCAL, (uint8_t)type->size()); // reserve space on the stack for the variable
				emitUint8((uint8_t)type->size());

				// If an initializer was defined, the global stack effect has already been calculated. 
				stackDepth += (int)type->size();
				if (stackDepth > maxStackDepth)
					maxStackDepth = stackDepth;
			}
		}
	}
}

//
// ===== If / While =====
//

void Compiler::compileIf(IfStmt* stmt)
{
	Type type = compileExpression(stmt->condition.get());
	if (!isScalarType(type) && !isBoolType(type))
		throw std::runtime_error("Condition expression must be of type bool.");

	int thenJump = emitJump(OpCode::JUMP_IF_FALSE);

	switch (type->kind)
	{
	case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
	case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
	default:
		throw std::runtime_error("Unsupported condition type in while statement");
	}

	compileStatement(stmt->thenBranch.get());

	int elseJump = emitJump(OpCode::JUMP);

	patchJump(thenJump);

	switch (type->kind)
	{
	case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
	case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
	default:
		throw std::runtime_error("Unsupported condition type in while statement");
	}

	if (stmt->elseBranch)
		compileStatement(stmt->elseBranch.get());

	patchJump(elseJump);
}

void Compiler::compileWhile(WhileStmt* stmt)
{
	int loopStart = (int)prg.code.size();

	Type type = compileExpression(stmt->condition.get());
	if (!isLogicalType(type))
		throw std::runtime_error("Condition expression must be of a numeric type or bool.");

	int exitJump = emitJump(OpCode::JUMP_IF_FALSE);
	emit(OpCode::POP_BOOL);
	compileStatement(stmt->body.get());

	emitLoop(loopStart);

	patchJump(exitJump);
	emit(OpCode::POP_BOOL);
}

void Compiler::compileFor(ForStmt* stmt)
{
	if (stmt->initializer) compileStatement(stmt->initializer.get());

	int loopStart = (int)prg.code.size();
	Type type = compileExpression(stmt->condition.get());
	if (!isScalarType(type) && type != prg.types.Bool())
		throw std::runtime_error("Condition expression must be of a numeric type or bool.");

	int exitJump = emitJump(OpCode::JUMP_IF_FALSE);
	emit(OpCode::POP_BOOL);

	compileStatement(stmt->body.get());
	Type incrType = compileExpression(stmt->increment.get());
	pop(incrType);

	emitLoop(loopStart);

	patchJump(exitJump);
	emit(OpCode::POP_BOOL);
}

void Compiler::compileReturn(ReturnStmt* stmt)
{
	if (stmt->value)
	{
		Type returnType = expectedReturnType;
		returnType = compileExpression(stmt->value.get());

		if (expectedReturnType)
		{
			Type coerceType = coerce(returnType, expectedReturnType);
			if (coerceType)
				returnType = coerceType;

			if (returnType != expectedReturnType)
				error(stmt->value->location, "Return type mismatch. Expected " + TypeToString(expectedReturnType) + " but got " + TypeToString(returnType));
		}

		switch (returnType->kind)
		{
		case TypeKind::Bool  : emit(OpCode::RETURN_BOOL  ); break;
		case TypeKind::Int   : emit(OpCode::RETURN_INT   ); break;
		case TypeKind::Double: emit(OpCode::RETURN_DOUBLE); break;
		case TypeKind::Vec2  : emit(OpCode::RETURN_VEC2  ); break;
		case TypeKind::Vec3  : emit(OpCode::RETURN_VEC3  ); break;
		case TypeKind::Mat2  : emit(OpCode::RETURN_MAT2  ); break;
		case TypeKind::Mat3  : emit(OpCode::RETURN_MAT3  ); break;
		case TypeKind::Array : 
			emit(OpCode::RETURN_ARRAY ); 
			emitUint8(returnType->typeIndex);
			break;
		case TypeKind::Struct: 
			emit(OpCode::RETURN_STRUCT); 
			emitUint8(returnType->typeIndex);
			break;
		default:
			throw std::runtime_error("Unsupported return type");
		};
	}
	else
	{
		if (expectedReturnType != nullptr && expectedReturnType != prg.types.Void())
			error(stmt, "Missing return value in function with non-void return type.");

		// return without value -> push monostate
		emit(OpCode::PUSH_VOID);
		emit(OpCode::RETURN_VOID);
	}

	// Mark that a return statement was encountered. 
	// This is used to determine whether we need to emit an implicit return at the end of a function/program.
	hasReturn = true;
}

void Compiler::compileStruct(StructStmt* stmt)
{
	// nothing to do here since struct definitions are handled during parsing and type registration.
}

//
// ===== Functions =====
//

void Compiler::compileFunction(FunctionStmt* fn)
{
	// Emit jump over function body
	int jumpOver = emitJump(OpCode::JUMP);

	std::vector<Type> argTypes;
	for (auto& p : fn->params)
		argTypes.push_back(p.first);

	int fnIndex = prg.resolveFunction(fn->name, argTypes);
	if (fnIndex < 0)
		throw std::runtime_error("Compiler error: function not found after resolution.");

	FunctionInfo& info = prg.functions[fnIndex];
	info.entry = prg.code.size();
	currentFunction = fnIndex;

	Type currentReturnType = expectedReturnType;
	expectedReturnType = fn->returnType;
	bool hasReturnBefore = hasReturn;
	hasReturn = false;

	beginScope();

	for (auto& p : fn->params)
	{
		m_locals.push_back({ p.second, p.first, m_scopeDepth, localStackSize });
		localStackSize += (int)p.first->size();
		info.argSize += (int)p.first->size();
		stackDepth += (int)p.first->size();
		if (stackDepth > maxStackDepth)
			maxStackDepth = stackDepth;
	}

	size_t currentStackSize = stackDepth;

	BlockStmt* body = dynamic_cast<BlockStmt*>(fn->body.get());
	for (auto& stmt : body->statements)
		compileStatement(stmt.get());

	if (!hasReturn)
	{
		if (fn->returnType != nullptr && fn->returnType != prg.types.Void())
			throw std::runtime_error("Missing return statement in function with non-void return type.");

		emit(OpCode::RETURN_VOID);
	}

	// We don't call endScope() here since the function body must end with a return, this will never pop the function's locals. 
	// Instead, just clear locals and reset the scope depth.
	while (!m_locals.empty() && m_locals.back().depth == m_scopeDepth)
	{
		Local& local = m_locals.back();
		switch (local.type->kind)
		{
		case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
		case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
		case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
		case TypeKind::Vec2  : emit(OpCode::POP_VEC2  ); break;
		case TypeKind::Vec3  : emit(OpCode::POP_VEC3  ); break;
		case TypeKind::Mat2  : emit(OpCode::POP_MAT2  ); break;
		case TypeKind::Mat3  : emit(OpCode::POP_MAT3  ); break;
		case TypeKind::Array : 
			if (local.type->size() > 255)
				throw std::runtime_error("Array size exceeds maximum of 255.");
			emit(OpCode::POP_ARRAY, local.type->size()); 
			emitUint8(local.type->size());
			break;
		case TypeKind::Struct:
			emit(OpCode::POP_STRUCT, local.type->size());
			emitUint8(local.type->typeIndex);
			break;
		default:
			throw std::runtime_error("Unsupported local variable type");
			break;
		}
		m_locals.pop_back();
		localStackSize -= (int)local.type->size();
	}
	m_scopeDepth--;

	// Patch jump so execution skips function body
	patchJump(jumpOver);

	expectedReturnType = currentReturnType;
	hasReturn = hasReturnBefore;

	info.maxStackSize = maxStackDepth - currentStackSize;

	stackDepth = currentStackSize;

	currentFunction = -1;
}

//
// ===== Expressions =====
//

Type Compiler::compileExpression(Expression* expr)
{
	Type type;
	if      (auto b = dynamic_cast<BinaryExpr*     >(expr)) type = compileBinary(b);
	else if (auto u = dynamic_cast<UnaryExpr*      >(expr)) type = compileUnary(u);
	else if (auto l = dynamic_cast<LiteralExpr*    >(expr)) type = compileLiteral(l);
	else if (auto v = dynamic_cast<VariableExpr*   >(expr)) type = compileVariable(v);
	else if (auto a = dynamic_cast<AssignExpr*     >(expr)) type = compileAssign(a);
	else if (auto c = dynamic_cast<CallExpr*       >(expr)) type = compileCall(c);
	else if (auto m = dynamic_cast<MemberExpr*     >(expr)) type = compileMember(m);
	else if (auto s = dynamic_cast<IndexExpr*      >(expr)) type = compileIndex(s);
	else if (auto i = dynamic_cast<InitExpr*       >(expr)) type = compileInitializer(i);
	else if (auto c = dynamic_cast<ConstructorExpr*>(expr)) type = compileConstructor(c);
	else
		throw std::runtime_error("Unsupported expression type");

	return type;
}

Type Compiler::compileLiteral(LiteralExpr* expr)
{
	Type type = prg.types.getBuiltinType(expr->value);
	switch (type->kind)
	{
	case TypeKind::Bool  : emit(OpCode::PUSH_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::PUSH_INT   ); break;
	case TypeKind::Double: emit(OpCode::PUSH_DOUBLE); break;
	case TypeKind::Vec2  : emit(OpCode::PUSH_VEC2  ); break;
	case TypeKind::Vec3  : emit(OpCode::PUSH_VEC3  ); break;
	case TypeKind::Mat2  : emit(OpCode::PUSH_MAT2  ); break;
	case TypeKind::Mat3  : emit(OpCode::PUSH_MAT3  ); break;
	default:
		throw std::runtime_error("Unsupported literal type");
	}

	uint8_t idx = prg.addConstant(expr->value);
	emitUint8(idx);
	return type;
}

Type Compiler::compileVariable(VariableExpr* expr)
{
	int local = resolveLocal(expr->name);
	if (local != -1)
	{
		Local& localInfo = m_locals[local];
		Type type = localInfo.type;

		switch (type->kind)
		{
		case TypeKind::Bool  : emit(OpCode::GET_LOCAL_BOOL  ); break;
		case TypeKind::Int   : emit(OpCode::GET_LOCAL_INT   ); break;
		case TypeKind::Double: emit(OpCode::GET_LOCAL_DOUBLE); break;
		case TypeKind::Vec2  : emit(OpCode::GET_LOCAL_VEC2  ); break;
		case TypeKind::Vec3  : emit(OpCode::GET_LOCAL_VEC3  ); break;
		case TypeKind::Mat2  : emit(OpCode::GET_LOCAL_MAT2  ); break;
		case TypeKind::Mat3  : emit(OpCode::GET_LOCAL_MAT3  ); break;
		case TypeKind::Array:
		{
			emit(OpCode::GET_LOCAL_ARRAY, type->size()); 
			emitUint8(type->size());
			break;
		}
		case TypeKind::Struct:
		{
			emit(OpCode::GET_LOCAL_STRUCT, type->size());
			emitUint8(type->size());
			break;
		}
		default:
			throw std::runtime_error("Unsupported local variable type");
		}

		emitUint8(localInfo.slot);

		return type;
	}

	int globalIndex = resolveGlobal(expr->name);
	if (globalIndex < 0)
		throw std::runtime_error("Undefined variable: " + expr->name);

	Program::Global& glob = prg.globals[globalIndex];

	glob.refcount++;

	switch (glob.type->kind)
	{
	case TypeKind::Bool  : emit(OpCode::GET_GLOBAL_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::GET_GLOBAL_INT   ); break;
	case TypeKind::Double: emit(OpCode::GET_GLOBAL_DOUBLE); break;
	case TypeKind::Vec2  : emit(OpCode::GET_GLOBAL_VEC2  ); break;
	case TypeKind::Vec3  : emit(OpCode::GET_GLOBAL_VEC3  ); break;
	case TypeKind::Mat2  : emit(OpCode::GET_GLOBAL_MAT2  ); break;
	case TypeKind::Mat3  : emit(OpCode::GET_GLOBAL_MAT3  ); break;
	case TypeKind::Array:
	{
		emit(OpCode::GET_GLOBAL_ARRAY, glob.type->size()); 
		emitUint8(glob.type->size());
		break;
	}
	case TypeKind::Struct:
	{
		emit(OpCode::GET_GLOBAL_STRUCT, glob.type->size());
		emitUint8(glob.type->size());
		break;
	}
	default:
		throw std::runtime_error("Unsupported global variable type");
	}

	emitUint8(glob.slot);
	return glob.type;
}

Type Compiler::compileAssign(AssignExpr* expr)
{
	Type l_type = compileLValue(expr->target.get());
	Type r_type = coerce(compileExpression(expr->value.get()), l_type);

	switch (l_type->kind)
	{
	case TypeKind::Bool  : emit(OpCode::STORE_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::STORE_INT   ); break;
	case TypeKind::Double: emit(OpCode::STORE_DOUBLE); break;
	case TypeKind::Vec2  : emit(OpCode::STORE_VEC2  ); break;
	case TypeKind::Vec3  : emit(OpCode::STORE_VEC3  ); break;
	case TypeKind::Mat2  : emit(OpCode::STORE_MAT2  ); break;
	case TypeKind::Mat3  : emit(OpCode::STORE_MAT3  ); break;
	case TypeKind::Array:
	{
		emit(OpCode::STORE_ARRAY);
		emitUint8(l_type->size());
		break;
	}
	case TypeKind::Struct:
	{
		emit(OpCode::STORE_STRUCT);
		emitUint8(l_type->size());
		break;
	}
	default:
		throw std::runtime_error("Unsupported type in assignment");
		break;
	}

	return l_type;
}

Type Compiler::compileLValue(Expression* expr)
{
	if (auto var = dynamic_cast<VariableExpr*>(expr))
	{
		return compileVariableRef(var);
	}

	if (auto member = dynamic_cast<MemberExpr*>(expr))
	{
		return compileMemberRef(member);
	}

	if (auto index = dynamic_cast<IndexExpr*>(expr))
	{
		return compileIndexRef(index);
	}

	throw std::runtime_error("Invalid assignment target.");
}

Type Compiler::compileVariableRef(VariableExpr* expr)
{
	Type returnType = nullptr;

	int index = resolveLocal(expr->name);
	if (index != -1)
	{
		returnType = m_locals[index].type;

		emit(OpCode::GET_LOCAL_REF);
		emitUint8((uint8_t)m_locals[index].slot);
	}
	else
	{
		index = resolveGlobal(expr->name);
		if (index < 0)
			throw std::runtime_error("Undefined variable: " + expr->name);

		if (prg.globals[index].immutable)
			throw std::runtime_error("Cannot assign to immutable global variable: " + expr->name);

		returnType = prg.globals[index].type;

		emit(OpCode::GET_GLOBAL_REF);
		emitUint8((uint8_t)prg.globals[index].slot);
	}

	return returnType;
}

Type Compiler::compileMemberRef(MemberExpr* expr)
{
	Type objectType = compileLValue(expr->object.get());
	if (isStructType(objectType))
	{
		int memberIndex  = resolveMember(objectType, expr->property);
		int memberOffset = resolveMemberOffset(objectType, expr->property);

		emit(OpCode::GET_MEMBER_REF);
		emitUint8((uint8_t)memberOffset);

		return memberType(objectType, memberIndex);
	}
	else if (isVec2Type(objectType))
	{
		if (expr->property == "x")
		{
			emit(OpCode::GET_VEC2_X_REF);
			return prg.types.Double();
		}
		else if (expr->property == "y")
		{
			emit(OpCode::GET_VEC2_Y_REF);
			return prg.types.Double();
		}
		else
			throw std::runtime_error("Vec2 has no member named '" + expr->property + "'.");
	}
	else if (isVec3Type(objectType))
	{
		if (expr->property == "x")
		{
			emit(OpCode::GET_VEC3_X_REF);
			return prg.types.Double();
		}
		else if (expr->property == "y")
		{
			emit(OpCode::GET_VEC3_Y_REF);
			return prg.types.Double();
		}
		else if (expr->property == "z")
		{
			emit(OpCode::GET_VEC3_Z_REF);
			return prg.types.Double();
		}
		else
			throw std::runtime_error("Vec3 has no member named '" + expr->property + "'.");
	}
	else
	{
		throw std::runtime_error("Left-hand side of member access must be a struct or vector.");
	}
}

int Compiler::resolveMember(Type type, const std::string& member)
{
	assert(type->kind == TypeKind::Struct);
	for (size_t i = 0; i < type->fields.size(); ++i)
	{
		if (type->fields[i].second == member)
			return (int)i;
	}
	throw std::runtime_error("Struct '" + TypeToString(type) + "' has no member named '" + member + "'.");
}

int Compiler::resolveMemberOffset(Type type, const std::string& member)
{
	assert(type->kind == TypeKind::Struct);
	int offset = 0;
	for (size_t i = 0; i < type->fields.size(); ++i)
	{
		if (type->fields[i].second == member)
			return offset;
		offset += (int)type->fields[i].first->size();
	}
	throw std::runtime_error("Struct '" + TypeToString(type) + "' has no member named '" + member + "'.");
}

Type Compiler::memberType(Type type, int memberIndex)
{
	assert(type->kind == TypeKind::Struct);
	if (memberIndex < 0 || memberIndex >= (int)type->fields.size())
		throw std::runtime_error("Invalid member index for struct '" + TypeToString(type) + "'.");
	return type->fields[memberIndex].first;
}

Type Compiler::compileIndexRef(IndexExpr* expr)
{
	Type objectType = compileLValue(expr->object.get());

	if (objectType->kind != TypeKind::Array)
		throw std::runtime_error("Left-hand side of index access must be an array.");

	compileExpression(expr->index.get());

	switch (objectType->elementType->kind)
	{
	case TypeKind::Bool  : emit(OpCode::GET_INDEX_REF_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::GET_INDEX_REF_INT   ); break;
	case TypeKind::Double: emit(OpCode::GET_INDEX_REF_DOUBLE); break;
	case TypeKind::Vec2  : emit(OpCode::GET_INDEX_REF_VEC2  ); break;
	case TypeKind::Vec3  : emit(OpCode::GET_INDEX_REF_VEC3  ); break;
	case TypeKind::Mat2  : emit(OpCode::GET_INDEX_REF_MAT2  ); break;
	case TypeKind::Mat3  : emit(OpCode::GET_INDEX_REF_MAT3  ); break;
	case TypeKind::Array : 
	case TypeKind::Struct: 
		emit(OpCode::GET_INDEX_REF);
		emitUint8(objectType->elementType->size());
		break;
	default:
		throw std::runtime_error("Unsupported array element type in index access");
	}

	return objectType->elementType;
}

//
// ===== Binary (includes &&, ||) =====
//

Type Compiler::compileBinary(BinaryExpr* expr)
{
	BinaryOp op = expr->op;

	// Short-circuit AND
	if (op == BinaryOp::AndAnd)
	{
		Type type = compileExpression(expr->left.get());

		if (!isScalarType(type) && !isBoolType(type))
			throw std::runtime_error("Cannot convert left operand of '&&' to boolean.");

		int endJump = emitJump(OpCode::JUMP_IF_FALSE);

		switch (type->kind)
		{
		case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
		case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
		case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
		default:
			throw std::runtime_error("Unsupported condition type in '&&' operator");
		}

		type = compileExpression(expr->right.get());
		if (!isScalarType(type) && !isBoolType(type))
			throw std::runtime_error("Cannot convert right operand of '&&' to boolean.");

		patchJump(endJump);
		return prg.types.Bool();
	}

	// Short-circuit OR
	if (op == BinaryOp::OrOr)
	{
		Type type = compileExpression(expr->left.get());

		if (!isScalarType(type) && !isBoolType(type))
			throw std::runtime_error("Cannot convert left operand of '||' to boolean.");

		int endJump = emitJump(OpCode::JUMP_IF_TRUE);
		switch (type->kind)
		{
		case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
		case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
		case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
		default:
			throw std::runtime_error("Unsupported condition type in || operator");
		}

		type = compileExpression(expr->right.get());
		if (!isScalarType(type) && !isBoolType(type))
			throw std::runtime_error("Cannot convert right operand of '||' to boolean.");

		patchJump(endJump);
		return prg.types.Bool();
	}

	// Handle special cases first
	// see if both expressions are global variables
	if (isVariable(expr->left) && isVariable(expr->right))
	{
		int globA = resolveGlobal(dynamic_cast<VariableExpr*>(expr->left.get())->name);
		int globB = resolveGlobal(dynamic_cast<VariableExpr*>(expr->right.get())->name);

		if ((globA >= 0) && (globB >= 0))
		{
			Program::Global& a = prg.globals[globA];
			Program::Global& b = prg.globals[globB];

			if (isMat3Type(a.type) && isMat3Type(b.type))
			{
				switch (op)
				{
				case BinaryOp::Plus    : emit(OpCode::ADD_GLOBAL_MAT3); break;
				case BinaryOp::Minus   : emit(OpCode::SUB_GLOBAL_MAT3); break;
				case BinaryOp::Multiply: emit(OpCode::MUL_GLOBAL_MAT3); break;
				default:
					throw std::runtime_error("Unsupported binary op for mat3 type.");
				}

				emitUint8((uint8_t)a.slot);
				emitUint8((uint8_t)b.slot);
				return a.type;
			}
		}
	}


	Type type_l = compileExpression(expr->left.get());

	if (isDoubleType(type_l) && isLiteral(expr->right) && (expr->op == BinaryOp::Exponent))
	{
		Value e = dynamic_cast<LiteralExpr*>(expr->right.get())->value;
		if (isDouble(e))
		{
			// If the exponent is a literal, we can optimize by using a fast exponentiation algorithm.
			double exponentValue = e.d;
			if (exponentValue == 0.5)
			{
				emit(OpCode::SQRT_DOUBLE);
				return type_l;
			}
			else if (exponentValue == 2.0)
			{
				emit(OpCode::SQR_DOUBLE);
				return type_l;
			}
		}
		else if (isInt(e))
		{
			int exponentValue = e.i;
			if (exponentValue == 2)
			{
				emit(OpCode::SQR_DOUBLE);
				return type_l;
			}
		}
	}

	Type type_r = compileExpression(expr->right.get());
	Type t = commonType(type_l, type_r);
	if (t)
	{
		type_l = t;
		type_r = t;
	}

	if (isBoolType(type_l) && isBoolType(type_r))
	{
		switch (op)
		{
		case BinaryOp::EqualEqual: emit(OpCode::EQUAL_BOOL); type_l = prg.types.Bool(); break;
		case BinaryOp::NotEqual  : emit(OpCode::NEQ_BOOL  ); type_l = prg.types.Bool(); break;
		default: throw std::runtime_error("Unsupported binary op for bool type.");
		}
		return type_l;
	}

	if (isIntType(type_l) && isIntType(type_r))
	{
		switch (op)
		{
		case BinaryOp::Plus    : emit(OpCode::ADD_INT); break;
		case BinaryOp::Minus   : emit(OpCode::SUB_INT); break;
		case BinaryOp::Multiply: emit(OpCode::MUL_INT); break;
		case BinaryOp::Divide  : emit(OpCode::DIV_INT); break;
		case BinaryOp::Exponent: emit(OpCode::EXP_INT); break;
		case BinaryOp::Greater     : emit(OpCode::GT_INT   ); type_l = prg.types.Bool(); break;
		case BinaryOp::Less        : emit(OpCode::LT_INT   ); type_l = prg.types.Bool(); break;
		case BinaryOp::GreaterEqual: emit(OpCode::GE_INT   ); type_l = prg.types.Bool(); break;
		case BinaryOp::LessEqual   : emit(OpCode::LE_INT   ); type_l = prg.types.Bool(); break;
		case BinaryOp::EqualEqual  : emit(OpCode::EQUAL_INT); type_l = prg.types.Bool(); break;
		case BinaryOp::NotEqual    : emit(OpCode::NEQ_INT  ); type_l = prg.types.Bool(); break;
		default: throw std::runtime_error("Unsupported binary op for int type.");
		}
		return type_l;
	}

	if (isDoubleType(type_l) && isDoubleType(type_r))
	{
		switch (op)
		{
		case BinaryOp::Plus    : emit(OpCode::ADD_DOUBLE); break;
		case BinaryOp::Minus   : emit(OpCode::SUB_DOUBLE); break;
		case BinaryOp::Multiply: emit(OpCode::MUL_DOUBLE); break;
		case BinaryOp::Divide  : emit(OpCode::DIV_DOUBLE); break;
		case BinaryOp::Exponent: emit(OpCode::EXP_DOUBLE); break;
		case BinaryOp::Greater     : emit(OpCode::GT_DOUBLE   ); type_l = prg.types.Bool(); break;
		case BinaryOp::Less        : emit(OpCode::LT_DOUBLE   ); type_l = prg.types.Bool(); break;
		case BinaryOp::GreaterEqual: emit(OpCode::GE_DOUBLE   ); type_l = prg.types.Bool(); break;
		case BinaryOp::LessEqual   : emit(OpCode::LE_DOUBLE   ); type_l = prg.types.Bool(); break;
		case BinaryOp::EqualEqual  : emit(OpCode::EQUAL_DOUBLE); type_l = prg.types.Bool(); break;
		case BinaryOp::NotEqual    : emit(OpCode::NEQ_DOUBLE  ); type_l = prg.types.Bool(); break;
		default: throw std::runtime_error("Unsupported binary op for double type.");
		}
		return type_l;
	}

	// vec2 operations
	if (isVec2Type(type_l) && isVec2Type(type_r))
	{
		Type returnType = type_l; // default return type is vec2, but some ops like dot product will return double
		switch (op)
		{
		case BinaryOp::Plus    : emit(OpCode::ADD_VEC2); break;
		case BinaryOp::Minus   : emit(OpCode::SUB_VEC2); break;
		case BinaryOp::Multiply: emit(OpCode::DOT_VEC2); returnType = prg.types.Double(); break;
		default: throw std::runtime_error("Unsupported binary op for vec2 type.");
		}
		return returnType;
	}

	// vec2 op scalar
	if (isVec2Type(type_l) && isScalarType(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_VEC2_DOUBLE); break;
		case BinaryOp::Divide : emit(OpCode::DIV_VEC2_DOUBLE); break;
		default: throw std::runtime_error("Unsupported binary op for vec2 and double types.");
		}
		return type_l;
	}

	// scalar * vec2
	if (isScalarType(type_l) && isVec2Type(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_DOUBLE_VEC2); break;
		default: throw std::runtime_error("Unsupported binary op for double and vec2 types.");
		}
		return type_r;
	}

	// vec3 operations
	if (isVec3Type(type_l) && isVec3Type(type_r))
	{
		Type returnType = type_l; // default return type is vec3, but some ops like dot product will return double
		switch (op)
		{
		case BinaryOp::Plus : emit(OpCode::ADD_VEC3); break;
		case BinaryOp::Minus: emit(OpCode::SUB_VEC3); break;
		case BinaryOp::Multiply: emit(OpCode::DOT_VEC3); returnType = prg.types.Double(); break;
		default: throw std::runtime_error("Unsupported binary op for vec3 type.");
		}
		return returnType;
	}

	// vec3 op scalar
	if (isVec3Type(type_l) && isScalarType(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_VEC3_DOUBLE); break;
		case BinaryOp::Divide  : emit(OpCode::DIV_VEC3_DOUBLE); break;
		default: throw std::runtime_error("Unsupported binary op for vec3 and double types.");
		}
		return type_l;
	}

	// scalar * vec3
	if (isScalarType(type_l) && isVec3Type(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_DOUBLE_VEC3); break;
		default: throw std::runtime_error("Unsupported binary op for double and vec3 types.");
		}
		return type_r;
	}

	// mat2 operations
	if (isMat2Type(type_l) && isMat2Type(type_r))
	{
		Type returnType = type_l;
		switch (op)
		{
		case BinaryOp::Plus    : emit(OpCode::ADD_MAT2); break;
		case BinaryOp::Minus   : emit(OpCode::SUB_MAT2); break;
		case BinaryOp::Multiply: emit(OpCode::MUL_MAT2); break;
		default: throw std::runtime_error("Unsupported binary op for mat2 type.");
		}
		return returnType;
	}

	// mat2 * vec2
	if (isMat2Type(type_l) && isVec2Type(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_MAT2_VEC2); break;
		default: throw std::runtime_error("Unsupported binary op for mat2 and vec2 types.");
		}
		return type_r;
	}

	// mat2 op scalar
	if (isMat2Type(type_l) && isScalarType(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_MAT2_DOUBLE); break;
		case BinaryOp::Divide  : emit(OpCode::DIV_MAT2_DOUBLE); break;
		default: throw std::runtime_error("Unsupported binary op for mat2 and double types.");
		}
		return type_l;
	}

	// scalar * mat2
	if (isScalarType(type_l) && isMat2Type(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_DOUBLE_MAT2); break;
		default: throw std::runtime_error("Unsupported binary op for double and mat2 types.");
		}
		return type_r;
	}

	// mat3 operations
	if (isMat3Type(type_l) && isMat3Type(type_r))
	{
		Type returnType = type_l;
		switch (op)
		{
		case BinaryOp::Plus    : emit(OpCode::ADD_MAT3); break;
		case BinaryOp::Minus   : emit(OpCode::SUB_MAT3); break;
		case BinaryOp::Multiply: emit(OpCode::MUL_MAT3); break;
		default: throw std::runtime_error("Unsupported binary op for mat3 type.");
		}
		return returnType;
	}

	// mat3 * scalar
	if (isMat3Type(type_l) && isScalarType(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_MAT3_DOUBLE); break;
		case BinaryOp::Divide  : emit(OpCode::DIV_MAT3_DOUBLE); break;
		default: throw std::runtime_error("Unsupported binary op for mat3 and double types.");
		}
		return type_l;
	}

	// scalar * mat3
	if (isScalarType(type_l) && isMat3Type(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_DOUBLE_MAT3); break;
		default: throw std::runtime_error("Unsupported binary op for double and mat3 types.");
		}
		return type_r;
	}

	// mat3 * vec3
	if (isMat3Type(type_l) && isVec3Type(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_MAT3_VEC3); break;
		default: throw std::runtime_error("Unsupported binary op for mat3 and vec3 types.");
		}
		return type_r;
	}

	throw std::runtime_error("Compiler error: Unsupported binary operator '" + opToString(op) + "' with operand types '" + TypeToString(type_l) + "' and '" + TypeToString(type_r) + "'.");

	return type_l;
}

Type Compiler::compileUnary(UnaryExpr* expr)
{
	Type type = compileExpression(expr->right.get());
	switch (expr->op)
	{
	case UnaryOp::Negate:
	{
		if      (isIntType   (type)) emit(OpCode::NEG_INT);
		else if (isDoubleType(type)) emit(OpCode::NEG_DOUBLE);
		else if (isVec2Type  (type)) emit(OpCode::NEG_VEC2);
		else if (isVec3Type  (type)) emit(OpCode::NEG_VEC3);
		else if (isMat2Type  (type)) emit(OpCode::NEG_MAT2);
		else if (isMat3Type  (type)) emit(OpCode::NEG_MAT3);
		else
			throw std::runtime_error("Invalid operand type for unary '-'.");
		break;
	}
	case UnaryOp::Not:
	{
		if (isBoolType(type))
		{
			emit(OpCode::NOT);
			type = prg.types.Bool();
		}
		else
			throw std::runtime_error("Invalid operand type for unary '!'.");
		break;
	}
	default: throw std::runtime_error("Unsupported unary op");
	}
	return type;
}

bool Compiler::isNativeFunction(const std::string& name)
{
	for (int i = 0; i < (int)prg.functions.size(); ++i)
	{
		if (prg.functions[i].name == name)
		{
			return prg.functions[i].isNative;
		}
	}
	return false;
}

std::vector<Type> Compiler::compileFncArgs(std::vector<std::unique_ptr<Expression>>& args)
{
	std::vector<Type> argTypes;
	for (auto& arg : args)
	{
		Type type = compileExpression(arg.get());
		argTypes.push_back(type);
	}
	return argTypes;
}

Type Compiler::compileCall(CallExpr* call)
{
	// don't copy args for native functions
	std::vector<Type> argTypes = compileFncArgs(call->arguments);
	int fnIndex = prg.resolveFunction(call->name, argTypes);
	if (fnIndex == -1)
		throw std::runtime_error("Undefined function: " + call->name);

	// We don't allow recursive calls.
	if (currentFunction == fnIndex)
		throw std::runtime_error("Recursive calls are not allowed: " + call->name);

	const FunctionInfo& fi = prg.functions[fnIndex];

	// check arg count
	if (fi.args.size() != (int)call->arguments.size())
		throw std::runtime_error("Argument count mismatch in call to function: " + call->name);

	// check arg types
	for (int i = 0; i < argTypes.size(); ++i)
	{
		if (argTypes[i] != fi.args[i])
			throw std::runtime_error("Argument type mismatch in call to function: " + call->name);
	}

	stackDepth += (int)fi.maxStackSize;

	// calculate stack size needed for arguments
	int argStackSize = 0;
	for (Type argType : argTypes)
	{
		argStackSize += (int)argType->size();
	}

	// stack size needed for return value
	int returnStackSize = (int)prg.functions[fnIndex].returnType->size();

	emit(OpCode::CALL, returnStackSize - argStackSize);
	emitUint8(fnIndex);
	emitUint8((uint8_t)call->arguments.size());

	return fi.returnType;
}

Type Compiler::compileMember(MemberExpr* expr)
{
	Type type = compileExpression(expr->object.get());

	if (isStructType(type))
	{
		auto it = std::find_if(type->fields.begin(), type->fields.end(),
			[&](const auto& field) { return field.second == expr->property; });
		if (it == type->fields.end())
			throw std::runtime_error("Unknown property: " + expr->property);

		std::size_t index = it - type->fields.begin();

		switch (it->first->kind)
		{
		case TypeKind::Bool  : emit(OpCode::GET_PROPERTY_BOOL  ); break;
		case TypeKind::Int   : emit(OpCode::GET_PROPERTY_INT   ); break;
		case TypeKind::Double: emit(OpCode::GET_PROPERTY_DOUBLE); break;
		case TypeKind::Vec2  : emit(OpCode::GET_PROPERTY_VEC2  ); break;
		case TypeKind::Vec3  : emit(OpCode::GET_PROPERTY_VEC3  ); break;
		case TypeKind::Mat2  : emit(OpCode::GET_PROPERTY_MAT2  ); break;
		case TypeKind::Mat3  : emit(OpCode::GET_PROPERTY_MAT3  ); break;
		case TypeKind::Array : emit(OpCode::GET_PROPERTY_ARRAY ); break;
		case TypeKind::Struct: emit(OpCode::GET_PROPERTY_STRUCT); break;
		default:
			throw std::runtime_error("Unsupported struct member type");
		}

		emitUint8((uint8_t)type->typeIndex);
		emitUint8((uint8_t)index);

		return it->first;
	}
	else if (isVec2Type(type))
	{
		if (expr->property == "x")
		{
			emit(OpCode::GET_VEC2_X);
			return prg.types.Double();
		}
		else if (expr->property == "y")
		{
			emit(OpCode::GET_VEC2_Y);
			return prg.types.Double();
		}
		else
		{
			std::string swizzle = expr->property;
			if (((swizzle.size() == 2) || (swizzle.size() == 3)) &&
				swizzle.find_first_not_of("xy") == std::string::npos)
			{
				uint8_t size = (uint8_t)swizzle.size();
				uint8_t mask = 0;
				for (char c : swizzle)
				{
					switch (c)
					{
					case 'x': mask = (mask << 2) | 0b00; break;
					case 'y': mask = (mask << 2) | 0b01; break;
					default:
						break;
					}

				}
				emit(OpCode::GET_VEC2_SWIZZLE);
				emitUint8(mask);
				emitUint8(size);
				return (size == 2 ? prg.types.Vec2() : prg.types.Vec3());
			}
			else
				throw std::runtime_error("Unknown property: " + expr->property);
		}
	}
	else if (isVec3Type(type))
	{
		if (expr->property == "x")
		{
			emit(OpCode::GET_VEC3_X);
			return prg.types.Double();
		}
		else if (expr->property == "y")
		{
			emit(OpCode::GET_VEC3_Y);
			return prg.types.Double();
		}
		else if (expr->property == "z")
		{
			emit(OpCode::GET_VEC3_Z);
			return prg.types.Double();
		}
		else
		{
			std::string swizzle = expr->property;
			if (((swizzle.size() == 2) || (swizzle.size() == 3)) &&
				swizzle.find_first_not_of("xyz") == std::string::npos)
			{
				uint8_t size = (uint8_t)swizzle.size();
				uint8_t mask = 0;
				for (char c : swizzle)
				{
					switch (c)
					{
					case 'x': mask = (mask << 2) | 0b00; break;
					case 'y': mask = (mask << 2) | 0b01; break;
					case 'z': mask = (mask << 2) | 0b10; break;
					default:
						break;
					}

				}
				emit(OpCode::GET_VEC3_SWIZZLE);
				emitUint8(mask);
				emitUint8(size);
				return (size == 2 ? prg.types.Vec2() : prg.types.Vec3());
			}
			else
				throw std::runtime_error("Unknown property: " + expr->property);
		}
	}
	else
	{
		throw std::runtime_error("Cannot access property of non-struct type.");
	}
}

Type Compiler::compileIndex(IndexExpr* expr)
{
	// handle special cases first
	if (isVariable(expr->object))
	{
		int globIndex = resolveGlobal(dynamic_cast<VariableExpr*>(expr->object.get())->name);
		if (globIndex >= 0)
		{
			Type globType = prg.globals[globIndex].type;
			if (isArrayType(globType) && isDoubleType(globType->elementType))
			{
				Type indxType = compileExpression(expr->index.get());
				if (indxType->kind != TypeKind::Int)
					throw std::runtime_error("array index must be an integer.");
				emit(OpCode::GET_GLOBAL_INDEX_DOUBLE);
				emitUint8((uint8_t)prg.globals[globIndex].slot);
				return prg.types.Double();
			}
		}
	}

	Type exprType = compileExpression(expr->object.get());
	Type indxType = compileExpression(expr->index.get());

	// vec2 indexing
	if (exprType->kind == TypeKind::Vec2)
	{
		if (indxType->kind != TypeKind::Int)
			throw std::runtime_error("vec2 index must be an integer.");

		emit(OpCode::GET_VEC2_INDEX);
		return prg.types.Double();
	}

	// vec3 indexing
	if (exprType->kind == TypeKind::Vec3)
	{
		if (indxType->kind != TypeKind::Int)
			throw std::runtime_error("vec3 index must be an integer.");

		emit(OpCode::GET_VEC3_INDEX);
		return prg.types.Double();
	}

	// mat2 indexing
	if (exprType->kind == TypeKind::Mat2)
	{
		if (indxType->kind != TypeKind::Int)
			throw std::runtime_error("mat2 index must be an integer.");

		emit(OpCode::GET_MAT2_INDEX);
		return prg.types.Vec2();
	}

	// mat3 indexing
	if (exprType->kind == TypeKind::Mat3)
	{
		if (indxType->kind != TypeKind::Int)
			throw std::runtime_error("mat3 index must be an integer.");

		emit(OpCode::GET_MAT3_INDEX);
		return prg.types.Vec3();
	}

	if (exprType->kind != TypeKind::Array)
		throw std::runtime_error("Cannot index non-array type.");
	if (indxType->kind != TypeKind::Int)
		throw std::runtime_error("Array index must be a number.");

	int stackEffect = -(int)exprType->size() + exprType->elementType->size();

	switch (exprType->elementType->kind)
	{
	case TypeKind::Bool  : emit(OpCode::GET_INDEX_BOOL  , stackEffect); break;
	case TypeKind::Int   : emit(OpCode::GET_INDEX_INT   , stackEffect); break;
	case TypeKind::Double: emit(OpCode::GET_INDEX_DOUBLE, stackEffect); break;
	case TypeKind::Vec2  : emit(OpCode::GET_INDEX_VEC2  , stackEffect); break;
	case TypeKind::Vec3  : emit(OpCode::GET_INDEX_VEC3  , stackEffect); break;
	case TypeKind::Mat2  : emit(OpCode::GET_INDEX_MAT2  , stackEffect); break;
	case TypeKind::Mat3  : emit(OpCode::GET_INDEX_MAT3  , stackEffect); break;
	case TypeKind::Array : emit(OpCode::GET_INDEX_ARRAY , stackEffect); break;
	case TypeKind::Struct: emit(OpCode::GET_INDEX_STRUCT, stackEffect); break;
		default:
		throw std::runtime_error("Unsupported array element type");
	}

	emitUint8(exprType->typeIndex);

	return exprType->elementType;
}


const char* febcode::OpCodeToString(febcode::OpCode op)
{
	switch (op)
	{
	case OpCode::PUSH_VOID  : 
	case OpCode::PUSH_BOOL  :
	case OpCode::PUSH_INT   :
	case OpCode::PUSH_DOUBLE:
	case OpCode::PUSH_VEC2  :
	case OpCode::PUSH_VEC3  :
	case OpCode::PUSH_MAT2  :
	case OpCode::PUSH_MAT3  :
		return "MOV ";

	case OpCode::PUSH_LOCAL:
		return "PSHL";

	case OpCode::GET_GLOBAL_BOOL  :
	case OpCode::GET_GLOBAL_INT   :
	case OpCode::GET_GLOBAL_DOUBLE:
	case OpCode::GET_GLOBAL_VEC2  :
	case OpCode::GET_GLOBAL_VEC3  : 
	case OpCode::GET_GLOBAL_MAT2  : 
	case OpCode::GET_GLOBAL_MAT3  : 
	case OpCode::GET_GLOBAL_ARRAY : 
	case OpCode::GET_GLOBAL_STRUCT: return "GTG ";

	case OpCode::SET_GLOBAL_BOOL  :
	case OpCode::SET_GLOBAL_INT   :
	case OpCode::SET_GLOBAL_DOUBLE:
	case OpCode::SET_GLOBAL_VEC2  :
	case OpCode::SET_GLOBAL_VEC3  :
	case OpCode::SET_GLOBAL_MAT2  :
	case OpCode::SET_GLOBAL_MAT3  :
	case OpCode::SET_GLOBAL_ARRAY :
	case OpCode::SET_GLOBAL_STRUCT: return "STG ";

	case OpCode::STORE_BOOL  :
	case OpCode::STORE_INT   :
	case OpCode::STORE_DOUBLE:
	case OpCode::STORE_VEC2  :
	case OpCode::STORE_VEC3  :
	case OpCode::STORE_MAT2:
	case OpCode::STORE_MAT3:
	case OpCode::STORE_ARRAY :
	case OpCode::STORE_STRUCT:
		return "STOR";

	case OpCode::GET_GLOBAL_REF: 
		return "GREF";

	case OpCode::GET_LOCAL_REF       :
		return "LREF";

	case OpCode::GET_INDEX_REF       : 
	case OpCode::GET_INDEX_REF_BOOL  : 
	case OpCode::GET_INDEX_REF_INT   : 
	case OpCode::GET_INDEX_REF_DOUBLE: 
	case OpCode::GET_INDEX_REF_VEC2  : 
	case OpCode::GET_INDEX_REF_VEC3  : 
	case OpCode::GET_INDEX_REF_MAT2  :
	case OpCode::GET_INDEX_REF_MAT3  :
		return "IREF";

	case OpCode::GET_MEMBER_REF: return "MREF";

	case OpCode::GET_LOCAL_BOOL  :
	case OpCode::GET_LOCAL_INT   :
	case OpCode::GET_LOCAL_DOUBLE:
	case OpCode::GET_LOCAL_VEC2  :
	case OpCode::GET_LOCAL_VEC3  :
	case OpCode::GET_LOCAL_MAT2  :
	case OpCode::GET_LOCAL_MAT3  :
	case OpCode::GET_LOCAL_ARRAY :
	case OpCode::GET_LOCAL_STRUCT: return "GETL";

	case OpCode::ADD_INT       : return "ADDI";
	case OpCode::ADD_DOUBLE    : return "ADDF";
	case OpCode::SUB_INT       : return "SUBI";
	case OpCode::SUB_DOUBLE    : return "SUBF";
	case OpCode::MUL_INT       : return "MULI";
	case OpCode::MUL_DOUBLE    : return "MULF";
	case OpCode::DIV_INT       : return "DIVI";
	case OpCode::DIV_DOUBLE    : return "DIVF";
	case OpCode::EXP_INT       : return "EXPI";
	case OpCode::EXP_DOUBLE    : return "EXPF";
	case OpCode::SQR_DOUBLE    : return "SQR ";
	case OpCode::SQRT_DOUBLE   : return "SQRT";

	case OpCode::EQUAL_BOOL  :
	case OpCode::EQUAL_INT   :
	case OpCode::EQUAL_DOUBLE: return "EQ  ";

	case OpCode::NEQ_BOOL  :
	case OpCode::NEQ_INT   :
	case OpCode::NEQ_DOUBLE: return "NEQ ";

	case OpCode::GT_INT        : return "GTI ";
	case OpCode::GT_DOUBLE     : return "GTF ";
	case OpCode::LT_INT        : return "LTI ";
	case OpCode::LT_DOUBLE     : return "LTF ";
	case OpCode::GE_INT        : return "GEI ";
	case OpCode::GE_DOUBLE     : return "GEF ";
	case OpCode::LE_INT        : return "LEI ";
	case OpCode::LE_DOUBLE     : return "LEF ";
	case OpCode::NEG_INT       : return "NEGI";
	case OpCode::NEG_DOUBLE    : return "NEGF";

	case OpCode::CREATE_VEC2_1ARG: return "V2_1";
	case OpCode::GET_VEC2_X    : return "GV2X";
	case OpCode::GET_VEC2_Y    : return "GV2Y";
	case OpCode::GET_VEC2_SWIZZLE: return "G2SW";
	case OpCode::GET_VEC2_INDEX: return "G2I ";
	case OpCode::ADD_VEC2      : return "ADD2";
	case OpCode::SUB_VEC2      : return "SUB2";
	case OpCode::DOT_VEC2      : return "DOT2";
	case OpCode::MUL_VEC2_DOUBLE: return "ML2F";
	case OpCode::MUL_DOUBLE_VEC2: return "MLF2";
	case OpCode::DIV_VEC2_DOUBLE: return "DV2F";
	case OpCode::NEG_VEC2      : return "NEG2";

	case OpCode::CREATE_VEC3_1ARG: return "V3_1";
	case OpCode::GET_VEC3_X    : return "GV3X";
	case OpCode::GET_VEC3_Y    : return "GV3Y";
	case OpCode::GET_VEC3_Z    : return "GV3Z";
	case OpCode::GET_VEC3_SWIZZLE: return "G3SW";
	case OpCode::GET_VEC3_INDEX: return "G3I ";
	case OpCode::ADD_VEC3      : return "ADD3";
	case OpCode::SUB_VEC3      : return "SUB3";
	case OpCode::DOT_VEC3      : return "DOT3";
	case OpCode::MUL_VEC3_DOUBLE: return "ML3F";
	case OpCode::MUL_DOUBLE_VEC3: return "MLF3";
	case OpCode::DIV_VEC3_DOUBLE: return "DV3F";
	case OpCode::NEG_VEC3      : return "NEG3";

	case OpCode::ADD_MAT2       : return "ADM2";
	case OpCode::SUB_MAT2       : return "SBM2";
	case OpCode::MUL_MAT2       : return "MLM2";
	case OpCode::MUL_MAT2_DOUBLE: return "MM2F";
	case OpCode::MUL_DOUBLE_MAT2: return "FM2 ";
	case OpCode::MUL_MAT2_VEC2  : return "M2V2";
	case OpCode::DIV_MAT2_DOUBLE: return "DM2F";
	case OpCode::NEG_MAT2       : return "NGM2";
	case OpCode::GET_MAT2_INDEX : return "GM2I";
	case OpCode::CREATE_MAT2_DIAG: return "MAT2";

	case OpCode::ADD_MAT3       : return "ADM3";
	case OpCode::SUB_MAT3       : return "SBM3";
	case OpCode::MUL_MAT3       : return "MLM3";
	case OpCode::MUL_MAT3_DOUBLE: return "MM3F";
	case OpCode::DIV_MAT3_DOUBLE: return "DM3F";
	case OpCode::MUL_DOUBLE_MAT3: return "FM3 ";
	case OpCode::MUL_MAT3_VEC3  : return "M3V3";
	case OpCode::NEG_MAT3       : return "NGM3";
	case OpCode::GET_MAT3_INDEX : return "GM3I";
	case OpCode::ADD_GLOBAL_MAT3: return "ADM3";
	case OpCode::SUB_GLOBAL_MAT3: return "SBM3";
	case OpCode::MUL_GLOBAL_MAT3: return "MLM3";
	case OpCode::CREATE_MAT3_DIAG: return "MAT3";
	case OpCode::CREATE_MAT3_VEC3: return "C3V3";

	case OpCode::NOT           : return "NOT ";

	case OpCode::GET_PROPERTY_BOOL:
	case OpCode::GET_PROPERTY_INT:
	case OpCode::GET_PROPERTY_DOUBLE:
	case OpCode::GET_PROPERTY_VEC2:
	case OpCode::GET_PROPERTY_VEC3:
	case OpCode::GET_PROPERTY_MAT2:
	case OpCode::GET_PROPERTY_MAT3:
	case OpCode::GET_PROPERTY_ARRAY:
	case OpCode::GET_PROPERTY_STRUCT:
		return "GETP";

	case OpCode::GET_INDEX_BOOL  :
	case OpCode::GET_INDEX_INT   :
	case OpCode::GET_INDEX_DOUBLE:
	case OpCode::GET_INDEX_VEC2  :
	case OpCode::GET_INDEX_VEC3  :
	case OpCode::GET_INDEX_ARRAY :
	case OpCode::GET_INDEX_STRUCT:
		return "GETI";

	case OpCode::GET_GLOBAL_INDEX_DOUBLE:
		return "GETI";

	case OpCode::JUMP          : return "JMP ";
	case OpCode::JUMP_IF_FALSE : return "JMPF";
	case OpCode::JUMP_IF_TRUE  : return "JMPT";
	case OpCode::LOOP          : return "LOOP";
	case OpCode::CALL          : return "CALL";

	case OpCode::POP_VOID  : return "POPV";
	case OpCode::POP_BOOL  : return "POPB";
	case OpCode::POP_INT   : return "POPI";
	case OpCode::POP_DOUBLE: return "POPD";
	case OpCode::POP_VEC2  : return "POP2";
	case OpCode::POP_VEC3  : return "POP3";
	case OpCode::POP_MAT2  : return "PPM2";
	case OpCode::POP_MAT3  : return "PPM3";
	case OpCode::POP_ARRAY : return "POPA";
	case OpCode::POP_STRUCT: return "POPS";	

	case OpCode::GET_VEC2_X_REF: return "RV2X";
	case OpCode::GET_VEC2_Y_REF: return "RV2Y";
	case OpCode::GET_VEC3_X_REF: return "RV3X";
	case OpCode::GET_VEC3_Y_REF: return "RV3Y";
	case OpCode::GET_VEC3_Z_REF: return "RV3Z";

	case OpCode::RETURN_VOID: 
	case OpCode::RETURN_BOOL: 
	case OpCode::RETURN_INT: 
	case OpCode::RETURN_DOUBLE: 
	case OpCode::RETURN_VEC2: 
	case OpCode::RETURN_VEC3: 
	case OpCode::RETURN_MAT2:
	case OpCode::RETURN_MAT3:
	case OpCode::RETURN_ARRAY: 
	case OpCode::RETURN_STRUCT: 
		return "RET ";
	}
	return "(UNKNOWN)";
}
