#include "compiler.h"
#include "scanner.h"
#include "parser.h"

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

	// 3. Compile -> bytecode
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
		case TypeKind::Array : 
			emit(OpCode::POP_ARRAY, local.type->size());
			emitUint8(local.type->typeIndex);
			break;
		case TypeKind::Struct:
			emit(OpCode::POP_STRUCT, local.type->size());
			emitUint8(local.type->typeIndex);
		break;		
		default:
			throw std::runtime_error("Unknown type kind in endScope.");
			break;
		}

		m_locals.pop_back();
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

	throw std::runtime_error("Undefined global variable: " + name);
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

	case OpCode::GET_GLOBAL_BOOL  : return +1;
	case OpCode::GET_GLOBAL_INT   : return +1;
	case OpCode::GET_GLOBAL_DOUBLE: return +1;
	case OpCode::GET_GLOBAL_VEC2  : return +2;
	case OpCode::GET_GLOBAL_VEC3  : return +3;
	case OpCode::GET_GLOBAL_ARRAY : return +arg;
	case OpCode::GET_GLOBAL_STRUCT: return +arg;

	case OpCode::SET_GLOBAL_BOOL: 
	case OpCode::SET_GLOBAL_INT: 
	case OpCode::SET_GLOBAL_DOUBLE: 
	case OpCode::SET_GLOBAL_VEC2: 
	case OpCode::SET_GLOBAL_VEC3: 
	case OpCode::SET_GLOBAL_ARRAY: 
	case OpCode::SET_GLOBAL_STRUCT: 
		return 0;

	case OpCode::GET_GLOBAL_REF       : return +1;

	case OpCode::GET_LOCAL_BOOL  : return +1;
	case OpCode::GET_LOCAL_INT   : return +1;
	case OpCode::GET_LOCAL_DOUBLE: return +1;
	case OpCode::GET_LOCAL_VEC2  : return +2;
	case OpCode::GET_LOCAL_VEC3  : return +3;
	case OpCode::GET_LOCAL_ARRAY : return +arg;
	case OpCode::GET_LOCAL_STRUCT: return +arg;

	case OpCode::GET_LOCAL_REF       : return +1;

	case OpCode::CREATE_STRUCT: return 0;
	case OpCode::COPY_STRUCT: return 0;

	case OpCode::GET_PROPERTY_BOOL  : return +1;
	case OpCode::GET_PROPERTY_INT   : return +1;
	case OpCode::GET_PROPERTY_DOUBLE: return +1;
	case OpCode::GET_PROPERTY_VEC2  : return +2;
	case OpCode::GET_PROPERTY_VEC3  : return +3;
	case OpCode::GET_PROPERTY_ARRAY : return +1;
	case OpCode::GET_PROPERTY_STRUCT: return +1;

	case OpCode::GET_MEMBER_REF: return 0;

	case OpCode::GET_INDEX_BOOL  : return +arg;
	case OpCode::GET_INDEX_INT   : return +arg;
	case OpCode::GET_INDEX_DOUBLE: return +arg;
	case OpCode::GET_INDEX_VEC2  : return +arg;
	case OpCode::GET_INDEX_VEC3  : return +arg;
	case OpCode::GET_INDEX_ARRAY : return +arg;
	case OpCode::GET_INDEX_STRUCT: return +arg;

	case OpCode::GET_INDEX_REF:
	case OpCode::GET_INDEX_REF_BOOL:
	case OpCode::GET_INDEX_REF_INT:
	case OpCode::GET_INDEX_REF_DOUBLE:
	case OpCode::GET_INDEX_REF_VEC2:
	case OpCode::GET_INDEX_REF_VEC3:
		return 0;

	case OpCode::GET_VEC2_X: return 0;
	case OpCode::GET_VEC2_Y: return 0;
	case OpCode::GET_VEC2_X_REF: return 0;
	case OpCode::GET_VEC2_Y_REF: return 0;
	case OpCode::GET_VEC2_SWIZZLE: return +1; // in case the swizzle returns a vec3
	case OpCode::GET_VEC3_X: return 0;
	case OpCode::GET_VEC3_Y: return 0;
	case OpCode::GET_VEC3_Z: return 0;
	case OpCode::GET_VEC3_X_REF: return 0;
	case OpCode::GET_VEC3_Y_REF: return 0;
	case OpCode::GET_VEC3_Z_REF: return 0;
	case OpCode::GET_VEC3_SWIZZLE: return +1;
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
	case OpCode::ADD_DOUBLE: return -1;
	case OpCode::SUB_DOUBLE: return -1;
	case OpCode::MUL_DOUBLE: return -1;
	case OpCode::DIV_DOUBLE: return -1;
	case OpCode::EXP_DOUBLE: return -1;
	case OpCode::GT_DOUBLE: return -1;
	case OpCode::LT_DOUBLE: return -1;
	case OpCode::GE_DOUBLE: return -1;
	case OpCode::LE_DOUBLE: return -1;
	case OpCode::NEG_VEC2: return 0;
	case OpCode::ADD_VEC2: return -2;
	case OpCode::SUB_VEC2: return -2;
	case OpCode::DOT_VEC2: return -3;
	case OpCode::MUL_VEC2_DOUBLE: return -1;
	case OpCode::MUL_DOUBLE_VEC2: return -1;
	case OpCode::NEG_VEC3: return 0;
	case OpCode::ADD_VEC3: return -3;
	case OpCode::SUB_VEC3: return -3;
	case OpCode::DOT_VEC3: return -5;
	case OpCode::MUL_VEC3_DOUBLE: return -1;
	case OpCode::MUL_DOUBLE_VEC3: return -1;
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
	case OpCode::STORE_ARRAY:
	case OpCode::STORE_STRUCT:
		return -1;

	case OpCode::POP_VOID  : return -1;
	case OpCode::POP_BOOL  : return -1;
	case OpCode::POP_INT   : return -1;
	case OpCode::POP_DOUBLE: return -1;
	case OpCode::POP_VEC2  : return -2;
	case OpCode::POP_VEC3  : return -3;
	case OpCode::POP_ARRAY : return -arg;
	case OpCode::POP_STRUCT: return -arg;

	case OpCode::CALL: return -arg + 1; // TODO: This is not correct anymore!

	case OpCode::RETURN_VOID:
	case OpCode::RETURN_BOOL: 
	case OpCode::RETURN_INT: 
	case OpCode::RETURN_DOUBLE: 
	case OpCode::RETURN_VEC2: 
	case OpCode::RETURN_VEC3: 
	case OpCode::RETURN_ARRAY: 
	case OpCode::RETURN_STRUCT: 
		return 0;
	default:
		assert(false);
		return 0;
	}
}

uint8_t Compiler::addConstant(const Value& v)
{
	prg.constants.push_back(v);
	return (uint8_t)(prg.constants.size() - 1);
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
	else if (auto r = dynamic_cast<ReturnStmt*    >(stmt)) compileReturn(r);
	else if (auto s = dynamic_cast<StructStmt*    >(stmt)) compileStruct(s);
	else
		throw std::runtime_error("Unsupported statement type");
}

void Compiler::compileExprStmt(ExpressionStmt* stmt)
{
	Type type = compileExpression(stmt->expr.get());

	switch (type->kind)
	{
	case TypeKind::Void  : emit(OpCode::POP_VOID  ); break;
	case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
	case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
	case TypeKind::Vec2  : emit(OpCode::POP_VEC2  ); break;
	case TypeKind::Vec3  : emit(OpCode::POP_VEC3  ); break;
	case TypeKind::Array : 
		emit(OpCode::POP_ARRAY, type->size());
		emitUint8(type->typeIndex);
		break;
	case TypeKind::Struct: 
		emit(OpCode::POP_STRUCT, type->size());
		emitUint8(type->typeIndex);
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

void Compiler::compileInitializer(Expression* expr, Type expectedType)
{
	if (expectedType->kind == TypeKind::Struct)
	{
		InitializerExpr* init = dynamic_cast<InitializerExpr*>(expr);
		if (init)
		{
			if (init->elements.size() != expectedType->fields.size())
				throw std::runtime_error("Struct initializer has incorrect number of fields.");

			for (size_t i = 0; i < init->elements.size(); ++i)
			{
				Type fieldType = expectedType->fields[i].first;
				compileInitializer(init->elements[i].get(), fieldType);
			}

			emit(OpCode::CREATE_STRUCT);
		}
		else
		{
			Type returnType = compileExpression(expr);
			if (returnType != expectedType)
				throw std::runtime_error("Invalid initializer type for struct.");

			emit(OpCode::COPY_STRUCT);
		}
	}
	else if (expectedType->kind == TypeKind::Array)
	{
		InitializerExpr* init = dynamic_cast<InitializerExpr*>(expr);
		if (init)
		{
			if (init->elements.size() != expectedType->arraySize)
				throw std::runtime_error("Array initializer has incorrect number of elements.");

			Type arrayType = expectedType->elementType;
			for (size_t i = 0; i < init->elements.size(); ++i)
			{
				compileInitializer(init->elements[i].get(), arrayType);
			}
		}
		else
		{
			Type returnType = compileExpression(expr);
			if (returnType != expectedType)
				throw std::runtime_error("Invalid initializer type for array.");
		}
	}
	else if (expectedType->kind == TypeKind::Vec2)
	{
		InitializerExpr* init = dynamic_cast<InitializerExpr*>(expr);
		if (init)
		{
			if (init->elements.size() != 2)
				throw std::runtime_error("Vec2 initializer must have exactly 2 elements.");
			for (size_t i = 0; i < 2; ++i)
			{
				Type elemType = compileExpression(init->elements[i].get());
				if (elemType != prg.types.getDoubleType())
					throw std::runtime_error("Vec2 initializer elements must be of type double.");
			}
		}
		else
		{
			Type returnType = compileExpression(expr);
			if (returnType != expectedType)
				throw std::runtime_error("Invalid initializer type for vec2.");
		}
	}
	else if (expectedType->kind == TypeKind::Vec3)
	{
		InitializerExpr* init = dynamic_cast<InitializerExpr*>(expr);
		if (init)
		{
			if (init->elements.size() != 3)
				throw std::runtime_error("Vec3 initializer must have exactly 3 elements.");
			for (size_t i = 0; i < 3; ++i)
			{
				Type elemType = compileExpression(init->elements[i].get());
				if (elemType != prg.types.getDoubleType())
					throw std::runtime_error("Vec3 initializer elements must be of type double.");
			}
		}
		else
		{
			Type returnType = compileExpression(expr);
			if (returnType != expectedType)
				throw std::runtime_error("Invalid initializer type for vec3.");
		}
	}
	else
	{
		Type returnType = compileExpression(expr);
		if (returnType != expectedType)
			throw std::runtime_error("Invalid initializer type");
	}
}

void Compiler::compileVarDecl(VarDeclStmt* decl)
{
	Type baseType = decl->type;
	for (auto& var : decl->vars)
	{
		Type type = baseType;
		if (var.arraySizes.size() > 0)
		{
			type = prg.types.getArrayType(baseType, var.arraySizes);
		}

		if (var.initializer)
		{
			compileInitializer(var.initializer.get(), type);
		}

		if (m_scopeDepth == 0)
		{
			prg.addGlobal(var.name, type);

			if (var.initializer)
			{
				auto it = prg.globalIndices.find(var.name);

				Program::Global& global = prg.globals[it->second];
				global.isInitialized = true;

				switch (type->kind)
				{
				case TypeKind::Bool  : emit(OpCode::SET_GLOBAL_BOOL  ); break;
				case TypeKind::Int   : emit(OpCode::SET_GLOBAL_INT   ); break;
				case TypeKind::Double: emit(OpCode::SET_GLOBAL_DOUBLE); break;
				case TypeKind::Vec2  : emit(OpCode::SET_GLOBAL_VEC2  ); break;
				case TypeKind::Vec3  : emit(OpCode::SET_GLOBAL_VEC3  ); break;
				case TypeKind::Array : 
					emit(OpCode::SET_GLOBAL_ARRAY);
					emitUint8(type->size());
					break;
				case TypeKind::Struct:
					emit(OpCode::SET_GLOBAL_STRUCT);
					emitUint8(type->size());
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
				case TypeKind::Array : 
					emit(OpCode::POP_ARRAY, type->size());
					emitUint8(type->typeIndex);
					break;
				case TypeKind::Struct:
					emit(OpCode::POP_STRUCT, type->size());
					emitUint8(type->typeIndex);
					break;
				default:
					throw std::runtime_error("Unsupported global variable type");
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
				if (m_locals[i].name == var.name)
					throw std::runtime_error("Variable '" + var.name + "' is already declared in this scope.");
			}

			m_locals.push_back({ var.name, type, m_scopeDepth, localStackSize });

			if (var.initializer)
			{
				m_locals.back().isInitialized = true;
			}

			localStackSize += (int)type->size();
		}
	}
}

//
// ===== If / While =====
//

void Compiler::compileIf(IfStmt* stmt)
{
	Type type = compileExpression(stmt->condition.get());
	if (!isNumericType(type) && !isBoolType(type))
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
	if (!isNumericType(type) && type != prg.types.getBoolType())
		throw std::runtime_error("Condition expression must be of a numeric type or bool.");

	int exitJump = emitJump(OpCode::JUMP_IF_FALSE);

	switch (type->kind)
	{
	case TypeKind::Bool  : emit(OpCode::POP_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::POP_INT   ); break;
	case TypeKind::Double: emit(OpCode::POP_DOUBLE); break;
	default:
		throw std::runtime_error("Unsupported condition type in while statement");
	}

	compileStatement(stmt->body.get());

	emitLoop(loopStart);

	patchJump(exitJump);
	emit(OpCode::POP_BOOL);
}

void Compiler::compileReturn(ReturnStmt* stmt)
{
	if (stmt->value)
	{
		Type returnType = expectedReturnType;
		auto init = dynamic_cast<InitializerExpr*>(stmt->value.get());
		if (init && expectedReturnType != nullptr)
		{
			compileInitializer(init, expectedReturnType);
		}
		else
		{
			returnType = compileExpression(stmt->value.get());

			if (expectedReturnType != nullptr && returnType != expectedReturnType)
				throw std::runtime_error("Return type mismatch: expected " + TypeToString(expectedReturnType) + ", got " + TypeToString(returnType));
		}

		switch (returnType->kind)
		{
		case TypeKind::Bool  : emit(OpCode::RETURN_BOOL  ); break;
		case TypeKind::Int   : emit(OpCode::RETURN_INT   ); break;
		case TypeKind::Double: emit(OpCode::RETURN_DOUBLE); break;
		case TypeKind::Vec2  : emit(OpCode::RETURN_VEC2  ); break;
		case TypeKind::Vec3  : emit(OpCode::RETURN_VEC3  ); break;
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
		if (expectedReturnType != nullptr && expectedReturnType != prg.types.getVoidType())
			throw std::runtime_error("Missing return value in function with non-void return type.");

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

	FunctionInfo info;
	info.name = fn->name;
	info.entry = prg.code.size();
	info.returnType = fn->returnType;

	for (auto& p : fn->params)
		info.args.push_back(p.first);

	int fnIndex = (int)prg.functions.size();

	currentFunction = fnIndex;

	Type currentReturnType = expectedReturnType;
	expectedReturnType = fn->returnType;
	bool hasReturnBefore = hasReturn;
	hasReturn = false;

	beginScope();

	for (auto& p : fn->params)
	{
		m_locals.push_back({ p.second, p.first, m_scopeDepth, localStackSize, true });
		localStackSize += (int)p.first->size();
		info.argSize += (int)p.first->size();
		stackDepth += (int)p.first->size();
	}
	prg.functions.push_back(info);

	size_t currentStackSize = stackDepth;

	BlockStmt* body = dynamic_cast<BlockStmt*>(fn->body.get());
	for (auto& stmt : body->statements)
		compileStatement(stmt.get());

	if (!hasReturn)
	{
		if (fn->returnType != nullptr && fn->returnType != prg.types.getVoidType())
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
		case TypeKind::Array : 
			emit(OpCode::POP_ARRAY, local.type->size()); 
			emitUint8(local.type->typeIndex);
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

	prg.functions[fnIndex].maxStackSize = maxStackDepth - currentStackSize;

	currentFunction = -1;
}

//
// ===== Expressions =====
//

Type Compiler::compileExpression(Expression* expr)
{
	Type type;
	if      (auto b = dynamic_cast<BinaryExpr*       >(expr)) type = compileBinary(b);
	else if (auto u = dynamic_cast<UnaryExpr*        >(expr)) type = compileUnary(u);
	else if (auto l = dynamic_cast<LiteralExpr*      >(expr)) type = compileLiteral(l);
	else if (auto v = dynamic_cast<VariableExpr*     >(expr)) type = compileVariable(v);
	else if (auto a = dynamic_cast<AssignExpr*       >(expr)) type = compileAssign(a);
	else if (auto c = dynamic_cast<CallExpr*         >(expr)) type = compileCall(c);
	else if (auto m = dynamic_cast<MemberExpr*       >(expr)) type = compileMember(m);
	else if (auto s = dynamic_cast<IndexExpr*        >(expr)) type = compileIndex(s);
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
	default:
		throw std::runtime_error("Unsupported literal type");
	}

	uint8_t idx = addConstant(expr->value);
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

		if (!m_locals[local].isInitialized)
			throw std::runtime_error("Cannot read uninitialized local variable: " + expr->name);

		return type;
	}

	int globalIndex = resolveGlobal(expr->name);
	Program::Global& glob = prg.globals[globalIndex];

	if (!glob.isInitialized)
		throw std::runtime_error("Cannot read uninitialized global variable: " + expr->name);

	glob.refcount++;

	switch (glob.type->kind)
	{
	case TypeKind::Bool  : emit(OpCode::GET_GLOBAL_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::GET_GLOBAL_INT   ); break;
	case TypeKind::Double: emit(OpCode::GET_GLOBAL_DOUBLE); break;
	case TypeKind::Vec2  : emit(OpCode::GET_GLOBAL_VEC2  ); break;
	case TypeKind::Vec3  : emit(OpCode::GET_GLOBAL_VEC3  ); break;
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
	Type r_type = compileExpression(expr->value.get());

	if (l_type != r_type)
		throw std::runtime_error("Type mismatch in assignment: cannot assign " + TypeToString(r_type) + " to " + TypeToString(l_type));

	switch (l_type->kind)
	{
	case TypeKind::Bool  : emit(OpCode::STORE_BOOL  ); break;
	case TypeKind::Int   : emit(OpCode::STORE_INT   ); break;
	case TypeKind::Double: emit(OpCode::STORE_DOUBLE); break;
	case TypeKind::Vec2  : emit(OpCode::STORE_VEC2  ); break;
	case TypeKind::Vec3  : emit(OpCode::STORE_VEC3  ); break;
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
			return prg.types.getDoubleType();
		}
		else if (expr->property == "y")
		{
			emit(OpCode::GET_VEC2_Y_REF);
			return prg.types.getDoubleType();
		}
		else
			throw std::runtime_error("Vec2 has no member named '" + expr->property + "'.");
	}
	else if (isVec3Type(objectType))
	{
		if (expr->property == "x")
		{
			emit(OpCode::GET_VEC3_X_REF);
			return prg.types.getDoubleType();
		}
		else if (expr->property == "y")
		{
			emit(OpCode::GET_VEC3_Y_REF);
			return prg.types.getDoubleType();
		}
		else if (expr->property == "z")
		{
			emit(OpCode::GET_VEC3_Z_REF);
			return prg.types.getDoubleType();
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

		if (!isNumericType(type) && !isBoolType(type))
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
		if (!isNumericType(type) && !isBoolType(type))
			throw std::runtime_error("Cannot convert right operand of '&&' to boolean.");

		patchJump(endJump);
		return prg.types.getBoolType();
	}

	// Short-circuit OR
	if (op == BinaryOp::OrOr)
	{
		Type type = compileExpression(expr->left.get());

		if (!isNumericType(type) && !isBoolType(type))
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
		if (!isNumericType(type) && !isBoolType(type))
			throw std::runtime_error("Cannot convert right operand of '||' to boolean.");

		patchJump(endJump);
		return prg.types.getBoolType();
	}

	Type type_l = compileExpression(expr->left.get());
	Type type_r = compileExpression(expr->right.get());

	if (isBoolType(type_l) && isBoolType(type_r))
	{
		switch (op)
		{
		case BinaryOp::EqualEqual: emit(OpCode::EQUAL_BOOL); type_l = prg.types.getBoolType(); break;
		case BinaryOp::NotEqual  : emit(OpCode::NEQ_BOOL  ); type_l = prg.types.getBoolType(); break;
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
		case BinaryOp::Greater     : emit(OpCode::GT_INT   ); type_l = prg.types.getBoolType(); break;
		case BinaryOp::Less        : emit(OpCode::LT_INT   ); type_l = prg.types.getBoolType(); break;
		case BinaryOp::GreaterEqual: emit(OpCode::GE_INT   ); type_l = prg.types.getBoolType(); break;
		case BinaryOp::LessEqual   : emit(OpCode::LE_INT   ); type_l = prg.types.getBoolType(); break;
		case BinaryOp::EqualEqual  : emit(OpCode::EQUAL_INT); type_l = prg.types.getBoolType(); break;
		case BinaryOp::NotEqual    : emit(OpCode::NEQ_INT  ); type_l = prg.types.getBoolType(); break;
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
		case BinaryOp::Greater     : emit(OpCode::GT_DOUBLE); type_l = prg.types.getBoolType(); break;
		case BinaryOp::Less        : emit(OpCode::LT_DOUBLE); type_l = prg.types.getBoolType(); break;
		case BinaryOp::GreaterEqual: emit(OpCode::GE_DOUBLE); type_l = prg.types.getBoolType(); break;
		case BinaryOp::LessEqual   : emit(OpCode::LE_DOUBLE); type_l = prg.types.getBoolType(); break;
		case BinaryOp::EqualEqual  : emit(OpCode::EQUAL_DOUBLE); type_l = prg.types.getBoolType(); break;
		case BinaryOp::NotEqual    : emit(OpCode::NEQ_DOUBLE  ); type_l = prg.types.getBoolType(); break;
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
		case BinaryOp::Multiply: emit(OpCode::DOT_VEC2); returnType = prg.types.getDoubleType(); break;
		default: throw std::runtime_error("Unsupported binary op for vec2 type.");
		}
		return returnType;
	}

	// vec2 * scalar
	if (isVec2Type(type_l) && isDoubleType(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_VEC2_DOUBLE); break;
		default: throw std::runtime_error("Unsupported binary op for vec2 and double types.");
		}
		return type_l;
	}

	// scalar * vec2
	if (isDoubleType(type_l) && isVec2Type(type_r))
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
		case BinaryOp::Multiply: emit(OpCode::DOT_VEC3); returnType = prg.types.getDoubleType(); break;
		default: throw std::runtime_error("Unsupported binary op for vec3 type.");
		}
		return returnType;
	}

	// vec3 * scalar
	if (isVec3Type(type_l) && isDoubleType(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_VEC3_DOUBLE); break;
		default: throw std::runtime_error("Unsupported binary op for vec3 and double types.");
		}
		return type_l;
	}

	// scalar * vec3
	if (isDoubleType(type_l) && isVec3Type(type_r))
	{
		switch (op)
		{
		case BinaryOp::Multiply: emit(OpCode::MUL_DOUBLE_VEC3); break;
		default: throw std::runtime_error("Unsupported binary op for double and vec3 types.");
		}
		return type_r;
	}

	throw std::runtime_error("Unsupported binary op for given operand types.");

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
		else
			throw std::runtime_error("Invalid operand type for unary '-'.");
		break;
	}
	case UnaryOp::Not:
	{
		if (isBoolType(type))
		{
			emit(OpCode::NOT);
			type = prg.types.getBoolType();
		}
		else
			throw std::runtime_error("Invalid operand type for unary '!'.");
		break;
	}
	default: throw std::runtime_error("Unsupported unary op");
	}
	return type;
}

int Compiler::resolveFunction(const std::string& name, std::vector<Type> args)
{
	for (int i = 0; i < (int)prg.functions.size(); ++i)
	{
		if ((prg.functions[i].name == name) && (prg.functions[i].args == args))
		{
			return i ;
		}
	}
	return -1;
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
	// Callee must be a variable (function name)
	if (auto calleeVar = dynamic_cast<VariableExpr*>(call->callee.get()))
	{
		// don't copy args for native functions
		std::vector<Type> argTypes = compileFncArgs(call->arguments);
		int fnIndex = resolveFunction(calleeVar->name, argTypes);
		if (fnIndex == -1)
			throw std::runtime_error("Undefined function: " + calleeVar->name);

		// We don't allow recursive calls.
		if (currentFunction == fnIndex)
			throw std::runtime_error("Recursive calls are not allowed: " + calleeVar->name);

		// check arg count
		if (prg.functions[fnIndex].args.size() != (int)call->arguments.size())
			throw std::runtime_error("Argument count mismatch in call to function: " + calleeVar->name);

		// check arg types
		for (int i = 0; i < argTypes.size(); ++i)
		{
			if (argTypes[i] != prg.functions[fnIndex].args[i])
				throw std::runtime_error("Argument type mismatch in call to function: " + calleeVar->name);
		}

		stackDepth += (int)prg.functions[fnIndex].maxStackSize;

		emit(OpCode::CALL, (int)prg.functions[fnIndex].args.size());
		emitUint8(fnIndex);
		emitUint8((uint8_t)call->arguments.size());

		return prg.functions[fnIndex].returnType;
	}
	else
	{
		throw std::runtime_error("Invalid function call.");
	}
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
			return prg.types.getDoubleType();
		}
		else if (expr->property == "y")
		{
			emit(OpCode::GET_VEC2_Y);
			return prg.types.getDoubleType();
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
				return (size == 2 ? prg.types.getVec2Type() : prg.types.getVec3Type());
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
			return prg.types.getDoubleType();
		}
		else if (expr->property == "y")
		{
			emit(OpCode::GET_VEC3_Y);
			return prg.types.getDoubleType();
		}
		else if (expr->property == "z")
		{
			emit(OpCode::GET_VEC3_Z);
			return prg.types.getDoubleType();
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
				return (size == 2 ? prg.types.getVec2Type() : prg.types.getVec3Type());
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
	Type exprType = compileExpression(expr->object.get());
	Type indxType = compileExpression(expr->index.get());
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
	case TypeKind::Array : emit(OpCode::GET_INDEX_ARRAY , stackEffect); break;
	case TypeKind::Struct: emit(OpCode::GET_INDEX_STRUCT, stackEffect); break;
		default:
		throw std::runtime_error("Unsupported array element type");
	}

	emitUint8(exprType->typeIndex);

	return exprType->elementType;
}
