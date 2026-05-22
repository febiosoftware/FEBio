#include "resolver.h"

using namespace febcode;

Resolver::Resolver(Program& prg) : prg(prg)
{

}

Resolver::~Resolver()
{

}

void Resolver::resolve()
{
	// add all global variables to the variables map so they can be resolved in expressions
	pushScope();
	for (const auto& var : prg.globalIndices)
	{
		auto& global = prg.globals[var.second];
		declare(var.first, global.type);
	}

	// resolve all statements in the program
	resolveBlockStmt(&prg.ast->root);
}

void Resolver::resolveStatement(Statement* stmt)
{
	if (stmt == nullptr)
		throw std::runtime_error("Null statement cannot be resolved.");

	if      (auto e = dynamic_cast<ExpressionStmt*>(stmt)) resolveExprStmt    (e);
	else if (auto v = dynamic_cast<VarDeclStmt*   >(stmt)) resolveVarDecl     (v);
	else if (auto i = dynamic_cast<IfStmt*        >(stmt)) resolveIfStmt      (i);
	else if (auto w = dynamic_cast<WhileStmt*     >(stmt)) resolveWhileStmt   (w);
	else if (auto f = dynamic_cast<ForStmt*       >(stmt)) resolveForStmt     (f);
	else if (auto b = dynamic_cast<BlockStmt*     >(stmt)) resolveBlockStmt   (b);
	else if (auto r = dynamic_cast<ReturnStmt*    >(stmt)) resolveReturnStmt  (r);
	else if (auto c = dynamic_cast<FunctionStmt*  >(stmt)) resolveFunctionStmt(c);
	else if (auto s = dynamic_cast<StructStmt*    >(stmt)) resolveStructStmt  (s);
	else 
		error(stmt, "cannot resolve statement");
}

void Resolver::resolveBlockStmt(BlockStmt* stmt)
{
	pushScope();
	for (auto& s : stmt->statements)
	{
		resolveStatement(s.get());
	}
	popScope();
}

void Resolver::resolveExprStmt(ExpressionStmt* stmt)
{
	resolveExpression(stmt->expr.get());
}

void Resolver::resolveVarDecl(VarDeclStmt* stmt)
{
	for (auto& var : stmt->vars)
	{
		// resolve initializer first!
		if (var->initializer)
		{
			resolveExpression(var->initializer.get());

			// make sure initializer type can be converted to variable type
			if (commonType(var->initializer->valType, var->type) == nullptr)
				error(var->initializer.get(), "Type mismatch in variable initializer for variable '" + var->name + "'. Expected type: " + TypeToString(var->type) + ", but got: " + TypeToString(var->initializer->valType));
		}

		// now we can declare the variable in the current scope
		declare(var->name, var->type);
	}
}

void Resolver::resolveIfStmt(IfStmt* stmt)
{
	resolveExpression(stmt->condition.get());
	if (!isLogicalType(stmt->condition->valType))
		error(stmt->condition->location, "Condition expression in if statement must be of type bool. Got: " + TypeToString(stmt->condition->valType));

	resolveStatement(stmt->thenBranch.get());
	if (stmt->elseBranch)
		resolveStatement(stmt->elseBranch.get());
}

void Resolver::resolveWhileStmt(WhileStmt* stmt)
{
	resolveExpression(stmt->condition.get());
	if (!isLogicalType(stmt->condition->valType))
		error(stmt->condition->location, "Condition expression in while statement must be of type bool. Got: " + TypeToString(stmt->condition->valType));
	resolveStatement(stmt->body.get());
}

void Resolver::resolveForStmt(ForStmt* stmt)
{
	if (stmt->initializer)
		resolveStatement(stmt->initializer.get());
	resolveExpression(stmt->condition.get());
	if (stmt->condition->valType != prg.types.Bool())
		error(stmt->condition->location, "Condition expression in for statement must be of type bool. Got: " + TypeToString(stmt->condition->valType));
	resolveExpression(stmt->increment.get());
	resolveStatement(stmt->body.get());
}

void Resolver::resolveReturnStmt(ReturnStmt* stmt)
{
	if (stmt->value)
	{
		resolveExpression(stmt->value.get());
		if (prg.returnType && (stmt->value->valType != prg.returnType))
			error(stmt->value->location, "Type mismatch in return statement. Expected type: " + TypeToString(prg.returnType) + ", but got: " + TypeToString(stmt->value->valType));
	}
	else
	{
		if (prg.returnType && (prg.returnType != prg.types.Void()))
			error(stmt->location, "Missing return value in function with non-void return type.");
	}
}

void Resolver::resolveFunctionStmt(FunctionStmt* fn)
{
	FunctionInfo info;
	info.name = fn->name;
	info.entry = prg.code.size();
	info.returnType = fn->returnType;

	pushScope();
	// add function parameters to variables map
	for (const auto& param : fn->params)
	{
		info.args.push_back(param.first);
		declare(param.second, param.first);
	}
	resolveStatement(fn->body.get());
	popScope();

	// make sure the function is not defined yet
	int index = prg.resolveFunction(fn->name, info.args);
	if (index != -1)
		error(fn, "Function '" + fn->name + "' with the same parameter types is already defined.");

	// add it to the program's function list
	prg.functions.push_back(info);
}

void Resolver::resolveStructStmt(StructStmt* stmt)
{
	// nothing to do here since struct definitions are handled during parsing and type registration.
}

void Resolver::resolveExpression(Expression* expr)
{
	if (expr == nullptr)
		throw std::runtime_error("Null expression cannot be resolved.");

	if      (auto l = dynamic_cast<LiteralExpr*    >(expr)) resolveLiteral    (l);
	else if (auto v = dynamic_cast<VariableExpr*   >(expr)) resolveVariable   (v);
	else if (auto a = dynamic_cast<AssignExpr*     >(expr)) resolveAssignment (a);
	else if (auto u = dynamic_cast<UnaryExpr*      >(expr)) resolveUnary      (u);
	else if (auto b = dynamic_cast<BinaryExpr*     >(expr)) resolveBinary     (b);
	else if (auto m = dynamic_cast<MemberExpr*     >(expr)) resolveMember     (m);
	else if (auto i = dynamic_cast<IndexExpr*      >(expr)) resolveIndex      (i);
	else if (auto i = dynamic_cast<InitExpr*       >(expr)) resolveInitializer(i);
	else if (auto c = dynamic_cast<CallExpr*       >(expr)) resolveCall       (c);
	else if (auto c = dynamic_cast<ConstructorExpr*>(expr)) resolveConstructor(c);
	else
		error(expr, "cannot resolve expression");
}

void Resolver::resolveLiteral(LiteralExpr* expr)
{
	expr->valType = prg.types.getBuiltinType(expr->value);
}

void Resolver::resolveVariable(VariableExpr* expr)
{
	Type type = lookup(expr->name);
	if (type == nullptr)
		error(expr, "Undefined variable: " + expr->name);
	expr->valType = type;
}

void Resolver::resolveAssignment(AssignExpr* expr)
{
	resolveExpression(expr->target.get());
	resolveExpression(expr->value.get());
	if (expr->target->valType == nullptr)
		error(expr->target.get(), "Invalid assignment target.");

	Type type = coerce(expr->value->valType, expr->target->valType); // try to coerce value type to target type
	if (expr->target->valType != type)
		error(expr->value.get(), "Type mismatch in assignment. Expected type: " + TypeToString(expr->target->valType) + ", but got: " + TypeToString(expr->value->valType));
	expr->valType = expr->target->valType;
}

void Resolver::resolveUnary(UnaryExpr* expr)
{
	resolveExpression(expr->right.get());
	expr->valType = expr->right->valType;

	switch (expr->op)
	{
	case UnaryOp::Negate:
		if (!isNumericType(expr->valType))
			error(expr, "Unary '-' operator can only be applied to numeric types. Got: " + TypeToString(expr->valType));
		break;
	case UnaryOp::Not:
		if (expr->valType != prg.types.Bool())
			error(expr, "Unary '!' operator can only be applied to bool type. Got: " + TypeToString(expr->valType));
		break;
	default:
		error(expr, "Unsupported unary operator.");
	}
}

void Resolver::resolveBinary(BinaryExpr* expr)
{
	resolveExpression(expr->left.get());
	resolveExpression(expr->right.get());

	BinaryOpSignature sig = prg.resolveBinaryOp(expr->op, expr->left->valType, expr->right->valType);
	if (sig.resultType == nullptr)
	{
		error(expr, "Invalid binary operator \'" + opToString(expr->op) + "\' for types " + TypeToString(expr->left->valType) + " and " + TypeToString(expr->right->valType));
	}

	expr->valType = sig.resultType;
}

void Resolver::resolveMember(MemberExpr* expr)
{
	resolveExpression(expr->object.get());
	Type objType = expr->object->valType;
	if (objType == nullptr)
		error(expr->object.get(), "Invalid member access on unresolved type.");

	if (objType->kind == TypeKind::Vec2)
	{
		int l = (int)expr->property.length();
		// make sure it's a valid swizzle (only 'x', 'y' allowed, and length must be 1, 2, 3)
		for (char c : expr->property)
		{
			if ((c != 'x') && (c != 'y'))
				error(expr, "Invalid character '" + std::string(1, c) + "' in swizzle. Only 'x' and 'y' are allowed.");
		}

		if      (l == 1) expr->valType = prg.types.Double();
		else if (l == 2) expr->valType = prg.types.Vec2();
		else if (l == 3) expr->valType = prg.types.Vec3();
		else
			error(expr, "Invalid property '" + expr->property + "' for type vec2.");

		return;
	}
	else if (objType->kind == TypeKind::Vec3)
	{
		int l = (int)expr->property.length();
		// make sure it's a valid swizzle (only 'x', 'y', 'z' allowed, and length must be 1, 2, or 3)
		for (char c : expr->property)
		{
			if ((c != 'x') && (c != 'y') && (c != 'z'))
				error(expr, "Invalid character '" + std::string(1, c) + "' in swizzle. Only 'x', 'y', and 'z' are allowed.");
		}

		if      (l == 1) expr->valType = prg.types.Double();
		else if (l == 2) expr->valType = prg.types.Vec2();
		else if (l == 3) expr->valType = prg.types.Vec3();
		else
			error(expr, "Invalid property '" + expr->property + "' for type vec3.");

		return;
	}
	else if (objType->kind != TypeKind::Struct)
		error(expr, "Member access on non-struct type.");

	const auto& fields = objType->fields;
	for (const auto& field : fields)
	{
		if (field.second == expr->property)
		{
			expr->valType = field.first;
			return;
		}
	}

	error(expr, "Field '" + expr->property + "' not found in struct.");
}

void Resolver::resolveIndex(IndexExpr* expr)
{
	resolveExpression(expr->object.get());
	resolveExpression(expr->index.get());

	if (expr->index->valType != prg.types.Int())
		error(expr->index.get(), "Index must be of type int. Got: " + TypeToString(expr->index->valType));

	Type arrayType = expr->object->valType;
	if (arrayType == nullptr)
		error(expr->object.get(), "Invalid indexing on unresolved type.");
	if (arrayType->kind == TypeKind::Array)
	{
		expr->valType = arrayType->elementType;
	}
	else if (arrayType->kind == TypeKind::Mat2)
	{
		expr->valType = prg.types.Vec2();
	}
	else if (arrayType->kind == TypeKind::Mat3)
	{
		expr->valType = prg.types.Vec3();
	}
	else if (arrayType->kind == TypeKind::Vec2)
	{
		expr->valType = prg.types.Double();
	}
	else if (arrayType->kind == TypeKind::Vec3)
	{
		expr->valType = prg.types.Double();
	}
	else
	{
		error(expr, "Indexing can only be applied to array, vector, or matrix types. Got: " + TypeToString(arrayType));
	}
}

void Resolver::resolveInitializer(InitExpr* expr)
{
	Type elemType = nullptr;
	for (auto& element : expr->elements)
	{
		resolveExpression(element.get());
		if (element->valType == nullptr)
			error(element.get(), "Invalid initializer element with unresolved type.");

		if (elemType == nullptr)
			{
			elemType = element->valType;
		}
		else if (element->valType != elemType)
		{
			error(element.get(), "All elements in an initializer must have the same type. Got: " + TypeToString(elemType) + " and " + TypeToString(element->valType));
		}
	}

	expr->valType = prg.types.getArrayType(elemType, expr->elements.size());
}

void Resolver::resolveConstructor(ConstructorExpr* expr)
{
	if (expr->valType == nullptr)
		error(expr, "Invalid constructor with unresolved type.");

	for (auto& arg : expr->args)
	{
		resolveExpression(arg.get());
		if (arg->valType == nullptr)
			error(expr, "Invalid constructor argument with unresolved type.");
	}
}

void Resolver::resolveCall(CallExpr* expr)
{
	std::vector<Type> argTypes;
	for (size_t i = 0; i < expr->arguments.size(); ++i)
	{
		resolveExpression(expr->arguments[i].get());
		if (expr->arguments[i]->valType == nullptr)
			error(expr->arguments[i].get(), "Invalid function call argument with unresolved type.");
		argTypes.push_back(expr->arguments[i]->valType);
	}

	int index = prg.resolveFunction(expr->name, argTypes);
	if (index < 0)
		error(expr, "Undefined function: " + expr->name + " (or can't match arguments)");

	FunctionInfo& func = prg.functions[index];
	if (func.returnType == nullptr)
		error(expr, "Invalid function with unresolved return type.");

	expr->valType = func.returnType;
}
