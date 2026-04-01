#include "differentiator.h"
#include "simplifier.h"
#include <math.h>

using namespace febcode;

Type Differentiator::getDerivativeType(Type varType, TypeKind derivType)
{
	if (derivType == TypeKind::Double)
	{
		return varType;
	}
	else if (derivType == TypeKind::Vec3)
	{
		if (varType->kind == TypeKind::Double) return prg.types.getVec3Type();
		if (varType->kind == TypeKind::Vec3  ) return prg.types.getMat3Type();
	}

	throw std::runtime_error("Can't determine type of derivative.");
}

std::unique_ptr<AST> Differentiator::differentiate(const AST& ast, const deriveVar& var)
{
	auto derivativeAst = std::make_unique<AST>();

	// first, see if the source AST actually depends on the variable we're differentiating with respect to. 
	// If not, we can just return an empty AST
	dependencyFound = false;
	for (const auto& stmt : ast.root.statements)
	{
		if (febcode::dependsOn(stmt.get(), var.name))
		{
			dependencyFound = true;
			break;
		}
	}

	if (!dependencyFound)
	{
		// the derivative of a constant is zero, so we can just return an AST with a single statement that returns zero.
		switch (prg.returnType->kind)
		{
		case TypeKind::Double:
			derivativeAst->addStatement(std::make_unique<ReturnStmt>(Literal(Value(0.0))));
			return derivativeAst;
		case TypeKind::Vec2:
			derivativeAst->addStatement(std::make_unique<ReturnStmt>(Literal(Value(vec2(0., 0.)))));
			return derivativeAst;
		case TypeKind::Vec3:
			derivativeAst->addStatement(std::make_unique<ReturnStmt>(Literal(Value(vec3(0., 0., 0.)))));
			return derivativeAst;
		case TypeKind::Mat2:
			derivativeAst->addStatement(std::make_unique<ReturnStmt>(Literal(Value(mat2(0., 0., 0., 0.)))));
			return derivativeAst;
		case TypeKind::Mat3:
			derivativeAst->addStatement(std::make_unique<ReturnStmt>(Literal(Value(mat3()))));
			return derivativeAst;
		default:
			throw std::runtime_error("Don't know how to return zero for this type.");
		}
	}

	for (const auto& stmt : ast.root.statements)
	{
		differentiateStmt(derivativeAst->root, stmt.get(), var);
	}
	return derivativeAst;
}

void Differentiator::differentiateStmt(BlockStmt& ast, Statement* stmt, const deriveVar& var)
{
	if      (auto exprStmt   = dynamic_cast<ExpressionStmt*>(stmt)) diffExpressionStmt(ast, exprStmt  , var);
	else if (auto returnStmt = dynamic_cast<ReturnStmt*    >(stmt)) diffReturnStmt    (ast, returnStmt, var); 
	else if (auto structStmt = dynamic_cast<StructStmt*    >(stmt)) diffStructStmt    (ast, structStmt, var);
	else if (auto varStmt    = dynamic_cast<VarDeclStmt*   >(stmt)) diffVarDeclStmt   (ast, varStmt   , var);
	else if (auto ifStmt     = dynamic_cast<IfStmt*        >(stmt)) diffIfStmt        (ast, ifStmt    , var);
	else if (auto blockStmt  = dynamic_cast<BlockStmt*     >(stmt)) diffBlockStmt     (ast, blockStmt , var);
	else
	{
		throw std::runtime_error("Unsupported statement type for differentiation");
	}
}

void Differentiator::diffExpressionStmt(BlockStmt& ast, ExpressionStmt* stmt, const deriveVar& var)
{
	auto derivativeExpr = simplify(differentiate(stmt->expr.get(), var).get());
	ast.addStatement(std::make_unique<ExpressionStmt>(std::move(derivativeExpr)));
}

void Differentiator::diffReturnStmt(BlockStmt& ast, ReturnStmt* stmt, const deriveVar& var)
{
	auto derivativeExpr = simplify(differentiate(stmt->value.get(), var).get());
	ast.addStatement(std::make_unique<ReturnStmt>(std::move(derivativeExpr)));
}

void Differentiator::diffStructStmt(BlockStmt& ast, StructStmt* stmt, const deriveVar& var)
{
	// copy the original struct declaration to the derivative AST
	ast.addStatement(std::make_unique<StructStmt>(stmt->name, stmt->type, stmt->fields));
}

void Differentiator::diffVarDeclStmt(BlockStmt& ast, VarDeclStmt* stmt, const deriveVar& var)
{
	std::vector<Var> copyVars; // copy of the original variables for the derivative AST
	std::vector<Var> newVars; // new variables for the derivatives in the derivative AST

	Type baseType = stmt->type;

	// determine the type of the derivative variable based on the type of the original variable and the derivative variable.
	Type derivType = getDerivativeType(baseType, var.type);

	for (auto& var_i : stmt->vars)
	{
		// handle array types by getting the base type and then reconstructing the array type
		Type type = baseType;
		if (var_i.arraySizes.size() > 0)
		{
			type = prg.types.getArrayType(baseType, var_i.arraySizes);
		}

		// copy the original variable declaration to the derivative AST
		copyVars.push_back({ var_i.name, var_i.arraySizes, copy_expression(var_i.initializer.get()), var_i.input });

		// store the type of this variable in the map for later use when creating derivative variables for it.
		varTypes[var_i.name] = type;

		// No need to differentiate non-numeric types, since they don't contribute to the derivative. 
		// We can just copy them to the derivative AST without creating a derivative variable for them.
		if (type->kind == TypeKind::Bool || type->kind == TypeKind::Int) continue;

		// also, no need to differentiate input variables, since they are treated as constants with respect to differentiation.
		if (var_i.input) continue;

		// create a new variable for the derivative of this variable
		std::string derivName = "__d" + var_i.name + "_d" + var.name;

		ExprPtr init = nullptr;
		if (var_i.initializer)
		{
			// lets differentiate the initializer.
			init = simplify(differentiate(var_i.initializer.get(), var));

			// if it's zero, then this variable didn't contribute to the derivative, so we can skip creating a derivative variable for it.
			if (isZero(init))
				continue;
		}
		else
		{
			// If the variable was not initialized, it could be assigned later. 
			// To be safe, we should create a derivative variable for it and initialize it to zero.
			switch (derivType->kind)
			{
			case TypeKind::Double: init = Literal(Value(0.0)); break;
			case TypeKind::Vec2  : init = Literal(Value(vec2(0., 0.))); break;
			case TypeKind::Vec3  : init = Literal(Value(vec3(0., 0., 0.))); break;
			case TypeKind::Mat2  : init = Literal(Value(mat2(0.))); break;
			case TypeKind::Mat3  : init = Literal(Value(mat3(0.))); break;
			case TypeKind::Array : init = Initializer(type->arraySize, 0.0); break;
			case TypeKind::Struct: init = Initializer(type->fields); break;
			default:
				throw std::runtime_error("Don't know how to differentiate variable.");
			}
		}

		// add the new derivative variable to the map and the list of new variables for the derivative AST
		std::vector<size_t> arraySizes;
		if (type->kind == TypeKind::Array)
		{
			arraySizes.push_back(type->arraySize);
		}
		deriveVars[var_i.name] = derivName;
		varTypes[derivName] = derivType;
		newVars.push_back({ derivName, arraySizes, std::move(init)});
	}

	// create new variable declaration statements for the derivatives and add it to the derivative AST
	ast.addStatement(std::make_unique<VarDeclStmt>(stmt->type, copyVars));
	if (!newVars.empty()) ast.addStatement(std::make_unique<VarDeclStmt>(derivType, newVars));
}

void Differentiator::diffIfStmt(BlockStmt& ast, IfStmt* stmt, const deriveVar& var)
{
	// copy the condition
	auto newIf = std::make_unique<IfStmt>();
	newIf->condition = std::move(copy_expression(stmt->condition.get()));

	// differentiate the then branch
	std::unique_ptr<BlockStmt> thenStmt = std::make_unique<BlockStmt>();
	differentiateStmt(*thenStmt, stmt->thenBranch.get(), var);
	newIf->thenBranch = std::move(thenStmt);

	// differentiate the else branch if it exists
	if (stmt->elseBranch)
	{
		std::unique_ptr<BlockStmt> elseStmt = std::make_unique<BlockStmt>();
		differentiateStmt(*elseStmt, stmt->elseBranch.get(), var);
		newIf->elseBranch = std::move(elseStmt);
	}

	// create the new if statement with the differentiated branches
	ast.addStatement(std::move(newIf));
}

void Differentiator::diffBlockStmt(BlockStmt& ast, BlockStmt* stmt, const deriveVar& var)
{
	for (const auto& s : stmt->statements)
	{
		differentiateStmt(ast, s.get(), var);
	}
}

ExprPtr Differentiator::differentiate(const Expression* expr, const deriveVar& var)
{
	if      (auto literal  = dynamic_cast<const LiteralExpr*    >(expr)) return diffLiteral    (literal , var);
	else if (auto variable = dynamic_cast<const VariableExpr*   >(expr)) return diffVariable   (variable, var);
	else if (auto unary    = dynamic_cast<const UnaryExpr*      >(expr)) return diffUnary      (unary   , var);
	else if (auto binary   = dynamic_cast<const BinaryExpr*     >(expr)) return diffBinary     (binary  , var);
	else if (auto call     = dynamic_cast<const CallExpr*       >(expr)) return diffCall       (call    , var);
	else if (auto init     = dynamic_cast<const InitializerExpr*>(expr)) return diffInit       (init    , var);
	else if (auto ctor     = dynamic_cast<const ConstructorExpr*>(expr)) return diffConstructor(ctor, var);
	else if (auto assign   = dynamic_cast<const AssignExpr*     >(expr)) return diffAssign     (assign  , var);
	else if (auto index    = dynamic_cast<const IndexExpr*      >(expr)) return diffIndex      (index   , var);
	else if (auto member   = dynamic_cast<const MemberExpr*     >(expr)) return diffMember     (member  , var);
	else
		throw std::runtime_error("Unsupported expression type for differentiation");

	return nullptr;
}

ExprPtr Differentiator::diffLiteral(const LiteralExpr* literal, const deriveVar& var)
{
	// The derivative of a constant is zero
	return Literal(Value(0.0));
}

ExprPtr Differentiator::diffVariable(const VariableExpr* variable, const deriveVar& var)
{
	// The derivative of a variable with respect to itself is 1
	if (variable->name == var.name)
	{
		Type derivType = getDerivativeType(prg.types.getTypeFromKind(var.type), var.type);

		switch (derivType->kind)
		{
		case TypeKind::Double: return Literal(Value(1.0));
		case TypeKind::Mat2  : return Constructor(prg.types.getMat2Type(), 1.0);
		case TypeKind::Mat3  : return Constructor(prg.types.getMat3Type(), 1.0);
		default:
			throw std::runtime_error("Don't know how to make literal of derivative type.");
		}
	}

	// see if we have a derivative for this variable
	auto it = deriveVars.find(variable->name);
	if (it != deriveVars.end())
	{
		return Variable(it->second);
	}
	else
	{
		// we may get here if a derivative was not created for this variable 
		// (e.g. if it's a non-numeric type or has a literal initializer),
		// in which case we treat it as a constant and return zero
		return Literal(Value(0.0));
	}
}

ExprPtr Differentiator::diffUnary(const UnaryExpr* unary, const deriveVar& var)
{
	// For unary negation, the derivative is the negation of the derivative of the operand
	if (unary->op == UnaryOp::Negate) {
		auto operandDerivative = differentiate(unary->right.get(), var);
		return Negate(operandDerivative);
	}
	// For logical NOT, the derivative is zero (since it's not a mathematical operation)
	else if (unary->op == UnaryOp::Not) {
		return Literal(Value(0.0));
	}
	else
		throw std::runtime_error("Unsupported unary operator for differentiation");
}

std::unique_ptr<Expression> Differentiator::diffBinary(const BinaryExpr* binary, const deriveVar& var)
{
	const auto& left  = binary->left;
	const auto& right = simplify(binary->right);

	auto dleft  = differentiate(left.get(), var);
	auto dright = differentiate(right.get(), var);

	switch (binary->op)
	{
	case BinaryOp::Plus    : return dleft + dright;
	case BinaryOp::Minus   : return dleft - dright;
	case BinaryOp::Multiply: return dleft * right + left * dright;
	case BinaryOp::Divide  : return (dleft * right - left * dright) / (right * right);
	case BinaryOp::Exponent:
	{
		// some special cases first
		if (auto lit = dynamic_cast<LiteralExpr*>(right.get()))
		{
			const Value& e = lit->value;
			if (isIntNumber(e))
			{
				double p = (double)toIntNumber(e); 
				if (p == 1) return copy_expression(dleft.get());
				else if (p != 0)
				{
					return Literal(p) * Pow(left, Literal(p - 1.0)) * dleft; // d(x^p) = p * x^(p-1)
				}
			}
		}
		break;
	}
	}
	throw std::runtime_error("Unsupported binary operator for differentiation");
}

ExprPtr Differentiator::diffCall(const CallExpr* call, const deriveVar& var)
{
	if (auto calleeVar = dynamic_cast<VariableExpr*>(call->callee.get()))
	{
		std::string& fnc = calleeVar->name;
		auto& args = call->arguments;
		if (args.size() == 1)
		{
			ExprPtr dfnc;
			if      (fnc == "sin") dfnc =  Call("cos", args);
			else if (fnc == "cos") dfnc = -Call("sin", args);
			else if (fnc == "sqrt") dfnc = Literal(1.0) / (Literal(2.0) * Call("sqrt", args));
			else
				throw std::runtime_error("Don't know how to differentiate function " + fnc + ".");

			// differentiate argument
			auto diffArg = differentiate(args[0].get(), var);

			return dfnc * diffArg;
		}
	}

	throw std::runtime_error("Don't know how to compile function call.");
}

ExprPtr Differentiator::diffInit(const InitializerExpr* init, const deriveVar& var)
{
	std::vector<ExprPtr> diffElements;
	for (const auto& elem : init->elements)
	{
		diffElements.push_back(differentiate(elem.get(), var));
	}
	return std::make_unique<InitializerExpr>(std::move(diffElements));
}

ExprPtr Differentiator::diffConstructor(const ConstructorExpr* ctor, const deriveVar& var)
{
	std::vector<ExprPtr> diffArgs;
	for (const auto& arg : ctor->args)
	{
		diffArgs.push_back(simplify(differentiate(arg.get(), var)));
	}
	return std::make_unique<ConstructorExpr>(ctor->type, std::move(diffArgs));
}

std::unique_ptr<Expression> Differentiator::diffAssign(const AssignExpr* assign, const deriveVar& var)
{
	// For an assignment expression, we can use the rule: d( y = expr ) --> dy = d(expr)
	auto du = differentiate(assign->target.get(), var);
	auto dv = differentiate(assign->value.get(), var);

	// if the left is zero, that meants that it did not depend on the variable.
	// In that case, we just copy the original expression, since it doesn't contribute to the derivative.
	if (isZero(du))
		return copy_expression(assign);

	return Assign(du, dv);
}

std::unique_ptr<Expression> Differentiator::diffMember(const MemberExpr* member, const deriveVar& var)
{
	// For a member access expression, we can use the rule: d( obj.field ) --> dobj.field
	auto dobj = differentiate(member->object.get(), var);
	return Member(dobj, member->property);
}

std::unique_ptr<Expression> Differentiator::diffIndex(const IndexExpr* index, const deriveVar& var)
{
	// For an index access expression, we can use the rule: d( obj.[i] ) --> dobj[i]
	auto dobj = differentiate(index->object.get(), var);
	return Index(dobj, index->index);
}

static bool dependsOn(const Expression* expr, const std::string& var)
{
	if (auto literal = dynamic_cast<const LiteralExpr*>(expr))
	{
		return false;
	}
	else if (auto variable = dynamic_cast<const VariableExpr*>(expr))
	{
		return variable->name == var;
	}
	else if (auto binary = dynamic_cast<const BinaryExpr*>(expr))
	{
		return dependsOn(binary->left.get(), var) || dependsOn(binary->right.get(), var);
	}
	else if (auto unary = dynamic_cast<const UnaryExpr*>(expr))
	{
		return dependsOn(unary->right.get(), var);
	}
	else if (auto call = dynamic_cast<const CallExpr*>(expr))
	{
		for (const auto& arg : call->arguments)
		{
			if (dependsOn(arg.get(), var))
				return true;
		}
		return false;
	}
	else if (auto member = dynamic_cast<const MemberExpr*>(expr))
	{
		return dependsOn(member->object.get(), var);
	}
	else if (auto init = dynamic_cast<const InitializerExpr*>(expr))
	{
		for (const auto& elem : init->elements)
		{
			if (dependsOn(elem.get(), var))
				return true;
		}
		return false;
	}
	else if (auto assign = dynamic_cast<const AssignExpr*>(expr))
	{
		return dependsOn(assign->target.get(), var) || dependsOn(assign->value.get(), var);
	}
	else if (auto constructor = dynamic_cast<const ConstructorExpr*>(expr))
	{
		for (const auto& elem : constructor->args)
		{
			if (dependsOn(elem.get(), var))
				return true;
		}
		return false;
	}
	else
		throw std::runtime_error("Unsupported expression type for dependency analysis");
}

bool febcode::dependsOn(const Statement* stmt, const std::string& varName)
{
	if (auto exprStmt = dynamic_cast<const ExpressionStmt*>(stmt))
	{
		return ::dependsOn(exprStmt->expr.get(), varName);
	}
	else if (auto returnStmt = dynamic_cast<const ReturnStmt*>(stmt))
	{
		return ::dependsOn(returnStmt->value.get(), varName);
	}
	else if (auto structStmt = dynamic_cast<const StructStmt*>(stmt))
	{
		return false;
	}
	else if (auto varDeclStmt = dynamic_cast<const VarDeclStmt*>(stmt))
	{
		for (const auto& var : varDeclStmt->vars)
		{
			if (::dependsOn(var.initializer.get(), varName))
				return true;
		}
		return false;
	}
	else if (auto ifStmt = dynamic_cast<const IfStmt*>(stmt))
	{
		if (::dependsOn(ifStmt->condition.get(), varName)) return true;
		if (dependsOn(ifStmt->thenBranch.get(), varName)) return true;
		if (ifStmt->elseBranch && dependsOn(ifStmt->elseBranch.get(), varName)) return true;
		return false;
	}
	else if (auto block = dynamic_cast<const BlockStmt*>(stmt))
	{
		for (const StmtPtr& s : block->statements)
		{
			if (dependsOn(s.get(), varName))
				return true;
		}
		return false;
	}
	else
		throw std::runtime_error("Unsupported statement type for dependency analysis");
}
