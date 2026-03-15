#include "differentiator.h"
#include "simplifier.h"
#include <math.h>

using namespace febcode;

std::unique_ptr<AST> Differentiator::differentiate(const AST& ast, const std::string& var)
{
	auto derivativeAst = std::make_unique<AST>();
	for (const auto& stmt : ast.statements)
	{
		differentiateStmt(*derivativeAst, stmt.get(), var);
	}
	return derivativeAst;
}

void Differentiator::differentiateStmt(AST& ast, Statement* stmt, const std::string& var)
{
	if      (auto exprStmt   = dynamic_cast<ExpressionStmt*>(stmt)) diffExpressionStmt(ast, exprStmt  , var);
	else if (auto returnStmt = dynamic_cast<ReturnStmt*    >(stmt)) diffReturnStmt    (ast, returnStmt, var); 
	else if (auto structStmt = dynamic_cast<StructStmt*    >(stmt)) diffStructStmt    (ast, structStmt, var);
	else if (auto varStmt    = dynamic_cast<VarDeclStmt*   >(stmt)) diffVarDeclStmt   (ast, varStmt   , var);
	else
	{
		throw std::runtime_error("Unsupported statement type for differentiation");
	}
}

void Differentiator::diffExpressionStmt(AST& ast, ExpressionStmt* stmt, const std::string& var)
{
	auto derivativeExpr = simplify(differentiate(stmt->expr.get(), var).get());
	ast.addStatement(std::make_unique<ExpressionStmt>(std::move(derivativeExpr)));
}

void Differentiator::diffReturnStmt(AST& ast, ReturnStmt* stmt, const std::string& var)
{
	auto derivativeExpr = simplify(differentiate(stmt->value.get(), var).get());
	ast.addStatement(std::make_unique<ReturnStmt>(std::move(derivativeExpr)));
}

void Differentiator::diffStructStmt(AST& ast, StructStmt* stmt, const std::string& var)
{
	// copy the original struct declaration to the derivative AST
	ast.addStatement(std::make_unique<StructStmt>(stmt->name, stmt->type, stmt->fields));
}

void Differentiator::diffVarDeclStmt(AST& ast, VarDeclStmt* stmt, const std::string& var)
{
	std::vector<Var> copyVars; // copy of the original variables for the derivative AST
	std::vector<Var> newVars; // new variables for the derivatives in the derivative AST

	Type baseType = stmt->type;
	for (auto& var_i : stmt->vars)
	{
		// handle array types by getting the base type and then reconstructing the array type
		Type type = baseType;
		if (var_i.arraySizes.size() > 0)
		{
			type = prg.types.getArrayType(baseType, var_i.arraySizes);
		}

		// copy the original variable declaration to the derivative AST
		copyVars.push_back({ var_i.name, var_i.arraySizes, copy_expression(var_i.initializer.get()) });

		// create a new variable for the derivative of this variable
		std::string derivName = "__d" + var_i.name + "_d" + var;

		// only differentiate numeric types
		if (type->kind == TypeKind::Double)
		{
			deriveVars[var_i.name] = derivName;
			ExprPtr init = (var_i.initializer ? simplify(differentiate(var_i.initializer.get(), var).get()) : Literal(Value(0.0)));
			newVars.push_back({ derivName, std::vector<size_t>(), std::move(init) });
		}
		else if (type->kind == TypeKind::Vec2)
		{
			deriveVars[var_i.name] = derivName;
			ExprPtr init = (var_i.initializer ? simplify(differentiate(var_i.initializer.get(), var).get()) : Literal(Value(vec2(0., 0.))));
			newVars.push_back({ derivName, std::vector<size_t>(), std::move(init) });
		}
		else if (type->kind == TypeKind::Vec3)
		{
			deriveVars[var_i.name] = derivName;
			ExprPtr init = (var_i.initializer ? simplify(differentiate(var_i.initializer.get(), var).get()) : Literal(Value(vec3(0., 0., 0.))));
			newVars.push_back({ derivName, std::vector<size_t>(), std::move(init) });
		}
		else if (type->kind == TypeKind::Bool || type->kind == TypeKind::Int || type->kind == TypeKind::String)
		{
			// don't need to derive these types
		}
		else if (type->kind == TypeKind::Array)
		{
			deriveVars[var_i.name] = derivName;

			ExprPtr init = (var_i.initializer ? simplify(differentiate(var_i.initializer.get(), var).get()) : Initializer(type->arraySize, 0.0));
			newVars.push_back({ derivName, {type->arraySize}, std::move(init)});
		}
		else if (type->kind == TypeKind::Struct)
		{
			deriveVars[var_i.name] = derivName;
			ExprPtr init = (var_i.initializer ? simplify(differentiate(var_i.initializer.get(), var).get()) : Initializer(type->fields));
			newVars.push_back({ derivName, std::vector<size_t>(), std::move(init)});
		}
		else
			throw std::runtime_error("Don't know how to differentiate variable.");
	}

	// create new variable declaration statements for the derivatives and add it to the derivative AST
	ast.addStatement(std::make_unique<VarDeclStmt>(stmt->type, copyVars));
	ast.addStatement(std::make_unique<VarDeclStmt>(stmt->type, newVars ));
}

ExprPtr Differentiator::differentiate(const Expression* expr, const std::string& var) 
{
	if      (auto literal  = dynamic_cast<const LiteralExpr*    >(expr)) return diffLiteral (literal , var);
	else if (auto variable = dynamic_cast<const VariableExpr*   >(expr)) return diffVariable(variable, var);
	else if (auto unary    = dynamic_cast<const UnaryExpr*      >(expr)) return diffUnary   (unary   , var);
	else if (auto binary   = dynamic_cast<const BinaryExpr*     >(expr)) return diffBinary  (binary  , var);
	else if (auto call     = dynamic_cast<const CallExpr*       >(expr)) return diffCall    (call    , var);
	else if (auto init     = dynamic_cast<const InitializerExpr*>(expr)) return diffInit    (init    , var);
	else if (auto assign   = dynamic_cast<const AssignExpr*     >(expr)) return diffAssign  (assign  , var);
	else if (auto member   = dynamic_cast<const MemberExpr*     >(expr)) return diffMember  (member  , var);
	else
		throw std::runtime_error("Unsupported expression type for differentiation");

	return nullptr;
}

ExprPtr Differentiator::diffLiteral(const LiteralExpr* literal, const std::string& var)
{
	// The derivative of a constant is zero
	return Literal(Value(0.0));
}

ExprPtr Differentiator::diffVariable(const VariableExpr* variable, const std::string& var)
{
	// The derivative of a variable with respect to itself is 1
	if (variable->name == var) return Literal(Value(1.0));

	// see if we have a derivative for this variable
	auto it = deriveVars.find(variable->name);
	if (it != deriveVars.end())
	{
		return Variable(it->second);
	}
	else
		return Literal(Value(0.0));
}

ExprPtr Differentiator::diffUnary(const UnaryExpr* unary, const std::string& var)
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

std::unique_ptr<Expression> Differentiator::diffBinary(const BinaryExpr* binary, const std::string& var)
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

ExprPtr Differentiator::diffCall(const CallExpr* call, const std::string& var)
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
			else
				throw std::runtime_error("Don't know how to differentiate function.");

			// differentiate argument
			auto diffArg = differentiate(args[0].get(), var);

			return dfnc * diffArg;
		}
	}

	throw std::runtime_error("Don't know how to compile function call.");
}

ExprPtr Differentiator::diffInit(const InitializerExpr* init, const std::string& var)
{
	std::vector<ExprPtr> diffElements;
	for (const auto& elem : init->elements)
	{
		diffElements.push_back(differentiate(elem.get(), var));
	}
	return std::make_unique<InitializerExpr>(std::move(diffElements));
}

std::unique_ptr<Expression> Differentiator::diffAssign(const AssignExpr* assign, const std::string& var)
{
	// For an assignment expression, we can use the rule: d( y = expr ) --> dy = d(expr)
	auto du = differentiate(assign->target.get(), var);
	auto dv = differentiate(assign->value.get(), var);
	return Assign(du, dv);
}

std::unique_ptr<Expression> Differentiator::diffMember(const MemberExpr* member, const std::string& var)
{
	// For a member access expression, we can use the rule: d( obj.field ) --> dobj.field
	auto dobj = differentiate(member->object.get(), var);
	return Member(dobj, member->property);
}
