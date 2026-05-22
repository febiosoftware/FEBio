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
	else if (derivType == TypeKind::Vec2)
	{
		if (varType->kind == TypeKind::Double) return prg.types.Vec2();
		if (varType->kind == TypeKind::Vec2) return prg.types.Mat2();
	}
	else if (derivType == TypeKind::Vec3)
	{
		if (varType->kind == TypeKind::Double) return prg.types.Vec3();
		if (varType->kind == TypeKind::Vec3  ) return prg.types.Mat3();
	}
	throw std::runtime_error("Can't determine type of derivative.");
}

void Differentiator::differentiate(const std::string& var)
{
	// get the variable's type
	Type varType = prg.globalType(var);
	if (varType == nullptr)
		throw std::runtime_error("Variable not found in program: " + var);

	// update the program's return type
	if (prg.returnType)
	{
		Type derivType = getDerivativeType(prg.returnType, varType->kind);
		prg.returnType = derivType;
	}

	DerivVar dvar{ var, varType };

	// differentiate the program's AST
	auto diffAST = differentiate(*prg.ast, dvar);
	prg.ast.reset(diffAST.release()); // replace original AST with derivative AST
}

std::unique_ptr<AST> Differentiator::differentiate(const AST& ast, const DerivVar& var)
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
		ExprPtr returnValue;
		if (prg.returnType)
		{
			returnValue = Zero(prg.returnType);
		}
		derivativeAst->addStatement(std::make_unique<ReturnStmt>(std::move(returnValue)));
		return derivativeAst;
	}

	// add all injected globals to vartypes
	for (const auto& global : prg.globalIndices)
	{
		std::string name = global.first;
		Type type = prg.globals[global.second].type;
		varTypes[name] = type;
	}

	for (const auto& stmt : ast.root.statements)
	{
		differentiateStmt(derivativeAst->root, stmt.get(), var);
	}
	return derivativeAst;
}

void Differentiator::differentiateStmt(BlockStmt& ast, Statement* stmt, const DerivVar& var)
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

void Differentiator::diffExpressionStmt(BlockStmt& ast, ExpressionStmt* stmt, const DerivVar& var)
{
	auto derivativeExpr = differentiate(stmt->expr.get(), var);
	ast.addStatement(std::make_unique<ExpressionStmt>(std::move(derivativeExpr)));
}

void Differentiator::diffReturnStmt(BlockStmt& ast, ReturnStmt* stmt, const DerivVar& var)
{
	auto derivativeExpr = differentiate(stmt->value.get(), var);
	ast.addStatement(std::make_unique<ReturnStmt>(std::move(derivativeExpr)));
}

void Differentiator::diffStructStmt(BlockStmt& ast, StructStmt* stmt, const DerivVar& var)
{
	// copy the original struct declaration to the derivative AST
	ast.addStatement(std::make_unique<StructStmt>(stmt->name, stmt->type, stmt->fields));
}

void Differentiator::diffVarDeclStmt(BlockStmt& ast, VarDeclStmt* stmt, const DerivVar& var)
{
	if (stmt->input)
	{
		// input variables are treated as constants with respect to differentiation, so we can just copy the variable declaration to the derivative AST without creating a derivative variable for it.
		ast.addStatement(std::make_unique<VarDeclStmt>(stmt->type, stmt->vars, stmt->input));
		return;
	}

	std::vector<VarPtr> copyVars; // copy of the original variables for the derivative AST
	std::vector<VarPtr> newVars; // new variables for the derivatives in the derivative AST

	Type baseType = stmt->type;

	// determine the type of the derivative variable based on the type of the original variable and the derivative variable.
	Type derivType = getDerivativeType(baseType, var.type->kind);

	for (auto& var_i : stmt->vars)
	{
		// handle array types by getting the base type and then reconstructing the array type
		Type type = baseType;
		if (var_i->arraySizes.size() > 0)
		{
			type = prg.types.getArrayType(baseType, var_i->arraySizes);
		}

		// copy the original variable declaration to the derivative AST
		copyVars.push_back(std::make_unique<Var>(Var{ var_i->name, var_i->arraySizes, clone(var_i->initializer.get())}));
		// store the type of this variable in the map for later use when creating derivative variables for it.
		varTypes[var_i->name] = type;

		// No need to differentiate non-numeric types, since they don't contribute to the derivative. 
		// We can just copy them to the derivative AST without creating a derivative variable for them.
		if (type->kind == TypeKind::Bool || type->kind == TypeKind::Int) continue;

		// create a new variable for the derivative of this variable
		std::string derivName = "__d" + var_i->name + "_d" + var.name;

		ExprPtr init = nullptr;
		if (var_i->initializer)
		{
			// lets differentiate the initializer.
			init = differentiate(var_i->initializer.get(), var);

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
			case TypeKind::Double: init = Literal(0.0); break;
			case TypeKind::Vec2  : init = Literal(vec2()); break;
			case TypeKind::Vec3  : init = Literal(vec3()); break;
			case TypeKind::Mat2  : init = Literal(mat2()); break;
			case TypeKind::Mat3  : init = Literal(mat3()); break;
			case TypeKind::Array : init = Initializer(type); break;
			case TypeKind::Struct: init = Constructor(type); break;
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
		deriveVars[var_i->name] = derivName;
		varTypes[derivName] = derivType;
		newVars.push_back(std::make_unique<Var>(Var{ derivName, arraySizes, std::move(init)}));
	}

	// create new variable declaration statements for the derivatives and add it to the derivative AST
	ast.addStatement(std::make_unique<VarDeclStmt>(stmt->type, copyVars));
	if (!newVars.empty()) ast.addStatement(std::make_unique<VarDeclStmt>(derivType, newVars));
}

void Differentiator::diffIfStmt(BlockStmt& ast, IfStmt* stmt, const DerivVar& var)
{
	// copy the condition
	auto newIf = std::make_unique<IfStmt>();
	newIf->condition = std::move(clone(stmt->condition.get()));

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

void Differentiator::diffBlockStmt(BlockStmt& ast, BlockStmt* stmt, const DerivVar& var)
{
	for (const auto& s : stmt->statements)
	{
		differentiateStmt(ast, s.get(), var);
	}
}

ExprPtr Differentiator::differentiate(const Expression* expr, const DerivVar& var)
{
	if      (auto literal  = dynamic_cast<const LiteralExpr*    >(expr)) return simplify(diffLiteral    (literal , var));
	else if (auto variable = dynamic_cast<const VariableExpr*   >(expr)) return simplify(diffVariable   (variable, var));
	else if (auto unary    = dynamic_cast<const UnaryExpr*      >(expr)) return simplify(diffUnary      (unary   , var));
	else if (auto binary   = dynamic_cast<const BinaryExpr*     >(expr)) return simplify(diffBinary     (binary  , var));
	else if (auto call     = dynamic_cast<const CallExpr*       >(expr)) return simplify(diffCall       (call    , var));
	else if (auto init     = dynamic_cast<const InitExpr*       >(expr)) return simplify(diffInit       (init    , var));
	else if (auto ctor     = dynamic_cast<const ConstructorExpr*>(expr)) return simplify(diffConstructor(ctor    , var));
	else if (auto assign   = dynamic_cast<const AssignExpr*     >(expr)) return simplify(diffAssign     (assign  , var));
	else if (auto index    = dynamic_cast<const IndexExpr*      >(expr)) return simplify(diffIndex      (index   , var));
	else if (auto member   = dynamic_cast<const MemberExpr*     >(expr)) return simplify(diffMember     (member  , var));
	else
		throw std::runtime_error("Unsupported expression type for differentiation");

	return nullptr;
}

ExprPtr Differentiator::diffLiteral(const LiteralExpr* literal, const DerivVar& var)
{
	// scalar derivation
	if (var.type->kind == TypeKind::Double)
	{
		// The derivative of a constant is zero
		switch (literal->value.index)
		{
		case ValueIndex::INT: return Literal(0);
		case ValueIndex::DOUBLE: return Literal(0.0);
		case ValueIndex::VEC2: return Literal(vec2());
		case ValueIndex::VEC3: return Literal(vec3());
		case ValueIndex::MAT2: return Literal(mat2());
		case ValueIndex::MAT3: return Literal(mat3());
		}
	}

	// 2d gradient
	if (var.type->kind == TypeKind::Vec2)
	{
		switch (literal->value.index)
		{
		case ValueIndex::INT:
		case ValueIndex::DOUBLE: return Literal(vec2());
		case ValueIndex::VEC2  : return Literal(mat2());
		}
	}

	// 3d gradient
	if (var.type->kind == TypeKind::Vec3)
	{
		switch (literal->value.index)
		{
		case ValueIndex::INT:
		case ValueIndex::DOUBLE: return Literal(vec3());
		case ValueIndex::VEC3  : return Literal(mat3());
		}
	}

	throw std::runtime_error("Don't know how to differentiate this literal type.");
}

ExprPtr Differentiator::diffVariable(const VariableExpr* variable, const DerivVar& var)
{
	Type derivType = getDerivativeType(variable->valType, var.type->kind);

	// The derivative of a variable with respect to itself is 1
	if (variable->name == var.name)
	{
		switch (derivType->kind)
		{
		case TypeKind::Double: return Literal(1.0);
		case TypeKind::Mat2  : return Literal(mat2(1.0));
		case TypeKind::Mat3  : return Literal(mat3(1.0));
		default:
			throw std::runtime_error("Don't know how to make literal of derivative type.");
		}
	}

	// see if we have a derivative for this variable
	auto it = deriveVars.find(variable->name);
	if (it != deriveVars.end())
	{
		return Variable(it->second, derivType);
	}
	else
	{
		// we may get here if a derivative was not created for this variable 
		// (e.g. if it's a non-numeric type or has a literal initializer),
		// in which case we treat it as a constant and return zero
		switch (derivType->kind)
		{
		case TypeKind::Double: return Literal(0.0);
		case TypeKind::Vec2: return Literal(vec2());
		case TypeKind::Vec3: return Literal(vec3());
		case TypeKind::Mat2: return Literal(mat2());
		case TypeKind::Mat3: return Literal(mat3());
		default:
			throw std::runtime_error("Don't know how to make literal of derivative type.");
		}
	}
}

ExprPtr Differentiator::diffUnary(const UnaryExpr* unary, const DerivVar& var)
{
	auto right = simplify(unary->right);

	// For unary negation, the derivative is the negation of the derivative of the operand
	// d(-f) = -df
	if (unary->op == UnaryOp::Negate) {
		auto operandDerivative = differentiate(right.get(), var);
		return Negate(operandDerivative);
	}

	throw std::runtime_error("Unsupported unary operator for differentiation");
}

std::unique_ptr<Expression> Differentiator::diffBinary(const BinaryExpr* binary, const DerivVar& var)
{
	auto left  = simplify(binary->left);
	auto right = simplify(binary->right);

	auto dleft = differentiate(left.get(), var);
	auto dright = differentiate(right.get(), var);

	// plus and minus are easy, so let's handle them first.
	if (binary->op == BinaryOp::Plus)
		return Add(dleft, dright);
	if (binary->op == BinaryOp::Minus)
		return Sub(dleft, dright);

	// The other operators are more complex. 
	// First, figure out the types of all the expressions involved, so we can determine how to apply the differentiation rules.
	Type ltype  = left->valType;
	Type rtype  = right->valType;
	Type dltype = dleft->valType;
	Type drtype = dright->valType;

	// scalar differentation
	if (var.type->kind == TypeKind::Double)
	{
		if (isScalarType(ltype) && isScalarType(rtype))
		{
			switch (binary->op)
			{
			case BinaryOp::Multiply: return Add(Mul(dleft, right), Mul(left, dright)); // d( f * g ) = df * g + f * dg
			case BinaryOp::Divide: return Div(Sub(Mul(dleft, right), Mul(left, dright)), Mul(right, right)); // d( f / g ) = (df * g - f * dg) / (g * g)
			case BinaryOp::Exponent:
			{
				// some special cases first
				if (auto lit = dynamic_cast<LiteralExpr*>(right.get()))
				{
					const Value& e = lit->value;
					if (isIntNumber(e))
					{
						double p = (double)toIntNumber(e);
						if (p == 1) return clone(dleft.get());
						else if (p != 0)
						{
							return Mul(Mul(Literal(p), Pow(left, Literal(p - 1.0))), dleft); // d(x^p) = p * x^(p-1)*dx
						}
					}
				}
				break;
			}
			}
		}
		else
		{
			switch (binary->op)
			{
			// The rule d( f * g ) = df * g + f * dg applies to all types, since it doesn't matter what the types of f and g are, as long as we can multiply them together.
			case BinaryOp::Multiply: return Add(Mul(dleft, right), Mul(left, dright));
			};
		}
	}
	else if (var.type->kind == TypeKind::Vec2) // gradient w.r.t. a vec2 variable
	{

	}
	else if (var.type->kind == TypeKind::Vec3) // gradient w.r.t. a vec3 variable
	{
		if (binary->op == BinaryOp::Multiply)
		{
			if (isScalarType(ltype) && isScalarType(rtype))
			{
				// the derivatives of scalars with respect to a vec3 variable should be vec3s)
				assert(isVec3Type(dltype) && isVec3Type(drtype));
				// grad( f * g ) = grad(f) * g + f * grad(g), where f and g are scalars
				return Add(Mul(dleft, right), Mul(left, dright));
			}
			if (isScalarType(ltype) && (rtype->kind == TypeKind::Vec3))
			{
				// grad( f * v ) = v * grad(f) + f * grad(v), where f is a scalar and v is a vector
				return Add(OuterProduct(right, dleft), Mul(left, dright));
			}
			else if ((ltype->kind == TypeKind::Vec3) && isScalarType(rtype))
			{
				// grad( v * f ) = grad(v) * f + v * grad(f), where v is a vector and f is a scalar
				return Add(Mul(dleft, right), OuterProduct(left, dright));
			}
			else if ((ltype->kind == TypeKind::Vec3) && (rtype->kind == TypeKind::Vec3))
			{
				// grad( v1 . v2 ) = transpose(grad(v1)) . v2 + transpose(grad(v2)) . v1, where v1 and v2 are vectors
				return Add(Mul(Transpose(dleft), right), Mul(Transpose(dright), left));
			}
		}
		if (binary->op == BinaryOp::Divide)
		{
			if (isScalarType(ltype) && isScalarType(rtype))
			{
				// the derivatives of scalars with respect to a vec3 variable should be vec3s)
				assert(isVec3Type(dltype) && isVec3Type(drtype));
				// grad( f / g ) = (grad(f) * g - f * grad(g)) / (g * g), where f and g are scalars
				return Div(Sub(Mul(dleft, right), Mul(left, dright)), Mul(right, right));
			}
			else if ((ltype->kind == TypeKind::Vec3) && isScalarType(rtype))
			{
				// grad( v / f ) = (grad(v) * f - v & grad(f)) / (f * f), where v is a vector and f is a scalar
				return Div(Sub(Mul(dleft, right), OuterProduct(left, dright)), Mul(right, right));
			}
		}
		if (binary->op == BinaryOp::Exponent)
		{
			if (isScalarType(ltype) && isLiteral(right))
			{
				LiteralExpr* exponent = dynamic_cast<LiteralExpr*>(right.get());
				if (isIntNumber(exponent->value))
				{
					int p = toIntNumber(exponent->value);
					if (p == 1) return clone(dleft.get());
					else if (p != 0)
					{
						// grad(x^p) = p * x^(p-1)*grad(x), where x is a vector and p is a scalar literal
						return Mul(Mul(Literal(p), Pow(left, Literal(p - 1.0))), dleft);
					}
				}
			}
		}
	}

	throw std::runtime_error("AD error: Unsupported binary operator '" + opToString(binary->op) + "' with operand types '" + TypeToString(ltype) + "' and '" + TypeToString(rtype) + "'.");
}

ExprPtr Differentiator::diffCall(const CallExpr* call, const DerivVar& var)
{
	std::vector<Type> argTypes;
	for (const auto& arg : call->arguments)
	{
		argTypes.push_back(arg->valType);
	}

	// find the function with the matching name and argument types in the program's function definitions
	int fnIndex = prg.resolveFunction(call->name, argTypes);
	if (fnIndex == -1)
		throw std::runtime_error("Function \"" + call->name + "\" not found for differentiation.");

	Type returnType = call->valType;
	Type derivType = getDerivativeType(returnType, var.type->kind);

	const std::string& fnc = call->name;
	auto& args = call->arguments;
	if (args.size() == 1)
	{
		// differentiate argument
		auto diffArg = differentiate(args[0].get(), var);

		if      (fnc == "sin") return Mul(Call("cos", args), diffArg);
		else if (fnc == "cos") return Negate(Mul(Call("sin", args), diffArg));
		else if (fnc == "exp") return Mul(Call("exp", args), diffArg);
		else if (fnc == "sqrt") return Mul(Div(Literal(0.5), Call("sqrt", args)), diffArg);
		else if (fnc == "length" && isVec3Type(args[0]->valType))
		{
			if (derivType->kind == TypeKind::Double)
			{
				assert(diffArg->valType->kind == TypeKind::Vec3);
				// d(length(v)) = (v . d(v)) / length(v)
				return Div(Mul(args[0], diffArg), Call("length", args));
			}
			else if (derivType->kind == TypeKind::Vec3)
			{
				assert(diffArg->valType->kind == TypeKind::Mat3);
				// grad(length(v)) = (transpose(grad(v)) . v) / length(v)
				return Div(Mul(Transpose(diffArg), args[0]), Call("length", args));
			}
		}
		else if (fnc == "normalize" && isVec3Type(args[0]->valType))
		{
			// rewrite normalize(v) as v / length(v), then derivate that
			ExprPtr normalize = Div(args[0], Call("length", args));
			return differentiate(normalize.get(), var);
		}
	}

	throw std::runtime_error("Don't know how to differentiate function \"" + fnc + "\".");
}

ExprPtr Differentiator::diffInit(const InitExpr* init, const DerivVar& var)
{
	std::vector<ExprPtr> diffElements;
	for (const auto& elem : init->elements)
	{
		diffElements.push_back(differentiate(elem.get(), var));
	}
	return std::make_unique<InitExpr>(std::move(diffElements));
}

ExprPtr Differentiator::diffConstructor(const ConstructorExpr* ctor, const DerivVar& var)
{
	// determine the type of the derivative variable based on the type of the original variable and the derivative variable.
	Type derivType = getDerivativeType(ctor->valType, var.type->kind);

	if (var.type->kind == TypeKind::Double)
	{
		std::vector<ExprPtr> diffArgs;
		for (const auto& arg : ctor->args)
		{
			diffArgs.push_back(differentiate(arg.get(), var));
		}
		return std::make_unique<ConstructorExpr>(ctor->valType, std::move(diffArgs));
	}
	else if (var.type->kind == TypeKind::Vec3)
	{
		if (ctor->args.size() == 3)
		{
			std::vector<ExprPtr> diffArgs;
			for (const auto& arg : ctor->args)
			{
				ExprPtr darg = differentiate(arg.get(), var);
				if (darg->valType != prg.types.Vec3())
				{
					throw std::runtime_error("Don't know how to differentiate this constructor for vec3 variable.");
				}
				diffArgs.push_back(std::move(darg));
			}

			// put all the components into a mat3 constructor
			// grad(vec3(x, y, z)) = mat3( grad(x),
			//                             grad(y),
			//                             grad(z) )
			return std::make_unique<ConstructorExpr>(prg.types.Mat3(), std::move(diffArgs));
		}
		else
			throw std::runtime_error("Don't know how to differentiate this constructor for vec3 variable.");
	}

	throw std::runtime_error("Don't know how to differentiate constructor.");
}

std::unique_ptr<Expression> Differentiator::diffAssign(const AssignExpr* assign, const DerivVar& var)
{
	// For an assignment expression, we can use the rule: d( y = expr ) --> dy = d(expr)
	auto du = differentiate(assign->target.get(), var);
	auto dv = differentiate(assign->value.get(), var);

	// if the left is zero, that meants that it did not depend on the variable.
	// In that case, we just copy the original expression, since it doesn't contribute to the derivative.
	if (isZero(du))
		return clone(assign);

	return Assign(du, dv);
}

std::unique_ptr<Expression> Differentiator::diffMember(const MemberExpr* member, const DerivVar& var)
{
	if (var.type->kind == TypeKind::Double)
	{
		// For a member access expression, we can use the rule: d( obj.field ) --> dobj.field
		auto dobj = differentiate(member->object.get(), var);
		return Member(dobj, member->property);
	}
	else if (var.type->kind == TypeKind::Vec2)
	{
		if (isVariable(member->object))
		{
			auto var = dynamic_cast<const VariableExpr*>(member->object.get());
			if (var->name == var->name)
			{
				if      (member->property == "x") return Literal(vec2(1, 0));
				else if (member->property == "y") return Literal(vec2(0, 1));
				else
					throw std::runtime_error("Don't know how to differentiate this member access for vec2 variable.");
			}
		}
	}
	else if (var.type->kind == TypeKind::Vec3)
	{
		if (isVariable(member->object))
		{
			auto var = dynamic_cast<const VariableExpr*>(member->object.get());
			if (var->name == var->name)
			{
				if      (member->property == "x") return Literal(vec3(1, 0, 0));
				else if (member->property == "y") return Literal(vec3(0, 1, 0));
				else if (member->property == "z") return Literal(vec3(0, 0, 1));
				else
					throw std::runtime_error("Don't know how to differentiate this member access for vec3 variable.");
			}
		}
	}

	throw std::runtime_error("Don't know how to differentiate member access for this variable type.");
}

std::unique_ptr<Expression> Differentiator::diffIndex(const IndexExpr* index, const DerivVar& var)
{
	if (var.type->kind == TypeKind::Double)
	{
		// For an index access expression, we can use the rule: d( obj[i] ) --> d(obj)[i]
		auto dobj = differentiate(index->object.get(), var);
		return Index(dobj, index->index);
	}

	throw std::runtime_error("Don't know how to differentiate index access for this variable type.");
}

bool febcode::dependsOn(const Expression* expr, const std::string& var)
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
	else if (auto index = dynamic_cast<const IndexExpr*>(expr))
	{
		return dependsOn(index->object.get(), var) || dependsOn(index->index.get(), var);
	}
	else if (auto assign = dynamic_cast<const AssignExpr*>(expr))
	{
		return dependsOn(assign->target.get(), var) || dependsOn(assign->value.get(), var);
	}
	else if (auto init = dynamic_cast<const InitExpr*>(expr))
	{
		for (const auto& elem : init->elements)
		{
			if (dependsOn(elem.get(), var))
				return true;
		}
		return false;
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
		if (!varDeclStmt->input)
		{
			for (const auto& var : varDeclStmt->vars)
			{
				if (::dependsOn(var->initializer.get(), varName))
					return true;
			}
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
