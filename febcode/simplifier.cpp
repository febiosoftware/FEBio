#include "simplifier.h"
using namespace febcode;

Simplifier::Simplifier(Program& prg) : Modifier(prg) 
{
	// --- constant folding rules ---
	// integer constant folding
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression *a, *b;
		int l, r;
		if (isAdd(e, a, b) && isInt(a, l) && isInt(b, r)) return Literal(l + r);
		if (isSub(e, a, b) && isInt(a, l) && isInt(b, r)) return Literal(l - r);
		if (isMul(e, a, b) && isInt(a, l) && isInt(b, r)) return Literal(l * r);
		if (isDiv(e, a, b) && isInt(a, l) && isInt(b, r)) return Literal(l / r);
		if (isExp(e, a, b) && isInt(a, l) && isInt(b, r)) return Literal(static_cast<int>(std::pow(l, r)));
		if (isNegate(e, a) && isInt(a, l)) return Literal(-l);

		return nullptr;
	});

	// scalar constant folding
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;
		double l, r;
		if (isAdd(e, a, b) && isScalar(a, l) && isScalar(b, r)) return Literal(l + r);
		if (isSub(e, a, b) && isScalar(a, l) && isScalar(b, r)) return Literal(l - r);
		if (isMul(e, a, b) && isScalar(a, l) && isScalar(b, r)) return Literal(l * r);
		if (isExp(e, a, b) && isScalar(a, l) && isScalar(b, r)) return Literal(std::pow(l, r));
		if (isDiv(e, a, b) && isScalar(a, l) && isScalar(b, r))
		{
			if (r == 0.0)
				throw std::runtime_error("Division by zero in constant folding");
			return Literal(l / r);
		}
		if (isNegate(e, a) && isDouble(a, l)) return Literal(-l);
		return nullptr;
	});

	// vec2 constant folding for binary operations
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;
		
		vec2 l, r;
		if (isAdd(e, a, b) && isVec2(a, l) && isVec2(b, r)) return Literal(l + r);
		if (isSub(e, a, b) && isVec2(a, l) && isVec2(b, r)) return Literal(l - r);
		if (isMul(e, a, b) && isVec2(a, l) && isVec2(b, r)) return Literal(l * r);

		double scalar;
		if (isMul(e, a, b) && isVec2(a, l) && isScalar(b, scalar)) return Literal(l * scalar);
		if (isDiv(e, a, b) && isVec2(a, l) && isScalar(b, scalar)) return Literal(l / scalar);

		if (isMul(e, a, b) && isScalar(a, scalar) && isVec2(b, r)) return Literal(r * scalar);

		return nullptr;
	});

	// vec3 constant folding for binary operations
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;

		vec3 l, r;
		if (isAdd(e, a, b) && isVec3(a, l) && isVec3(b, r)) return Literal(l + r);
		if (isSub(e, a, b) && isVec3(a, l) && isVec3(b, r)) return Literal(l - r);
		if (isMul(e, a, b) && isVec3(a, l) && isVec3(b, r)) return Literal(l * r);

		double scalar;
		if (isMul(e, a, b) && isVec3(a, l) && isScalar(b, scalar)) return Literal(l * scalar);
		if (isDiv(e, a, b) && isVec3(a, l) && isScalar(b, scalar)) return Literal(l / scalar);

		if (isMul(e, a, b) && isScalar(a, scalar) && isVec3(b, r)) return Literal(r * scalar);

		return nullptr;
	});

	// mat2 constant folding for binary operations
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;

		mat2 l, r;
		if (isAdd(e, a, b) && isMat2(a, l) && isMat2(b, r)) return Literal(l + r);
		if (isSub(e, a, b) && isMat2(a, l) && isMat2(b, r)) return Literal(l - r);
		if (isMul(e, a, b) && isMat2(a, l) && isMat2(b, r)) return Literal(l * r);

		double scalar;
		if (isMul(e, a, b) && isMat2(a, l) && isScalar(b, scalar)) return Literal(l * scalar);
		if (isDiv(e, a, b) && isMat2(a, l) && isScalar(b, scalar)) return Literal(l / scalar);
		if (isMul(e, a, b) && isScalar(a, scalar) && isMat2(b, r)) return Literal(r * scalar);

		vec2 v;
		if (isMul(e, a, b) && isMat2(a, l) && isVec2(b, v)) return Literal(l * v);

		return nullptr;
	});

	// mat3 constant folding for binary operations
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;

		mat3 l, r;
		if (isAdd(e, a, b) && isMat3(a, l) && isMat3(b, r)) return Literal(l + r);
		if (isSub(e, a, b) && isMat3(a, l) && isMat3(b, r)) return Literal(l - r);
		if (isMul(e, a, b) && isMat3(a, l) && isMat3(b, r)) return Literal(l * r);

		double scalar;
		if (isMul(e, a, b) && isMat3(a, l) && isScalar(b, scalar)) return Literal(l * scalar);
		if (isDiv(e, a, b) && isMat3(a, l) && isScalar(b, scalar)) return Literal(l / scalar);
		if (isMul(e, a, b) && isScalar(a, scalar) && isMat3(b, r)) return Literal(r * scalar);

		vec3 v;
		if (isMul(e, a, b) && isMat3(a, l) && isVec3(b, v)) return Literal(l * v);

		return nullptr;
	});

	// function evaluation for built-in functions with constant arguments
	rules.push_back([this](const Expression* e) -> ExprPtr {
		if (auto call = dynamic_cast<const CallExpr*>(e))
		{
			std::vector<double> scalarArgs;
			for (const auto& arg : call->arguments)
			{
				double value;
				if (!isScalar(arg.get(), value)) return nullptr; // not all arguments are scalar constants
				scalarArgs.push_back(value);
			}

			VariableExpr* calleeVar = dynamic_cast<VariableExpr*>(call->callee.get());
			if (!calleeVar) 
				febcode::error(call, "Unsupported function call in constant folding: callee is not a variable");

			const std::string& fname = calleeVar->name;
			if (scalarArgs.size() == 1)
			{
				// Make sure this list matches the built-in functions defined in the math module!
				if (fname == "abs"  ) return Literal(std::abs  (scalarArgs[0]));
				if (fname == "acos" ) return Literal(std::acos (scalarArgs[0]));
				if (fname == "acosh") return Literal(std::acosh(scalarArgs[0]));
				if (fname == "asin" ) return Literal(std::asin (scalarArgs[0]));
				if (fname == "asinh") return Literal(std::asinh(scalarArgs[0]));
				if (fname == "atan" ) return Literal(std::atan (scalarArgs[0]));
				if (fname == "atanh") return Literal(std::atanh(scalarArgs[0]));
				if (fname == "cos"  ) return Literal(std::cos  (scalarArgs[0]));
				if (fname == "cosh" ) return Literal(std::cosh (scalarArgs[0]));
				if (fname == "exp"  ) return Literal(std::exp  (scalarArgs[0]));
				if (fname == "log"  ) return Literal(std::log  (scalarArgs[0]));
				if (fname == "log10") return Literal(std::log10(scalarArgs[0]));
				if (fname == "sin"  ) return Literal(std::sin  (scalarArgs[0]));
				if (fname == "sinh" ) return Literal(std::sinh (scalarArgs[0]));
				if (fname == "sqrt" ) return Literal(std::sqrt (scalarArgs[0]));
				if (fname == "tan"  ) return Literal(std::tan  (scalarArgs[0]));
				if (fname == "tanh" ) return Literal(std::tanh (scalarArgs[0]));
			}
		}
		return nullptr;
	});

	// --- algebraic simplification rules ---
	
	// -0 = 0
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a;
		if (isNegate(e, a) && isZero(a)) return Zero(e->valType); // -0 = 0
		return nullptr;
	});

	// --x = x
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression *a, *b;
		if (isNegate(e, a) && isNegate(a, b)) return simplify(b); // --x = x
		return nullptr;
	});
	
	// x +/- 0 = x
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression *a, *b;
		if (isAdd(e, a, b))
		{
			if (isZero(a)) return simplify(b); // 0 + x = x
			if (isZero(b)) return simplify(a); // x + 0 = x
		}
		if (isSub(e, a, b))
		{
			if (isZero(a)) return simplify(Negate(b)); // 0 - x = -x
			if (isZero(b)) return simplify(a); // x - 0 = x
		}
		return nullptr;
	});

	// a + (-b) = a - b, a - (-b) = a + b
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression *a, *n, *b;
		if (isAdd(e, a, n) && isNegate(n, b))
		{
			 return simplify(Binary(BinaryOp::Minus, simplify(a), simplify(b))); // a + (-b) = a - b
		}
		if (isSub(e, a, n) && isNegate(n, b))
		{
			return simplify(Binary(BinaryOp::Plus, simplify(a), simplify(b))); // a - (-b) = a + b
		}
		return nullptr;
	});

	// x * 0 = 0
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression *a, *b;
		if (isMul(e, a, b))
		{
			if (isZero(a)) return Zero(e->valType); // 0 * x = 0
			if (isZero(b)) return Zero(e->valType); // x * 0 = 0
		}
		return nullptr;
		});

	// x * 1 = x
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression *a, *b;
		if (isMul(e, a, b))
		{
			if (isOne(a)) return simplify(b); // 1 * x = x
			if (isOne(b)) return simplify(a); // x * 1 = x
		}
		return nullptr;
	});

	// x / 1 = x
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;
		if (isDiv(e, a, b))
		{
			if (isOne(b)) return simplify(a); // x / 1 = x
		}
		return nullptr;
	});

	// 0 / x = 0
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;
		if (isDiv(e, a, b))
		{
			if (isZero(a) && !isZero(b)) return Zero(e->valType); // 0 / x = 0
		}
		return nullptr;
	});

	// x / 0 --> throws
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;
		if (isDiv(e, a, b))
		{
			if (!isZero(a) && isZero(b)) throw std::runtime_error("Division by zero in algebraic simplification");
		}
		return nullptr;
	});

	// x ^ 1 = x, x ^ 0 = 1
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;
		if (isExp(e, a, b))
		{
			if (isOne(b)) return simplify(a); // x ^ 1 = x

			double l, r;
			if (isScalar(a, l) && isScalar(b, r))
			{
				if ((l != 0.0) && (r == 0.0)) return Literal(1.0); // x ^ 0 = 1
				if ((l == 1.0) && (r != 0.0)) return Literal(1.0); // 1 ^ x = 1
			} 
		}
		return nullptr;
	});

	// I * x = x, x * I = x
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;
		if (isMul(e, a, b))
		{
			if (isIdentity(a) && isVector(b)) return simplify(b); // I * x = x
			if (isIdentity(a) && isMatrix(b)) return simplify(b); // I * x = x
			if (isIdentity(b) && isMatrix(a)) return simplify(a); // x * I = x
		}
		return nullptr;
	});

	// a op a
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression* a, * b;
		if (isAdd(e, a, b) && isEqual(a, b)) return simplify(Binary(BinaryOp::Multiply, Literal(2), simplify(a))); // a + a = 2 * a
		if (isSub(e, a, b) && isEqual(a, b)) return Zero(e->valType); // a - a = 0
		if (isMul(e, a, b) && isEqual(a, b) && isScalar(a)) return simplify(Binary(BinaryOp::Exponent, simplify(a), Literal(2))); // a * a = a ^ 2
		if (isDiv(e, a, b) && isEqual(a, b) && isScalar(a)) return Literal(1.0); // a / a = 1
		return nullptr;
	});

	// --- rules for built-in functions ---
	// dot(0,x) = 0
	rules.push_back([this](const Expression* e) -> ExprPtr {
		const Expression* a, * b;
		if (isCall2(e, "dot", a, b) && (isZero(a) || isZero(b))) return Literal(0.0);
		return nullptr;
	});

	// length(0) = 0
	rules.push_back([this](const Expression* e) -> ExprPtr {
		const Expression* arg;
		if (isCall1(e, "length", arg) && isZero(arg)) return Literal(0.0);
		return nullptr;
	});

	// outer(0,x) = 0
	rules.push_back([this](const Expression* e) -> ExprPtr {
		const Expression* a, * b;
		if (isCall2(e, "outer", a, b) && (isZero(a) || isZero(b))) return Zero(e->valType);
		return nullptr;
	});

	// cross(0,x) = 0
	rules.push_back([this](const Expression* e) -> ExprPtr {
		const Expression* a, * b;
		if (isCall2(e, "cross", a, b) && (isZero(a) || isZero(b))) return Zero(e->valType);
		return nullptr;
	});

	// transpose(transpose(x)) = x
	rules.push_back([this](const Expression* e) -> ExprPtr {
		const Expression* arg;
		if (isCall1(e, "transpose", arg) && isCall1(arg, "transpose", arg))
		{
			return simplify(arg);
		}
		return nullptr;
	});

	// transpose(symmetric) = symmetric
	rules.push_back([this](const Expression* e) -> ExprPtr {
		const Expression* arg;
		if (isCall1(e, "transpose", arg) && isSymmetric(arg))
		{
			return simplify(arg);
		}
		return nullptr;
	});

	// --- rules for moving towards canonical forms ---
	// scalar * (scalar * x) = (scalar * scalar) * x
	rules.push_back([this](const Expression* e) -> ExprPtr {
		Expression *a, *m, *b, *c;
		double l, r;
		if (isMul(e, a, m) && isScalar(a, l) && isMul(m, b, c) && isScalar(b, r))
		{
			return simplify(Mul(Literal(l * r).get(), c));
		}
		if (isMul(e, a, m) && isScalar(a, l) && isMul(m, b, c) && isScalar(c, r))
		{
			return simplify(Mul(Literal(l * r).get(), b));
		}
		return nullptr;
	});
}

ExprPtr Simplifier::applyRules(const Expression* expr)
{
	for (const auto& rule : rules)
	{
		if (auto result = rule(expr))
			return result;
	}
	return nullptr;
}

ExprPtr Simplifier::simplify(const Expression* expr)
{
	if (auto ctor   = dynamic_cast<const ConstructorExpr*>(expr)) return simplifyConstructor(ctor  );
	if (auto unary  = dynamic_cast<const UnaryExpr*      >(expr)) return simplifyUnary      (unary );
	if (auto init   = dynamic_cast<const InitExpr*       >(expr)) return simplifyInitializer(init  );
	if (auto call   = dynamic_cast<const CallExpr*       >(expr)) return simplifyCall       (call  );
	if (auto binary = dynamic_cast<const BinaryExpr*     >(expr)) return simplifyBinary     (binary);
	if (auto assign = dynamic_cast<const AssignExpr*     >(expr)) return simplifyAssign     (assign);
	return clone(expr);
}

ExprPtr Simplifier::simplifyConstructor(const ConstructorExpr* ctor)
{
	// simplify all arguments first
	std::vector<ExprPtr> simplifiedArgs;
	for (const auto& arg : ctor->args)
	{
		simplifiedArgs.push_back(simplify(arg.get()));
	}

	// if all args are numbers, we'll make this a literal too
	if (std::all_of(simplifiedArgs.begin(), simplifiedArgs.end(), [](const ExprPtr& arg) {
		auto lit = dynamic_cast<const LiteralExpr*>(arg.get());
		return lit && (isDouble(lit->value) || isInt(lit->value));
		}))
	{
		std::vector<double> v;
		for (const auto& arg : simplifiedArgs)
		{
			auto lit = dynamic_cast<const LiteralExpr*>(arg.get());
			if (isInt(lit->value))
				v.push_back(lit->value.i);
			else
				v.push_back(lit->value.d);
		}

		TypeKind type = ctor->valType->kind;
		if ((type == TypeKind::Vec2) && (v.size() == 1)) return Literal(vec2(v[0], v[0]));
		if ((type == TypeKind::Vec2) && (v.size() == 2)) return Literal(vec2(v[0], v[1]));
		if ((type == TypeKind::Vec3) && (v.size() == 1)) return Literal(vec3(v[0], v[0], v[0]));
		if ((type == TypeKind::Vec3) && (v.size() == 3)) return Literal(vec3(v[0], v[1], v[2]));
		if ((type == TypeKind::Mat2) && (v.size() == 1)) return Literal(mat2(v[0], 0.0, 0.0, v[0]));
		if ((type == TypeKind::Mat2) && (v.size() == 4)) return Literal(mat2(v[0], v[1], v[2], v[3]));
		if ((type == TypeKind::Mat3) && (v.size() == 1)) return Literal(mat3(v[0], 0.0, 0.0, 0.0, v[0], 0.0, 0.0, 0.0, v[0]));
		if ((type == TypeKind::Mat3) && (v.size() == 9)) return Literal(mat3(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]));
	}

	return std::make_unique<ConstructorExpr>(ctor->valType, std::move(simplifiedArgs));
}

ExprPtr Simplifier::simplifyUnary(const UnaryExpr* unary)
{
	auto simplifiedRight = simplify(unary->right.get());
	ExprPtr tmp = Unary(unary->op, simplifiedRight);
	if (auto result = applyRules(tmp.get()))
	{
		return result;
	}
	return std::move(tmp);
}

ExprPtr Simplifier::simplifyInitializer(const InitExpr* init)
{
	std::vector<ExprPtr> simplifiedElements;
	for (const auto& elem : init->elements)
	{
		simplifiedElements.push_back(simplify(elem.get()));
	}
	return std::make_unique<InitExpr>(std::move(simplifiedElements));
}

ExprPtr Simplifier::simplifyCall(const CallExpr* call)
{
	std::vector<ExprPtr> simplifiedArgs;
	for (const auto& arg : call->arguments)
	{
		simplifiedArgs.push_back(simplify(arg.get()));
	}

	VariableExpr* calleeVar = dynamic_cast<VariableExpr*>(call->callee.get());
	if (!calleeVar)
		febcode::error(call, "Unsupported function call in constant folding: callee is not a variable");

	ExprPtr tmp = Call(calleeVar->name, simplifiedArgs);
	if (auto result = applyRules(tmp.get()))
	{
		return result;
	}

	return std::move(tmp);
}

ExprPtr Simplifier::simplifyBinary(const BinaryExpr* binary)
{
	auto l = simplify(binary->left.get());
	auto r = simplify(binary->right.get());

	ExprPtr tmp = Binary(binary->op, l, r);
	if (auto result = applyRules(tmp.get()))
	{
		return result;
	}

	return std::move(tmp);
}

ExprPtr Simplifier::simplifyAssign(const AssignExpr* assign)
{
	auto l = simplify(assign->target.get());
	auto r = simplify(assign->value.get());

	ExprPtr tmp = Assign(l, r);
	if (auto result = applyRules(tmp.get()))
	{
		return result;
	}

	return std::move(tmp);
}