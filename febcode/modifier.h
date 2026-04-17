#pragma once
#include "program.h"
#include <vector>
#include <memory>

namespace febcode
{
	// base class for AST modifiers. This is used for both differentiation and simplification, since they both modify the AST in some way. 
	// It provides a common interface for applying modifications to an AST, and it can be extended to implement specific modifications like differentiation or simplification.
	class Modifier
	{
	public:
		Modifier(Program& prg) : prg(prg) {}

	protected:
		// helper functions for creating (typed) expressions more easily
		ExprPtr Literal(bool   b) { ExprPtr e = std::make_unique<LiteralExpr>(Value(b)); e->valType = prg.types.Bool  (); return e; }
		ExprPtr Literal(int    n) { ExprPtr e = std::make_unique<LiteralExpr>(Value(n)); e->valType = prg.types.Int   (); return e; }
		ExprPtr Literal(double a) { ExprPtr e = std::make_unique<LiteralExpr>(Value(a)); e->valType = prg.types.Double(); return e; }
		ExprPtr Literal(vec2   a) { ExprPtr e = std::make_unique<LiteralExpr>(Value(a)); e->valType = prg.types.Vec2  (); return e; }
		ExprPtr Literal(vec3   a) { ExprPtr e = std::make_unique<LiteralExpr>(Value(a)); e->valType = prg.types.Vec3  (); return e; }
		ExprPtr Literal(mat2   a) { ExprPtr e = std::make_unique<LiteralExpr>(Value(a)); e->valType = prg.types.Mat2  (); return e; }
		ExprPtr Literal(mat3   a) { ExprPtr e = std::make_unique<LiteralExpr>(Value(a)); e->valType = prg.types.Mat3  (); return e; }

		ExprPtr Variable(const std::string& name, Type type) { ExprPtr v = std::make_unique<VariableExpr>(name); v->valType = type; return v; }

		ExprPtr Assign(const ExprPtr& target, const ExprPtr& value) { ExprPtr a = std::make_unique<AssignExpr>(copy_expression(target.get()), copy_expression(value.get())); a->valType = target->valType; return a; }

		ExprPtr Call(const std::string& name, const std::vector<ExprPtr>& args);

		ExprPtr Index(const ExprPtr& object, const ExprPtr& index);

		ExprPtr Member(const ExprPtr& object, const std::string& property);

		ExprPtr Negate(const ExprPtr& arg) { ExprPtr n = std::make_unique<UnaryExpr>(UnaryOp::Negate, copy_expression(arg.get())); n->valType = arg->valType; return n; }

		// create zero initializer for array
		ExprPtr Initializer(Type type);

		// create zero constructor for struct
		ExprPtr Constructor(Type type);

		// make a literal expression of type with zero values.
		ExprPtr Zero(Type type);

		// create a binary expression with the given operator and operands.
		ExprPtr Binary(BinaryOp op, const ExprPtr& left, const ExprPtr& right);

		febcode::ExprPtr OuterProduct(const febcode::ExprPtr& left, const febcode::ExprPtr& right)
		{
			std::vector<febcode::ExprPtr> args(2);
			args[0] = copy_expression(left.get());
			args[1] = copy_expression(right.get());
			return Call("outer", std::move(args));
		}

		febcode::ExprPtr Transpose(const febcode::ExprPtr& arg)
		{
			std::vector<febcode::ExprPtr> args(1);
			args[0] = copy_expression(arg.get());
			return Call("transpose", std::move(args));
		}

		febcode::ExprPtr Inverse(const febcode::ExprPtr& arg)
		{
			std::vector<febcode::ExprPtr> args(1);
			args[0] = copy_expression(arg.get());
			return Call("inverse", std::move(args));
		}

		febcode::ExprPtr Pow(const febcode::ExprPtr& left, const febcode::ExprPtr& right)
		{ 
			return Binary(febcode::BinaryOp::Exponent, left, right); 
		}

		febcode::ExprPtr Add(const febcode::ExprPtr& left, const febcode::ExprPtr& right) { return Binary(febcode::BinaryOp::Plus    , left, right); }
		febcode::ExprPtr Sub(const febcode::ExprPtr& left, const febcode::ExprPtr& right) { return Binary(febcode::BinaryOp::Minus   , left, right); }
		febcode::ExprPtr Mul(const febcode::ExprPtr& left, const febcode::ExprPtr& right) { return Binary(febcode::BinaryOp::Multiply, left, right); }
		febcode::ExprPtr Div(const febcode::ExprPtr& left, const febcode::ExprPtr& right) { return Binary(febcode::BinaryOp::Divide  , left, right); }

	protected:
		Program& prg;
	};
}
