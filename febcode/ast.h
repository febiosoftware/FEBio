#pragma once
#include <string>
#include <vector>
#include "types.h"

namespace febcode {

	enum class BinaryOp {
		Plus, Minus, Multiply, Divide, Exponent,
		EqualEqual, NotEqual,
		Greater, GreaterEqual,
		Less, LessEqual, 
		AndAnd, OrOr
	};

	enum class UnaryOp {
		Negate, // "-"
		Not	 // "!"
	};

	// Forward declarations
	struct Expression;
	struct Statement;

	// Base classes
	struct Expression {
		virtual ~Expression() = default;
	};

	struct Statement {
		virtual ~Statement() = default;
	};

	// --- Expressions ---

	struct LiteralExpr : Expression {
		Value value;
		LiteralExpr(const Value& v) : value(v) {}
	};

	struct VariableExpr : Expression {
		std::string name;
		VariableExpr(const std::string& n) : name(n) {}
	};

	struct AssignExpr : Expression {
		std::string name;
		std::unique_ptr<Expression> value;
		AssignExpr(const std::string& n, std::unique_ptr<Expression> v)
			: name(n), value(std::move(v)) {
		}
	};

	struct BinaryExpr : Expression {
		BinaryOp op;
		std::unique_ptr<Expression> left;
		std::unique_ptr<Expression> right;

		BinaryExpr(std::unique_ptr<Expression> l, BinaryOp o, std::unique_ptr<Expression> r)
			: left(std::move(l)), op(o), right(std::move(r)) {
		}
	};

	struct UnaryExpr : Expression {
		UnaryOp op;
		std::unique_ptr<Expression> right;

		UnaryExpr(UnaryOp o, std::unique_ptr<Expression> r)
			: op(o), right(std::move(r)) {
		}
	};

	struct CallExpr : Expression {
		std::unique_ptr<Expression> callee;
		std::vector<std::unique_ptr<Expression>> arguments;

		CallExpr(std::unique_ptr<Expression> callee,
			std::vector<std::unique_ptr<Expression>> arguments)
			: callee(std::move(callee)),
			arguments(std::move(arguments)) {
		}
	};

	struct MemberExpr : Expression
	{
		std::unique_ptr<Expression> object;
		std::string property;

		MemberExpr(std::unique_ptr<Expression> object,
			std::string property)
			: object(std::move(object)),
			property(std::move(property)) {
		}
	};

	struct SetMemberExpr : Expression
	{
		std::unique_ptr<Expression> object;
		std::string property;
		std::unique_ptr<Expression> value;

		SetMemberExpr(std::unique_ptr<Expression> object,
			std::string property,
			std::unique_ptr<Expression> value)
			: object(std::move(object)),
			property(std::move(property)),
			value(std::move(value)) {
		}
	};

	struct InitializerExpr : Expression {
		std::vector<std::unique_ptr<Expression>> elements;
		InitializerExpr(std::vector<std::unique_ptr<Expression>> elements)
			: elements(std::move(elements)) {
		}
	};

	struct IndexExpr : Expression {
		std::unique_ptr<Expression> object;
		std::unique_ptr<Expression> index;
		IndexExpr(std::unique_ptr<Expression> object,
			std::unique_ptr<Expression> index)
			: object(std::move(object)),
			index(std::move(index)) {
		}
	};

	struct SetIndexExpr : Expression {
		std::unique_ptr<Expression> object;
		std::unique_ptr<Expression> index;
		std::unique_ptr<Expression> value;
		SetIndexExpr(std::unique_ptr<Expression> object,
			std::unique_ptr<Expression> index,
			std::unique_ptr<Expression> value)
			: object(std::move(object)),
			index(std::move(index)),
			value(std::move(value)) {
		}
	};

	// --- Statements ---

	struct ExpressionStmt : Statement {
		std::unique_ptr<Expression> expr;
		ExpressionStmt(std::unique_ptr<Expression> e) : expr(std::move(e)) {}
	};

	struct Var {
		Type type = nullptr;
		std::string name;
		std::unique_ptr<Expression> initializer; // can be null
	};

	struct VarDeclStmt : Statement {
		std::vector<Var> vars;
		VarDeclStmt(std::vector<Var>& vars)
			: vars(std::move(vars)) {
		}
	};

	struct ReturnStmt : Statement {
		std::unique_ptr<Expression> value; // can be null
		ReturnStmt(std::unique_ptr<Expression> v) : value(std::move(v)) {}
	};

	struct BlockStmt : Statement {
		std::vector<std::unique_ptr<Statement>> statements;
	};

	struct IfStmt : Statement {
		std::unique_ptr<Expression> condition;
		std::unique_ptr<Statement> thenBranch;
		std::unique_ptr<Statement> elseBranch; // optional

		IfStmt(std::unique_ptr<Expression> cond,
			std::unique_ptr<Statement> thenStmt,
			std::unique_ptr<Statement> elseStmt = nullptr)
			: condition(std::move(cond)),
			thenBranch(std::move(thenStmt)),
			elseBranch(std::move(elseStmt)) {
		}
	};

	struct WhileStmt : Statement {
		std::unique_ptr<Expression> condition;
		std::unique_ptr<Statement> body;

		WhileStmt(std::unique_ptr<Expression> cond,
			std::unique_ptr<Statement> bodyStmt)
			: condition(std::move(cond)),
			body(std::move(bodyStmt)) {
		}
	};

	struct FunctionStmt : Statement {
		std::string name;
		Type returnType = nullptr;
		std::vector<std::pair<Type,std::string>> params;
		std::unique_ptr<Statement> body;

		FunctionStmt(std::string name, Type returnType,
			std::vector<std::pair<Type,std::string>> parameters,
			std::unique_ptr<Statement> bdy)
			: name(std::move(name)),
			returnType(returnType),
			params(std::move(parameters)),
			body(std::move(bdy)) {
		}
	};

	struct StructStmt : Statement {
		std::string name;
		Type type;
		std::vector<std::pair<Type, std::string>> fields;
		StructStmt(std::string name, Type type,
			std::vector<std::pair<Type, std::string>> fields)
			: name(std::move(name)),
			type(type),
			fields(std::move(fields)) {
		}
	};

	struct AST {
		std::vector<std::unique_ptr<Statement>> statements;

		AST() {}

		void clear() { statements.clear(); }

		bool empty() const { return statements.empty(); }

		size_t size() const { return statements.size(); }

		Statement* operator[](size_t index) const {
			if (index >= statements.size()) {
				return nullptr;
			}
			return statements[index].get();
		}
	};
} // namespace febcode

std::ostream& operator << (std::ostream& o, const febcode::Value& v);

std::string ValueToString(const febcode::Value& v);

std::string ValueTypeToString(const febcode::Value& v);

