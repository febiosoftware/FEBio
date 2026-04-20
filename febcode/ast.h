#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cmath>
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

	enum class ExpressionType {
		Literal,
		Variable,
		Binary,
		Unary,
		Call,
		Member,
		Initializer,
		Index,
		Assignment,
		Constructor
	};

	// --- Expressions ---

	// Base classes
	struct Expression {
		Expression(ExpressionType type) : exprType(type) {}
		virtual ~Expression() = default;
		ExpressionType exprType;
		Type valType = nullptr; // Value type of expression. Determined during resolution
	};

	using ExprPtr = std::unique_ptr<Expression>;

	struct LiteralExpr : Expression {
		Value value;
		LiteralExpr(const Value& v) : Expression(ExpressionType::Literal), value(v) {}
	};

	struct VariableExpr : Expression {
		std::string name;
		VariableExpr(const std::string& n) : Expression(ExpressionType::Variable), name(n) {}
	};

	class AssignExpr : public Expression
	{
	public:
		ExprPtr target;
		ExprPtr value;

	public:
		AssignExpr(ExprPtr target, ExprPtr value)
			: Expression(ExpressionType::Assignment), target(std::move(target)), value(std::move(value)) {
		}
	};

	struct BinaryExpr : Expression {
		BinaryOp op;
		ExprPtr left;
		ExprPtr right;

		BinaryExpr(ExprPtr l, BinaryOp o, ExprPtr r)
			: Expression(ExpressionType::Binary), left(std::move(l)), op(o), right(std::move(r)) {
		}
	};

	struct UnaryExpr : Expression {
		UnaryOp op;
		ExprPtr right;

		UnaryExpr(UnaryOp o, ExprPtr r)
			: Expression(ExpressionType::Unary), op(o), right(std::move(r)) {}
	};

	struct CallExpr : Expression {
		std::string name;
		std::vector<ExprPtr> arguments;

		CallExpr(const std::string& name, std::vector<ExprPtr> arguments)
			: Expression(ExpressionType::Call), name(name),
			arguments(std::move(arguments)) {}
	};

	struct MemberExpr : Expression
	{
		ExprPtr object;
		std::string property;

		MemberExpr(ExprPtr object, std::string property)
			: Expression(ExpressionType::Member), object(std::move(object)), property(std::move(property)) {}
	};

	struct InitExpr : Expression {
		std::vector<ExprPtr> elements;
		InitExpr(std::vector<ExprPtr> elements)
			: Expression(ExpressionType::Initializer), elements(std::move(elements)) {}
	};

	struct IndexExpr : Expression {
		ExprPtr object;
		ExprPtr index;
		IndexExpr(ExprPtr object, ExprPtr index)
			: Expression(ExpressionType::Index), object(std::move(object)), index(std::move(index)) {}
	};

	struct ConstructorExpr : Expression {
		std::vector<ExprPtr> args;
		ConstructorExpr(Type type, std::vector<ExprPtr> arguments) :
			Expression(ExpressionType::Constructor), args(std::move(arguments)) {
			valType = type;
		}
	};

	// --- Statements ---

	struct Statement {
		virtual ~Statement() = default;
	};
	using StmtPtr = std::unique_ptr<Statement>;

	struct ExpressionStmt : Statement {
		ExprPtr expr;
		ExpressionStmt(ExprPtr e) : expr(std::move(e)) {}
	};

	struct Var {
		std::string name;
		std::vector<size_t> arraySizes; // empty if not an array
		ExprPtr initializer; // can be null
	};

	struct VarDeclStmt : Statement {
		Type type = nullptr;
		bool input = false; // declare input variables (treated as const, non-differentiable)
		std::vector<Var> vars;
		VarDeclStmt(Type type, const std::string& name, ExprPtr initializer, bool input = false) : type(type), input(input)
		{
			vars.push_back({ name, std::vector<size_t>(), std::move(initializer)});
		}
		VarDeclStmt(Type type, Var& var, bool input = false)
			: type(type), input(input) { 
			vars.emplace_back(std::move(var));
		}
		VarDeclStmt(Type type, std::vector<Var>& vars, bool input = false)
			: type(type), input(input), vars(std::move(vars)) {}
	};

	struct ReturnStmt : Statement {
		ExprPtr value; // can be null
		ReturnStmt(ExprPtr v) : value(std::move(v)) {}
	};

	struct BlockStmt : Statement {
		std::vector<StmtPtr> statements;

		void addStatement(StmtPtr stmt) {
			statements.push_back(std::move(stmt));
		}
	};

	struct IfStmt : Statement {
		ExprPtr condition;
		StmtPtr thenBranch;
		StmtPtr elseBranch; // optional

		IfStmt() {}

		IfStmt(ExprPtr cond,
			StmtPtr thenStmt,
			StmtPtr elseStmt = nullptr)
			: condition(std::move(cond)),
			thenBranch(std::move(thenStmt)),
			elseBranch(std::move(elseStmt)) {
		}
	};

	struct WhileStmt : Statement {
		ExprPtr condition;
		StmtPtr body;

		WhileStmt(ExprPtr cond,
			std::unique_ptr<Statement> bodyStmt)
			: condition(std::move(cond)),
			body(std::move(bodyStmt)) {
		}
	};

	struct ForStmt : Statement {
		StmtPtr initializer; // can be null
		ExprPtr condition;
		ExprPtr increment;
		StmtPtr body;

		ForStmt(std::unique_ptr<Statement> init,
			ExprPtr cond,
			ExprPtr incr,
			std::unique_ptr<Statement> bodyStmt)
			: initializer(std::move(init)),
			condition(std::move(cond)),
			increment(std::move(incr)),
			body(std::move(bodyStmt)) {
		}
	};

	struct FunctionStmt : Statement {
		std::string name;
		Type returnType = nullptr;
		std::vector<std::pair<Type,std::string>> params;
		StmtPtr body;

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
		BlockStmt root;

		AST() {}

		void clear() { root.statements.clear(); }

		bool empty() const { return root.statements.empty(); }

		size_t size() const { return root.statements.size(); }

		Statement* operator[](size_t index) const {
			if (index >= root.statements.size()) {
				return nullptr;
			}
			return root.statements[index].get();
		}

		void addStatement(std::unique_ptr<Statement> stmt) {
			root.statements.push_back(std::move(stmt));
		}
	};

	// some helper variables for statements
	inline bool isVarDecl (const StmtPtr& stmt) { return dynamic_cast<const VarDeclStmt*   >(stmt.get()) != nullptr; }
	inline bool isExprStmt(const StmtPtr& stmt) { return dynamic_cast<const ExpressionStmt*>(stmt.get()) != nullptr; }
	inline bool isReturn  (const StmtPtr& stmt) { return dynamic_cast<const ReturnStmt*    >(stmt.get()) != nullptr; }

	// use this to make a deep copy of an expression when constructing new expressions from existing ones
	ExprPtr clone(const Expression* expr);

	// expression checks
	inline bool isLiteral    (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Literal    ); }
	inline bool isVariable   (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Variable   ); }
	inline bool isBinary     (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Binary     ); }
	inline bool isUnary      (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Unary      ); }
	inline bool isCall       (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Call       ); }
	inline bool isMember     (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Member     ); }
	inline bool isInitializer(const ExprPtr& expr) { return (expr->exprType == ExpressionType::Initializer); }
	inline bool isIndex      (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Index      ); }

	inline bool isInt(const Expression* expr, int& value) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr)) {
			if (isInt(literal->value)) {
				value = literal->value.i;
				return true;
			}
		}
		return false;
	}

	inline bool isDouble(const Expression* expr, double& value) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr)) {
			if (isDouble(literal->value)) {
				value = literal->value.d;
				return true;
			}
		}
		return false;
	}

	inline bool isScalar(const Expression* expr, double& value) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr)) {
			if (isInt(literal->value)) {
				value = static_cast<double>(literal->value.i);
				return true;
			}
			if (isDouble(literal->value)) {
				value = literal->value.d;
				return true;
			}
		}
		return false;
	}

	inline bool isVec2(const Expression* expr, vec2& value) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr)) {
			if (isVec2(literal->value)) {
				value = literal->value.vec2Value;
				return true;
			}
		}
		return false;
	}

	inline bool isVec3(const Expression* expr, vec3& value) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr)) {
			if (isVec3(literal->value)) {
				value = literal->value.vec3Value;
				return true;
			}
		}
		return false;
	}

	inline bool isMat2(const Expression* expr, mat2& value) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr)) {
			if (isMat2(literal->value)) {
				value = literal->value.mat2Value;
				return true;
			}
		}
		return false;
	}

	inline bool isMat3(const Expression* expr, mat3& value) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr)) {
			if (isMat3(literal->value)) {
				value = literal->value.mat3Value;
				return true;
			}
		}
		return false;
	}

	inline bool isAdd(const Expression* expr, Expression*& left, Expression*& right) {
		if (expr->exprType == ExpressionType::Binary) {
			auto binary = dynamic_cast<const BinaryExpr*>(expr);
			if (binary->op == BinaryOp::Plus) {
				left  = binary->left.get();
				right = binary->right.get();
				return true;
			}
		}
		return false;
	}

	inline bool isSub(const Expression* expr, Expression*& left, Expression*& right) {
		if (expr->exprType == ExpressionType::Binary) {
			auto binary = dynamic_cast<const BinaryExpr*>(expr);
			if (binary->op == BinaryOp::Minus) {
				left  = binary->left.get();
				right = binary->right.get();
				return true;
			}
		}
		return false;
	}

	inline bool isMul(const Expression* expr, Expression*& left, Expression*& right) {
		if (expr->exprType == ExpressionType::Binary) {
			auto binary = dynamic_cast<const BinaryExpr*>(expr);
			if (binary->op == BinaryOp::Multiply) {
				left = binary->left.get();
				right = binary->right.get();
				return true;
			}
		}
		return false;
	}

	inline bool isDiv(const Expression* expr, Expression*& left, Expression*& right) {
		if (expr->exprType == ExpressionType::Binary) {
			auto binary = dynamic_cast<const BinaryExpr*>(expr);
			if (binary->op == BinaryOp::Divide) {
				left = binary->left.get();
				right = binary->right.get();
				return true;
			}
		}
		return false;
	}

	inline bool isExp(const Expression* expr, Expression*& left, Expression*& right) {
		if (expr->exprType == ExpressionType::Binary) {
			auto binary = dynamic_cast<const BinaryExpr*>(expr);
			if (binary->op == BinaryOp::Exponent) {
				left = binary->left.get();
				right = binary->right.get();
				return true;
			}
		}
		return false;
	}

	inline bool isNegate(const Expression* expr, Expression*& operand) {
		if (expr->exprType != ExpressionType::Unary) return false;
		auto unary = dynamic_cast<const UnaryExpr*>(expr);
		if (unary->op == UnaryOp::Negate) {
			operand = unary->right.get();
			return true;
		}
		return false;
	}

	inline bool isCall1(const Expression* expr, const std::string& fncName, const Expression*& arg) {
		if (expr->exprType != ExpressionType::Call) return false;
		auto call = dynamic_cast<const CallExpr*>(expr);
		if (call->arguments.size() != 1) return false;
		if (fncName != call->name) return false;
		arg = call->arguments[0].get();
		return true;
	}

	inline bool isCall2(const Expression* expr, const std::string& fncName, const Expression*& arg1, const Expression*& arg2) {
		if (expr->exprType != ExpressionType::Call) return false;
		auto call = dynamic_cast<const CallExpr*>(expr);
		if (call->arguments.size() != 2) return false;
		if (fncName != call->name) return false;
		arg1 = call->arguments[0].get();
		arg2 = call->arguments[1].get();
		return true;
	}

	inline bool isScalar(const Expression* expr) {
		return isScalarType(expr->valType);
	}

	inline bool isScalar(const ExprPtr& expr) {
		return isScalar(expr.get());
	}

	inline bool isZero(const Expression* expr) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr))
			return isZero(literal->value);
		if (auto init = dynamic_cast<const InitExpr*>(expr))
			return std::all_of(init->elements.begin(), init->elements.end(), [](const ExprPtr& arg) { return isZero(arg.get()); });
		if (auto ctor = dynamic_cast<const ConstructorExpr*>(expr))
			return std::all_of(ctor->args.begin(), ctor->args.end(), [](const ExprPtr& arg) { return isZero(arg.get()); });
		return false;
	}

	inline bool isZero(const ExprPtr& expr) {
		return isZero(expr.get());
	}

	inline bool isOne(const Expression* expr) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr))
			return isOne(literal->value);
		return false;
	}

	inline bool isOne(const ExprPtr& expr) {
		return isOne(expr.get());
	}

	inline bool isIdentity(const Expression* expr) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr))
		{
			const Value& v = literal->value;
			if (isMat2(v) && v.mat2Value == mat2(1.0)) return true;
			if (isMat3(v) && v.mat3Value == mat3(1.0)) return true;
		}
		return false;
	}

	inline bool isVector(const Expression* expr) {
		return (isVec2Type(expr->valType) || isVec3Type(expr->valType));
	}

	inline bool isMatrix(const Expression* expr) {
		return (isMat2Type(expr->valType) || isMat3Type(expr->valType));
	}

	inline bool isSymmetric(const Expression* expr) {
		if (auto literal = dynamic_cast<const LiteralExpr*>(expr))
		{
			const Value& v = literal->value;
			if (isMat2(v) && isSymmetric(v.mat2Value)) return true;
			if (isMat3(v) && isSymmetric(v.mat3Value)) return true;
		}
		return false;
	}

	inline bool isNegation(const ExprPtr& expr) {
		if (expr->exprType != ExpressionType::Unary) return false;
		auto unary = dynamic_cast<UnaryExpr*>(expr.get());
		return (unary->op == UnaryOp::Negate);
	}

	bool isEqual(const Expression* l, const Expression* r);
	bool isEqual(const ExprPtr& l, const ExprPtr& r);
	bool isEqual(const std::vector<ExprPtr>& l, const std::vector<ExprPtr>& r);

} // namespace febcode

std::ostream& operator << (std::ostream& o, const febcode::Value& v);

std::string ValueToString(const febcode::Value& v);

std::string ValueTypeToString(const febcode::Value& v);

std::string opToString(febcode::BinaryOp op);
