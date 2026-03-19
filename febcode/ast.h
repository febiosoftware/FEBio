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

	enum class ExpressionType {
		Literal,
		Variable,
		Binary,
		Unary,
		Call,
		Member,
		Initializer,
		Index,
		Assignment
	};

	// Base classes
	struct Expression {
		Expression(ExpressionType type) : exprType(type) {}
		virtual ~Expression() = default;
		ExpressionType exprType;
	};

	struct Statement {
		virtual ~Statement() = default;
	};

	using ExprPtr = std::unique_ptr<Expression>;

	// --- Expressions ---

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
		ExprPtr callee;
		std::vector<ExprPtr> arguments;

		CallExpr(ExprPtr callee, std::vector<ExprPtr> arguments)
			: Expression(ExpressionType::Call), callee(std::move(callee)),
			arguments(std::move(arguments)) {}
	};

	struct MemberExpr : Expression
	{
		ExprPtr object;
		std::string property;

		MemberExpr(ExprPtr object, std::string property)
			: Expression(ExpressionType::Member), object(std::move(object)), property(std::move(property)) {}
	};

	struct InitializerExpr : Expression {
		std::vector<ExprPtr> elements;
		InitializerExpr(std::vector<ExprPtr> elements)
			: Expression(ExpressionType::Initializer), elements(std::move(elements)) {}
	};

	struct IndexExpr : Expression {
		ExprPtr object;
		ExprPtr index;
		IndexExpr(ExprPtr object, ExprPtr index)
			: Expression(ExpressionType::Index), object(std::move(object)), index(std::move(index)) {}
	};

	// --- Statements ---

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
		std::vector<Var> vars;
		VarDeclStmt(Type type, const std::string& name, ExprPtr initializer) : type(type)
		{
			vars.push_back({ name, std::vector<size_t>(), std::move(initializer)});
		}
		VarDeclStmt(Type type, std::vector<Var>& vars)
			: type(type), vars(std::move(vars)) {}
	};

	struct ReturnStmt : Statement {
		ExprPtr value; // can be null
		ReturnStmt(ExprPtr v) : value(std::move(v)) {}
	};

	struct BlockStmt : Statement {
		std::vector<std::unique_ptr<Statement>> statements;
	};

	struct IfStmt : Statement {
		ExprPtr condition;
		std::unique_ptr<Statement> thenBranch;
		std::unique_ptr<Statement> elseBranch; // optional

		IfStmt(ExprPtr cond,
			std::unique_ptr<Statement> thenStmt,
			std::unique_ptr<Statement> elseStmt = nullptr)
			: condition(std::move(cond)),
			thenBranch(std::move(thenStmt)),
			elseBranch(std::move(elseStmt)) {
		}
	};

	struct WhileStmt : Statement {
		ExprPtr condition;
		std::unique_ptr<Statement> body;

		WhileStmt(ExprPtr cond,
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

		void addStatement(std::unique_ptr<Statement> stmt) {
			statements.push_back(std::move(stmt));
		}
	};

	// use this to make a deep copy of an expression when constructing new expressions from existing ones
	ExprPtr copy_expression(const Expression* expr);

	// expression checks
	inline bool isLiteral    (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Literal    ); }
	inline bool isVariable   (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Variable   ); }
	inline bool isBinary     (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Binary     ); }
	inline bool isUnary      (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Unary      ); }
	inline bool isCall       (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Call       ); }
	inline bool isMember     (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Member     ); }
	inline bool isInitializer(const ExprPtr& expr) { return (expr->exprType == ExpressionType::Initializer); }
	inline bool isIndex      (const ExprPtr& expr) { return (expr->exprType == ExpressionType::Index      ); }

	inline bool isZero(const ExprPtr& expr) {
		if (!isLiteral(expr)) return false;
		auto literal = dynamic_cast<LiteralExpr*>(expr.get());
		return isZero(literal->value);
	}

	inline bool isNegation   (const ExprPtr& expr) {
		if (expr->exprType != ExpressionType::Unary) return false;
		auto unary = dynamic_cast<UnaryExpr*>(expr.get());
		return (unary->op == UnaryOp::Negate);
	}

	// helper functions to create expressions more easily
	inline ExprPtr Literal(const Value& v) { return std::make_unique<LiteralExpr>(v); }
	inline ExprPtr Variable(const std::string& name) { return std::make_unique<VariableExpr>(name); }
	inline ExprPtr Unary(UnaryOp op, const ExprPtr& arg) { return std::make_unique<UnaryExpr>(op, std::move(copy_expression(arg.get()))); }
	inline ExprPtr Negate(const ExprPtr& arg) { return std::make_unique<UnaryExpr>(UnaryOp::Negate, std::move(copy_expression(arg.get()))); }
	inline ExprPtr Binary(const ExprPtr& left, BinaryOp op, const ExprPtr& right) { return std::make_unique<BinaryExpr>(std::move(copy_expression(left.get())), op, std::move(copy_expression(right.get()))); }
	inline ExprPtr Assign(const ExprPtr& target, const ExprPtr& value) { return std::make_unique<AssignExpr>(std::move(copy_expression(target.get())), std::move(copy_expression(value.get()))); }
	inline ExprPtr Call(const ExprPtr& callee, const std::vector<ExprPtr>& args)
	{
		std::vector<ExprPtr> copyArgs;
		for (auto& arg : args)
		{
			copyArgs.emplace_back(copy_expression(arg.get()));
		}
		return std::make_unique<CallExpr>(std::move(copy_expression(callee.get())), std::move(copyArgs));
	}

	inline ExprPtr Call(const std::string& fnc, const std::vector<ExprPtr>& args)
	{
		return Call(Variable(fnc), args);
	}

	inline ExprPtr Member(const ExprPtr& object, const std::string& property)
	{ 
		return std::make_unique<MemberExpr>(std::move(copy_expression(object.get())), property); 
	}

	inline ExprPtr Initializer(size_t n, Value v)
	{
		std::vector<ExprPtr> init(n);
		for (int i = 0; i < n; ++i) init[i] = Literal(v);
		return std::make_unique<InitializerExpr>(std::move(init));
	}

	ExprPtr Initializer(const std::vector<StructField>& fields);

	bool isEqual(const ExprPtr& l, const ExprPtr& r);
	bool isEqual(const std::vector<ExprPtr>& l, const std::vector<ExprPtr>& r);

} // namespace febcode

std::ostream& operator << (std::ostream& o, const febcode::Value& v);

std::string ValueToString(const febcode::Value& v);

std::string ValueTypeToString(const febcode::Value& v);

// Expression "algebra" for easier construction of new expressions from existing ones
inline febcode::ExprPtr operator - (const febcode::ExprPtr& expr) { return Negate(expr); }
inline febcode::ExprPtr operator + (const febcode::ExprPtr& left, const febcode::ExprPtr& right) { return Binary(left, febcode::BinaryOp::Plus, right); }
inline febcode::ExprPtr operator - (const febcode::ExprPtr& left, const febcode::ExprPtr& right) { return Binary(left, febcode::BinaryOp::Minus, right); }
inline febcode::ExprPtr operator * (const febcode::ExprPtr& left, const febcode::ExprPtr& right) { return Binary(left, febcode::BinaryOp::Multiply, right); }
inline febcode::ExprPtr operator / (const febcode::ExprPtr& left, const febcode::ExprPtr& right) { return Binary(left, febcode::BinaryOp::Divide, right); }
inline febcode::ExprPtr Pow(const febcode::ExprPtr& left, const febcode::ExprPtr& right) { return Binary(left, febcode::BinaryOp::Exponent, right); }

