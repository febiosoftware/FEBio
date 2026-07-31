#pragma once
#include "ast.h"
#include "program.h"
#include <unordered_map>
#include <string>

namespace febcode {

	// This class is resonsible for type inference and other semantic analysis. 
	// It should be called after the AST is generated, and it will traverse the AST to resolve types and check for semantic errors. It will also generate intermediate code for the AST, which will be used later for code generation.
	class Resolver
	{
	public:
		Resolver(Program& prg);
		~Resolver();
		void resolve();

		using Scope = std::unordered_map<std::string, Type>;

	private:
		void resolveStatement   (Statement*      stmt);
		void resolveExprStmt    (ExpressionStmt* stmt);
		void resolveVarDecl     (VarDeclStmt*    stmt);
		void resolveIfStmt      (IfStmt*         stmt);
		void resolveWhileStmt   (WhileStmt*      stmt);
		void resolveForStmt     (ForStmt*        stmt);
		void resolveBlockStmt   (BlockStmt*      stmt);
		void resolveReturnStmt  (ReturnStmt*     stmt);
		void resolveFunctionStmt(FunctionStmt*   stmt);
		void resolveStructStmt  (StructStmt*     stmt);

		void resolveExpression (Expression*      expr);
		void resolveLiteral    (LiteralExpr*     expr);
		void resolveVariable   (VariableExpr*    expr);
		void resolveAssignment (AssignExpr*      expr);
		void resolveUnary      (UnaryExpr*       expr);
		void resolveBinary     (BinaryExpr*      expr);
		void resolveMember     (MemberExpr*      expr);
		void resolveIndex      (IndexExpr*       expr);
		void resolveInitializer(InitExpr*        expr);
		void resolveConstructor(ConstructorExpr* expr);
		void resolveCall       (CallExpr*        expr);

	private:
		void pushScope() {
			scopeStack.emplace_back();
		}

		void popScope() {
			scopeStack.pop_back();
		}

		void declare(const std::string& name, Type type) {
			Scope& current = scopeStack.back();
			current[name] = { type };
		}

		Type lookup(const std::string& name) {
			for (int i = (int)scopeStack.size() - 1; i >= 0; --i) {
				auto it = scopeStack[i].find(name);
				if (it != scopeStack[i].end())
					return it->second;
			}
			return nullptr;
		}

	private:
		Program& prg;

		std::vector<Scope> scopeStack;
	};
}
