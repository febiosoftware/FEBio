#pragma once
#include <unordered_set>
#include "ast.h"
#include "program.h"

namespace febcode {
	class Optimizer {
	public:
		Optimizer(Program& program);

		void optimize();

	private:
		bool shouldRemove(Statement* stmt);
		void updateLiveness(Expression* expr);
		bool hasSideEffects(Expression* expr);

	private:
		bool shouldRemoveReturn   (ReturnStmt*     stmt);
		bool shouldRemoveVarDecl  (VarDeclStmt*    stmt);
		bool shouldRemoveExprStmt (ExpressionStmt* stmt);
		bool shouldRemoveBlockStmt(BlockStmt*      stmt);
		bool shouldRemoveIfStmt   (IfStmt*         stmt);

	private:
		Program& prg;
		std::unordered_set<std::string> live; // variables that are currently live (used in the future)
		std::unordered_set<std::string> assignedLater; // variables that are re-assigned later
	};
}
