#pragma once
#include <unordered_set>
#include "ast.h"
#include "modifier.h"

namespace febcode {
	class Optimizer : public Modifier {
	public:
		Optimizer(Program& program);

		void optimize();

		void setLogStream(std::ostream& logStream) { log = &logStream; }

		size_t getRemovedStatementCount() const { return removedStatements; }

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
		bool shouldRemoveFncStmt  (FunctionStmt*   stmt);

	private:
		std::unordered_set<std::string> live; // variables that are currently live (used in the future)
		std::unordered_set<std::string> assignedLater; // variables that are re-assigned later

		std::ostream* log = nullptr;
		size_t removedStatements = 0;
	};
}
