#include "optimizer.h"
#include <iostream>
using namespace febcode;

Optimizer::Optimizer(Program& program) : Modifier(program) {}

void Optimizer::optimize()
{
	removedStatements = 0;
	AST& ast = *prg.ast;
	shouldRemoveBlockStmt(&ast.root);
}

bool Optimizer::shouldRemove(Statement* stmt)
{
	if      (auto retStmt  = dynamic_cast<ReturnStmt*    >(stmt)) return shouldRemoveReturn   (retStmt);
	else if (auto varDecl  = dynamic_cast<VarDeclStmt*   >(stmt)) return shouldRemoveVarDecl  (varDecl);
	else if (auto exprStmt = dynamic_cast<ExpressionStmt*>(stmt)) return shouldRemoveExprStmt (exprStmt);
	else if (auto blckStmt = dynamic_cast<BlockStmt     *>(stmt)) return shouldRemoveBlockStmt(blckStmt);
	else if (auto ifStmt   = dynamic_cast<IfStmt        *>(stmt)) return shouldRemoveIfStmt   (ifStmt  );
	else if (auto fncStmt  = dynamic_cast<FunctionStmt  *>(stmt)) return shouldRemoveFncStmt  (fncStmt );
	else if (auto structStmt = dynamic_cast<StructStmt*>(stmt)) return false; // don't remove struct declarations, since they are needed for type information, even if they are not used directly.

	throw std::runtime_error("Unsupported statement type in optimizer");
}

bool Optimizer::shouldRemoveReturn(ReturnStmt* stmt)
{
	updateLiveness(stmt->value.get());
	// don't remove return statements.
	return false;
}

bool Optimizer::shouldRemoveVarDecl(VarDeclStmt* stmt)
{
	// remove unused variables.
	for (int i = (int)stmt->vars.size() - 1; i >= 0; i--)
	{
		auto& var = stmt->vars[i];
		if (live.find(var->name) == live.end())
		{
			if (assignedLater.find(var->name) != assignedLater.end())
			{
				// the variable is assigned later, so we can't remove the declaration, 
				// but we can remove the initializer if it exists.
				assignedLater.erase(var->name);
				var->initializer.reset();
			}
			else
				stmt->vars.erase(stmt->vars.begin() + i);
		}
		else
		{
			live.erase(var->name);
			updateLiveness(var->initializer.get());
		}
	}
	return stmt->vars.empty();
}

bool Optimizer::shouldRemoveExprStmt(ExpressionStmt* stmt)
{
	ExprPtr& expr = stmt->expr;

	if (expr->exprType == ExpressionType::Assignment)
	{
		const auto* assignExpr = static_cast<const AssignExpr*>(expr.get());

		if (assignExpr->target->exprType == ExpressionType::Variable)
		{
			const auto* varExpr = static_cast<const VariableExpr*>(assignExpr->target.get());
			if (live.find(varExpr->name) == live.end())
			{
				// the variable is not live, so we can remove the assignment, but
				// only if the value being assigned doesn't have side effects, otherwise we need to keep the assignment for the side effects.
				if (hasSideEffects(assignExpr->value.get()))
				{
					assignedLater.insert(varExpr->name);
					updateLiveness(assignExpr->value.get());
					return false; // keep the assignment for the side effects, even though the variable is not live.
				}
				updateLiveness(assignExpr->value.get());
				return !hasSideEffects(assignExpr->value.get()); // remove this statement if it doesn't have side effects.
			}
			else
			{
				assignedLater.insert(varExpr->name);
				live.erase(varExpr->name);
				updateLiveness(assignExpr->value.get());
				return false;
			}
		}
		if (assignExpr->target->exprType == ExpressionType::Member)
		{
			// for member assignments, we need to check if the object is live, since the assignment has side effects on the object.
			const auto* memberExpr = static_cast<const MemberExpr*>(assignExpr->target.get());
			updateLiveness(memberExpr->object.get());
			updateLiveness(assignExpr->value.get());
			return false; // don't remove member assignments, since they have side effects on the object.
		}
	}

	updateLiveness(expr.get());
	return !hasSideEffects(expr.get());
}

bool Optimizer::shouldRemoveBlockStmt(BlockStmt* blckStmt)
{
	for (int i = (int)blckStmt->statements.size() - 1; i >= 0; --i)
	{
		Statement* stmt = blckStmt->statements[i].get();

		if (shouldRemove(stmt))
		{
			if (log)
			{
				std::string stmtType = febcode::StatementTypeToString(stmt->stmtType);

				// log the removed statement with its source location
				SourceLocation loc = stmt->location;
				*log << "Removed " << stmtType << " at line " << loc.line << "\n";
				removedStatements++;
			}
			blckStmt->statements.erase(blckStmt->statements.begin() + i);
		}
	}
	return blckStmt->statements.empty();
}

bool Optimizer::shouldRemoveIfStmt(IfStmt* stmt)
{
	auto liveBefore = live; // save live variables before processing branches

	bool removeThen = shouldRemove(stmt->thenBranch.get());
	auto liveThen = live;

	bool removeElse = true;
	if (stmt->elseBranch)
	{
		live = liveBefore; // reset live variables before processing else branch
		removeElse = shouldRemove(stmt->elseBranch.get());
	}

	if (removeThen && removeElse && !hasSideEffects(stmt->condition.get()))
	{
		live = liveBefore; // restore live variables before processing branches
		return true; // remove the entire if statement
	}
	else
	{
		// merge live variables from both branches
		for (const auto& var : liveThen)
			live.insert(var);

		updateLiveness(stmt->condition.get());

		// remove else branch if not needed
		if (removeElse)
			stmt->elseBranch.reset();

		return false;
	}
}

bool Optimizer::shouldRemoveFncStmt(FunctionStmt* stmt)
{
	auto liveBefore = live; // save live variables before processing function body
	auto assignedBefore = assignedLater;
	live.clear();
	assignedLater.clear();
	shouldRemove(stmt->body.get());
	live = liveBefore; // restore live variables before processing function body
	assignedLater = assignedBefore; // restore assignedLater before processing function body

	return false; // don't remove function declarations. Even if they are not used directly, they don't cost anything.
}

void Optimizer::updateLiveness(Expression* expr) 
{
	if (!expr) return;
	switch (expr->exprType)
	{
	case ExpressionType::Literal:
		break;
	case ExpressionType::Variable:
	{
		const auto* varExpr = static_cast<const VariableExpr*>(expr);
		live.insert(varExpr->name);
		break;
	}
	case ExpressionType::Member:
	{
		const auto* memberExpr = static_cast<const MemberExpr*>(expr);
		updateLiveness(memberExpr->object.get());
		break;
	}
	case ExpressionType::Index:
	{
		const auto* indexExpr = static_cast<const IndexExpr*>(expr);
		updateLiveness(indexExpr->object.get());
		updateLiveness(indexExpr->index.get());
		break;
	}
	case ExpressionType::Binary:
	{
		const auto* binaryExpr = static_cast<const BinaryExpr*>(expr);
		updateLiveness(binaryExpr->left.get());
		updateLiveness(binaryExpr->right.get());
		break;
	}
	case ExpressionType::Unary:
	{
		const auto* unaryExpr = static_cast<const UnaryExpr*>(expr);
		updateLiveness(unaryExpr->right.get());
		break;
	}
	case ExpressionType::Call:
	{
		const auto* callExpr = static_cast<const CallExpr*>(expr);
		for (const auto& arg : callExpr->arguments)
			updateLiveness(arg.get());
		break;
	}
	case ExpressionType::Constructor:
	{
		const auto* ctorExpr = static_cast<const ConstructorExpr*>(expr);
		for (const auto& arg : ctorExpr->args)
			updateLiveness(arg.get());
		break;
	}
	case ExpressionType::Initializer:
	{
		const auto* initExpr = static_cast<const InitExpr*>(expr);
		for (const auto& element : initExpr->elements)
			updateLiveness(element.get());
		break;
	}
	default:
		throw std::runtime_error("Unsupported expression type in optimizer");
		break; // For now, we won't handle these expression types
	}
}


bool Optimizer::hasSideEffects(Expression* expr)
{
	if (!expr) return false;
	switch (expr->exprType)
	{
	case ExpressionType::Literal:
	case ExpressionType::Variable:
	case ExpressionType::Member:
		return false;
	case ExpressionType::Index:
	{
		const auto* indexExpr = static_cast<const IndexExpr*>(expr);
		return hasSideEffects(indexExpr->object.get()) || hasSideEffects(indexExpr->index.get());
	}
	case ExpressionType::Binary:
	{
		const auto* binaryExpr = static_cast<const BinaryExpr*>(expr);
		return hasSideEffects(binaryExpr->left.get()) || hasSideEffects(binaryExpr->right.get());
	}
	case ExpressionType::Unary:
	{
		const auto* unaryExpr = static_cast<const UnaryExpr*>(expr);
		return hasSideEffects(unaryExpr->right.get());
	}
	case ExpressionType::Initializer:
	{
		const auto* initExpr = static_cast<const InitExpr*>(expr);
		for (const auto& element : initExpr->elements)
		{
			if (hasSideEffects(element.get()))
				return true;
		}
		return false;
	}
	case ExpressionType::Constructor:
	{
		const auto* ctorExpr = static_cast<const ConstructorExpr*>(expr);
		for (const auto& arg : ctorExpr->args)
		{
			if (hasSideEffects(arg.get()))
				return true;
		}
		return false;
	}
	case ExpressionType::Assignment:
		return true; // assignments have side effects
	case ExpressionType::Call:
		return true; // assume all function calls have side effects for simplicity
	default:
		throw std::runtime_error("Unsupported expression type in optimizer");
	}
}
