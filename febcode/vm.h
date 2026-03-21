#pragma once
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include "compiler.h"
#include <assert.h>

namespace febcode
{
	class VM
	{
	public:
		enum { MAX_CALL_DEPTH = 8 };

	public:
		VM() : m_program(nullptr) {}

		void setProgram(const Program& program)
		{
			m_program = &program;
			globalCount = (int)m_program->globals.size();
			m_stack.resize(globalCount + program.maxStackSize);
			stackTop = globalCount; // stack starts after global region
		}

		VM(const Program& program)
		{
			setProgram(program);
		}

		Value run()
		{
			if (m_program == nullptr) return Value();
			callFunction(0, 0);              // call "main"
			return execute();
		}

		void setDebugMode(bool b) { m_debug = b; }

		bool stackEmpty() const { return stackTop == globalCount; }
		size_t stackSize() const { return stackTop - globalCount; }

		Value getGlobal(size_t n)
		{
			if (n >= globalCount)
				throw std::runtime_error("Invalid global index: " + std::to_string(n));
			return m_stack[n];
		}

		Value getGlobal(const std::string& name)
		{
			auto it = m_program->globals.find(name);
#ifndef NDEBUG
			if (it == m_program->globals.end())
				throw std::runtime_error("Undefined global variable: " + name);
#endif
			return m_stack[it->second.slot];
		}

		void setGlobal(int i, const Value& v)
		{
#ifndef NDEBUG
			if ((i<0) || (i >= globalCount))
				throw std::runtime_error("Invalid global index: " + std::to_string(i));
#endif

			if (isStruct(v)) m_stack[i] = copyStruct(getStruct(v));
			else if (isArray(v)) m_stack[i] = copyArray(getArray(v));
			else m_stack[i] = v;
		}

		void setGlobal(int i, std::initializer_list<Value> values)
		{
#ifndef NDEBUG
			if ((i<0) || (i >= globalCount))
				throw std::runtime_error("Invalid global index: " + std::to_string(i));
#endif

			Value& v = m_stack[i];
			if (isArray(v))
			{
				ArrayValue& arr = getArray(v);

				if (arr.size() != values.size())
					throw std::runtime_error("Array initializer has incorrect number of elements.");

				Type elemType = m_program->types.getBuiltinType(*values.begin());
				std::copy(values.begin(), values.end(), arr.elements.begin());
			}
			else if (isStruct(v))
			{
				StructValue& obj = getStruct(v);
				if (obj.fields.size() != values.size())
					throw std::runtime_error("Array initializer has incorrect number of elements.");

				for (size_t i = 0; i < values.size(); ++i)
				{
					Type fieldType = obj.type->fields[i].first;
					Type valType = m_program->types.getBuiltinType(values.begin()[i]);
					if (fieldType != valType)
						throw std::runtime_error("Type mismatch in struct initializer for field: " + obj.type->fields[i].second);
					obj.fields[i] = values.begin()[i];
				}
			}
			else
				throw std::runtime_error("Initializer lists can only be assigned to arrays and structs.");
		}

	private:

		struct CallFrame
		{
			int functionIndex;
			size_t ip;
			size_t base;
		};

		// ===== Execution Loop =====

		Value execute();
	
		// ===== Call Handling =====

		void callFunction(int fnIndex, int argCount);
		void callBinaryOperator(int opIndex);

		// ===== Helpers =====

		CallFrame& currentFrame()
		{
			return m_frames[frameCount - 1];
		}

		uint8_t readByte()
		{
			return m_program->code[currentFrame().ip++];
		}

		uint16_t readUint16()
		{
			uint16_t high = readByte();
			uint16_t low = readByte();
			return (high << 8) | low;
		}

		const Value& pop()
		{
#ifndef NDEBUG
			if (stackTop <= globalCount)
				throw std::runtime_error("Stack underflow.");
#endif
			return m_stack[--stackTop];
		}

		void push(const Value& v)
		{
#ifndef NDEBUG
			if (stackTop >= m_stack.size())
				throw std::runtime_error("Stack overflow.");
#endif
			m_stack[stackTop++] = v;
		}

		void pushBool(bool b)
		{
			push(b);
		}

		void pushInt(int n)
		{
			push(n);
		}

		void pushDouble(double d)
		{
			push(d);
		}

		void pushVec2(const vec2& v)
		{
			push(v);
		}

		void pushVec3(const vec3& v)
		{
			push(v);
		}

		bool popBool()
		{
			return pop().b;
		}

		int popInt()
		{
			return pop().i;
		}

		double popDouble()
		{
			return pop().d;
		}

		vec2 popVec2()
		{
			return pop().vec2Value;
		}

		vec3 popVec3()
		{
			return pop().vec3Value;
		}

		Value& peek()
		{
#ifndef NDEBUG
			if (stackTop <= globalCount)
				throw std::runtime_error("Stack underflow.");
#endif
			return m_stack[stackTop-1];
		}

		bool peekBool()
		{
			return peek().b;
		}

		int peekInt()
		{
			return peek().i;
		}

		double peekDouble()
		{
			return peek().d;
		}

		vec2 peekVec2()
		{
			return peek().vec2Value;
		}

		vec3 peekVec3()
		{
			return peek().vec3Value;
		}

		void setBool(int slot, bool b)
		{
			m_stack[slot] = b;
		}

		void setInt(int slot, int n)
		{
			m_stack[slot] = n;
		}

		void setDouble(int slot, double d)
		{
			m_stack[slot] = d;
		}

		void setVec2(int slot, const vec2& v)
		{
			m_stack[slot] = v;
		}

		void setVec3(int slot, const vec3& v)
		{
			m_stack[slot] = v;
		}

		bool isTruthy(const Value& v)
		{
			switch (v.index)
			{
				case ValueIndex::VOID  : return false;
				case ValueIndex::BOOL  : return getBool(v);
				case ValueIndex::INT   : return getInt(v) != 0;
				case ValueIndex::DOUBLE: return getDouble(v) != 0.0;
				case ValueIndex::ARRAY : return true;
				case ValueIndex::STRUCT: return true;
			}

			return false;
		}

		std::string toString(const Value& v)
		{
			switch (v.index)
			{
			case ValueIndex::VOID  : return "void";
			case ValueIndex::BOOL  : return getBool(v) ? "true" : "false";
			case ValueIndex::INT   : return std::to_string(getInt(v));
			case ValueIndex::DOUBLE: return std::to_string(getDouble(v));
			}
			return "";
		}

	private:
		const Program* m_program;
		int globalCount = 0;

		std::vector<Value> m_stack;
		size_t stackTop = 0;

		CallFrame m_frames[MAX_CALL_DEPTH];
		size_t frameCount = 0;

		bool m_debug = false;
	};

	Value runScript(const std::string& script, bool initModules = false);
}
