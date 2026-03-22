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
			globalStackSize = m_program->globalStackSize;
			m_stack.resize(globalStackSize + program.maxStackSize);
			stackTop = globalStackSize; // stack starts after global region
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

		bool stackEmpty() const { return stackTop == globalStackSize; }
		size_t stackSize() const { return stackTop - globalStackSize; }

		Value getGlobal(size_t n)
		{
			if (n >= globalStackSize)
				throw std::runtime_error("Invalid global index: " + std::to_string(n));

			const Program::Global& glob = m_program->globals[n];
			int slot = glob.slot;
			switch (glob.type->kind)
			{
			case TypeKind::Bool  : return getBoolAt  (slot);
			case TypeKind::Int   : return getIntAt   (slot);
			case TypeKind::Double: return getDoubleAt(slot);
			case TypeKind::Vec2  : return getVec2At  (slot);
			case TypeKind::Vec3  : return getVec3At  (slot);
			default:
				return m_stack[slot];
			};
		}

		Value getGlobal(const std::string& name)
		{
			auto it = m_program->globalIndices.find(name);
#ifndef NDEBUG
			if (it == m_program->globalIndices.end())
				throw std::runtime_error("Undefined global variable: " + name);
#endif
			return getGlobal(it->second);
		}

		void setGlobal(int i, const Value& v)
		{
#ifndef NDEBUG
			if ((i<0) || (i >= globalStackSize))
				throw std::runtime_error("Invalid global index: " + std::to_string(i));
#endif
			const Program::Global& glob = m_program->globals[i];
			int slot = glob.slot;

			if      (isStruct(v)) m_stack[slot] = copyStruct(getStruct(v));
			else if (isArray (v)) m_stack[slot] = copyArray(getArray(v));
			else
			{
				switch (glob.type->kind)
				{
				case TypeKind::Bool  : setBoolAt  (slot, getBool  (v)); break;
				case TypeKind::Int   : setIntAt   (slot, getInt   (v)); break;
				case TypeKind::Double: setDoubleAt(slot, getDouble(v)); break;
				case TypeKind::Vec2  : setVec2At  (slot, getVec2  (v)); break;
				case TypeKind::Vec3  : setVec3At  (slot, getVec3  (v)); break;
					default:
					throw std::runtime_error("Unsupported global variable type for assignment.");
				}
			}
		}

		void setGlobal(int n, std::initializer_list<Value> values)
		{
#ifndef NDEBUG
			if ((n<0) || (n >= globalStackSize))
				throw std::runtime_error("Invalid global index: " + std::to_string(n));
#endif

			const Program::Global& glob = m_program->globals[n];
			int slot = (int)glob.slot;

			Value& v = m_stack[slot];
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

		void push(const Value& v)
		{
#ifndef NDEBUG
			if (stackTop >= m_stack.size())
				throw std::runtime_error("Stack overflow.");
#endif
			m_stack[stackTop++] = v;
		}

		void pushVoid()
		{
			push(Value());
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
			push(v.x);
			push(v.y);
		}

		void pushVec3(const vec3& v)
		{
			push(v.x);
			push(v.y);
			push(v.z);
		}

		void pushArray(const ArrayValuePtr& arr)
		{
			push(arr);
		}

		void pushStruct(const StructValuePtr& arr)
		{
			push(arr);
		}

		const Value& pop()
		{
#ifndef NDEBUG
			if (stackTop <= globalStackSize)
				throw std::runtime_error("Stack underflow.");
#endif
			return m_stack[--stackTop];
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
			double y = pop().d;
			double x = pop().d;
			return vec2(x, y);
		}

		vec3 popVec3()
		{
			double z = pop().d;
			double y = pop().d;
			double x = pop().d;
			return vec3(x, y, z);
		}

		const ArrayValue& popArray()
		{
			return getArray(pop());
		}

		const StructValue& popStruct()
		{
			return getStruct(pop());
		}

		const Ref& popRef()
		{
			return pop().ref;
		}

		Value& peek()
		{
#ifndef NDEBUG
			if (stackTop <= globalStackSize)
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
			return vec2(m_stack[stackTop - 2].d, m_stack[stackTop - 1].d);
		}

		vec3 peekVec3()
		{
			return vec3(m_stack[stackTop - 3].d, m_stack[stackTop - 2].d, m_stack[stackTop - 1].d);
		}

		void setBoolAt(int slot, bool b)
		{
			m_stack[slot] = b;
		}

		void setIntAt(int slot, int n)
		{
			m_stack[slot] = n;
		}

		void setDoubleAt(int slot, double d)
		{
			m_stack[slot] = d;
		}

		void setVec2At(int slot, const vec2& v)
		{
			m_stack[slot  ] = v.x;
			m_stack[slot+1] = v.y;
		}

		void setVec3At(int slot, const vec3& v)
		{
			m_stack[slot    ] = v.x;
			m_stack[slot + 1] = v.y;
			m_stack[slot + 2] = v.z;
		}

		bool getBoolAt(int slot)
		{
			return m_stack[slot].b;
		}

		int getIntAt(int slot)
		{
			return m_stack[slot].i;
		}

		double getDoubleAt(int slot)
		{
			return m_stack[slot].d;
		}

		vec2 getVec2At(int slot)
		{
			return vec2(m_stack[slot].d, m_stack[slot + 1].d);
		}

		vec3 getVec3At(int slot)
		{
			return vec3(m_stack[slot].d, m_stack[slot + 1].d, m_stack[slot + 2].d);
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
		size_t globalStackSize = 0;

		std::vector<Value> m_stack;
		size_t stackTop = 0;

		CallFrame m_frames[MAX_CALL_DEPTH];
		size_t frameCount = 0;

		bool m_debug = false;
	};

	Value runScript(const std::string& script, bool initModules = false);
}
