#pragma once
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>
#include <cstring>
#include "compiler.h"
#include <assert.h>

namespace febcode
{
	class VM
	{
		struct Ref {
			double* ptr;
			Ref() : ptr(nullptr) {}
		};

	public:
		enum { MAX_CALL_DEPTH = 8 };

	public:
		VM();

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
		void setDebugOutput(std::ostream& os) { m_debugOutput = &os; }

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

		void setGlobal(int i, double v)
		{
#ifndef NDEBUG
			if ((i<0) || (i >= globalStackSize))
				throw std::runtime_error("Invalid global index: " + std::to_string(i));
#endif
			const Program::Global& glob = m_program->globals[i];
			int slot = glob.slot;
			setDoubleAt(slot, v);
		}

		void setGlobal(int i, vec2 v)
		{
#ifndef NDEBUG
			if ((i < 0) || (i >= globalStackSize))
				throw std::runtime_error("Invalid global index: " + std::to_string(i));
#endif
			const Program::Global& glob = m_program->globals[i];
			int slot = glob.slot;
			setVec2At(slot, v);
		}

		void setGlobal(int i, vec3 v)
		{
#ifndef NDEBUG
			if ((i<0) || (i >= globalStackSize))
				throw std::runtime_error("Invalid global index: " + std::to_string(i));
#endif
			const Program::Global& glob = m_program->globals[i];
			int slot = glob.slot;
			setVec3At(slot, v);
		}

		void setGlobal(int i, const Value& value)
		{
#ifndef NDEBUG
			if ((i < 0) || (i >= globalStackSize))
				throw std::runtime_error("Invalid global index: " + std::to_string(i));
#endif
			const Program::Global& glob = m_program->globals[i];
			int slot = glob.slot;

			switch (glob.type->kind)
			{
			case TypeKind::Bool: setBoolAt(slot, value.b); break;
			case TypeKind::Int: setIntAt(slot, value.i); break;
			case TypeKind::Double: setDoubleAt(slot, value.d); break;
			case TypeKind::Vec2: setVec2At(slot, value.vec2Value); break;
			case TypeKind::Vec3: setVec3At(slot, value.vec3Value); break;
			default:
				throw std::runtime_error("Unsupported global variable type for setGlobal.");
				break;
			}
		}

		void setInput(const std::string& name, const Value& value)
		{
			const auto it = m_program->inputIndices.find(name);

#ifndef NDEBUG
			if (it == m_program->inputIndices.end())
				throw std::runtime_error("Input variable '" + name + "' is not defined.");
#endif

			size_t index = it->second;

#ifndef NDEBUG
			if (m_program->inputs[index].type != m_program->types.getBuiltinType(value))
				throw std::runtime_error("Type mismatch when setting input variable '" + name + "'.");
#endif

			int globalSlot = m_program->inputs[index].slot;
			setGlobal(globalSlot, value);
		}

	private:
		struct CallFrame
		{
			int functionIndex;
			const uint8_t* ip;
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
			return *ip++;
		}

		uint16_t readUint16()
		{
			uint16_t high = readByte();
			uint16_t low = readByte();
			return (high << 8) | low;
		}

		void push(const double& v)
		{
#ifndef NDEBUG
			if (stackTop >= m_stack.size())
				throw std::runtime_error("Stack overflow.");
#endif
			m_stack[stackTop++] = v;
		}

		void push(double* d, size_t n)
		{
#ifndef NDEBUG
			if (stackTop + n > m_stack.size())
				throw std::runtime_error("Stack overflow.");
#endif
			memcpy(&m_stack[stackTop], d, n * sizeof(double));
			stackTop += n;
		}

		void pushVoid()
		{
			push(0.0);
		}

		void pushBool(bool b)
		{
			push((double)b);
		}

		void pushInt(int n)
		{
			push((double)n);
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

		void pushMat2(const mat2& m)
		{
			push(m(0,0));
			push(m(0,1));
			push(m(1,0));
			push(m(1,1));
		}

		void pushMat3(const mat3& m)
		{
			memcpy(&m_stack[stackTop], &(m(0,0)), 9 * sizeof(double));
			stackTop += 9;
		}

		void pushFrom(int slot, size_t size)
		{
			memcpy(&m_stack[stackTop], &m_stack[slot], size * sizeof(double));
			stackTop += size;
		}

		void pushArray(const ArrayValue& arr)
		{
			for (int i = 0; i < arr.size(); ++i)
			{
				const Value& v = arr.elements[i];
				switch (arr.type->elementType->kind)
				{
				case TypeKind::Bool  : pushBool  (v.b); break;
				case TypeKind::Int   : pushInt   (v.i); break;
				case TypeKind::Double: pushDouble(v.d); break;
				case TypeKind::Vec2  : pushVec2  (v.vec2Value); break;
				case TypeKind::Vec3  : pushVec3  (v.vec3Value); break;
				case TypeKind::Array : pushArray (*v.arrayValue); break;
				case TypeKind::Struct: pushStruct(*v.structValue); break;
				default:
					throw std::runtime_error("Unsupported array element type for push.");
				};
			}
		}

		void pushStruct(const StructValue& obj)
		{
			for (int i=0; i<obj.fields.size(); ++i)
			{
				const Value& field = obj.fields[i];
				switch (obj.type->fields[i].first->kind)
				{
				case TypeKind::Bool  : pushBool  (field.b); break;
				case TypeKind::Int   : pushInt   (field.i); break;
				case TypeKind::Double: pushDouble(field.d); break;
				case TypeKind::Vec2  : pushVec2  (field.vec2Value); break;
				case TypeKind::Vec3  : pushVec3  (field.vec3Value); break;
				case TypeKind::Array : pushArray (*field.arrayValue); break;
				case TypeKind::Struct: pushStruct(*field.structValue); break;
				default:
					throw std::runtime_error("Unsupported field type in pushStruct.");
				}
			}
		}

		const double& pop()
		{
#ifndef NDEBUG
			if (stackTop <= globalStackSize)
				throw std::runtime_error("Stack underflow.");
#endif
			return m_stack[--stackTop];
		}

		void pop(size_t n)
		{
#ifndef NDEBUG
			if (stackTop < globalStackSize + n)
				throw std::runtime_error("Stack underflow.");
#endif
			stackTop -= n;
		}

		void popVoid()
		{
			pop();
		}

		bool popBool()
		{
			return (pop() != 0.0);
		}

		int popInt()
		{
			return (int)pop();
		}

		double popDouble()
		{
			return pop();
		}

		vec2& popVec2()
		{
			stackTop -= 2;
			return *(vec2*)(&m_stack[stackTop]);
		}

		vec3& popVec3()
		{
			stackTop -= 3;
			return *(vec3*)(&m_stack[stackTop]);
		}

		mat2& popMat2()
		{
			stackTop -= 4;
			return *(mat2*)(&m_stack[stackTop]);
		}

		mat3& popMat3()
		{
			stackTop -= 9;
			return *(mat3*)(&m_stack[stackTop]);
		}

		double* popPtr(size_t n)
		{
#ifndef NDEBUG
			if (stackTop < globalStackSize + n)
				throw std::runtime_error("Stack underflow.");
#endif
			stackTop -= n;
			return &m_stack[stackTop];
		}

		mat3& getMat3At(int slot)
		{
			return *(mat3*)(&m_stack[slot]);
		}

		ArrayValue popArray(Type type)
		{
			assert(type->kind == TypeKind::Array);
			ArrayValue arr;
			arr.type = type;
			arr.elements.resize(type->arraySize);
			for (int i = (int)arr.size() - 1; i >= 0; --i)
			{
				switch (type->elementType->kind)
				{
				case TypeKind::Bool  : arr.elements[i] = popBool(); break;
				case TypeKind::Int   : arr.elements[i] = popInt(); break;
				case TypeKind::Double: arr.elements[i] = popDouble(); break;
				case TypeKind::Vec2  : arr.elements[i] = popVec2(); break;
				case TypeKind::Vec3  : arr.elements[i] = popVec3(); break;
				case TypeKind::Array : arr.elements[i] = std::make_shared<ArrayValue >(popArray(type->elementType)); break;
				case TypeKind::Struct: arr.elements[i] = std::make_shared<StructValue>(popStruct(type->elementType)); break;
				default:
					throw std::runtime_error("Unsupported array element type for pop.");
				};
			}
			return arr;
		}

		StructValue popStruct(Type type)
		{
			assert(type->kind == TypeKind::Struct);
			StructValue obj;
			obj.type = type;
			obj.fields.resize(type->fields.size());
			for (int i = (int)type->fields.size() - 1; i >= 0; --i)
			{
				Value& field = obj.fields[i];
				switch (type->fields[i].first->kind)
				{
				case TypeKind::Bool  : field = popBool(); break;
				case TypeKind::Int   : field = popInt(); break;
				case TypeKind::Double: field = popDouble(); break;
				case TypeKind::Vec2  : field = popVec2(); break;
				case TypeKind::Vec3  : field = popVec3(); break;
				case TypeKind::Array : field = std::make_shared<ArrayValue >(popArray (type->fields[i].first)); break;
				case TypeKind::Struct: field = std::make_shared<StructValue>(popStruct(type->fields[i].first)); break;
				default:
					throw std::runtime_error("Unsupported array element type for pop.");
				};
			}
			return obj;
		}

		const double& peek()
		{
#ifndef NDEBUG
			if (stackTop <= globalStackSize)
				throw std::runtime_error("Stack underflow.");
#endif
			return m_stack[stackTop-1];
		}

		double* peekPtr(size_t n)
		{
#ifndef NDEBUG
			if (stackTop < globalStackSize + n)
				throw std::runtime_error("Stack underflow.");
#endif
			return &m_stack[stackTop - n];
		}

		bool peekBool()
		{
			return (peek() != 0.0);
		}

		int peekInt()
		{
			return (int)peek();
		}

		vec2& peekVec2()
		{
			return *(vec2*)&m_stack[stackTop - 2];
		}

		vec3& peekVec3()
		{
			return *(vec3*)&m_stack[stackTop - 3];
		}

		mat2& peekMat2()
		{
			return *(mat2*)&m_stack[stackTop - 4];
		}

		mat3& peekMat3()
		{
			return *(mat3*)&m_stack[stackTop - 9];
		}

		ArrayValue peekArray(Type type)
		{
			assert(type->kind == TypeKind::Array);
			int c = (int)(stackTop - type->size());
			ArrayValue arr = getArrayAt(c, type);
			return arr;
		}

		StructValue peekStruct(Type type)
		{
			assert(type->kind == TypeKind::Struct);
			int c = (int)(stackTop - type->size());
			StructValue obj = getStructAt(c, type);
			return obj;
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

		void setMat2At(int slot)
		{
			double* s = &m_stack[stackTop - 4]; // last 4 slots
			m_stack[slot  ] = s[0];
			m_stack[slot+1] = s[1];
			m_stack[slot+2] = s[2];
			m_stack[slot+3] = s[3];
		}

		void setMat3At(int slot)
		{
			double* s = &m_stack[stackTop - 9]; // last 9 slots
			memcpy(&m_stack[slot], s, 9 * sizeof(double));
		}

		void setArrayAt(int slot, const ArrayValue& arr)
		{
			for (int i = 0; i < arr.size(); ++i)
			{
				const Value& v = arr.elements[i];
				switch (arr.type->elementType->kind)
				{
				case TypeKind::Bool  : setBoolAt  (slot, v.b); break;
				case TypeKind::Int   : setIntAt   (slot, v.i); break;
				case TypeKind::Double: setDoubleAt(slot, v.d); break;
				case TypeKind::Vec2  : setVec2At  (slot, v.vec2Value); break;
				case TypeKind::Vec3  : setVec3At  (slot, v.vec3Value); break;
				default:
					throw std::runtime_error("Unsupported array element type for setArrayAt.");
				};
				slot += (int)arr.type->elementType->size();
			}
		}

		void setStructAt(int slot, const StructValue& obj)
		{
			for (int i = 0; i < obj.fields.size(); ++i)
			{
				const Value& v = obj.fields[i];
				switch (obj.type->fields[i].first->kind)
				{
				case TypeKind::Bool  : setBoolAt  (slot, v.b); break;
				case TypeKind::Int   : setIntAt   (slot, v.i); break;
				case TypeKind::Double: setDoubleAt(slot, v.d); break;
				case TypeKind::Vec2  : setVec2At  (slot, v.vec2Value); break;
				case TypeKind::Vec3  : setVec3At  (slot, v.vec3Value); break;
				default:
					throw std::runtime_error("Unsupported array element type for setArrayAt.");
				};
				slot += (int)obj.type->fields[i].first->size();
			}
		}

		double* getPtrAt(int slot)
		{
			return &m_stack[slot];
		}

		bool getBoolAt(int slot)
		{
			return (m_stack[slot] != 0.0);
		}

		int getIntAt(int slot)
		{
			return (int)m_stack[slot];
		}

		double getDoubleAt(int slot)
		{
			return m_stack[slot];
		}

		vec2 getVec2At(int slot)
		{
			return vec2(m_stack[slot], m_stack[slot + 1]);
		}

		vec3 getVec3At(int slot)
		{
			return vec3(m_stack[slot], m_stack[slot + 1], m_stack[slot + 2]);
		}

		ArrayValue getArrayAt(int slot, Type type)
		{
			assert(type->kind == TypeKind::Array);
			ArrayValue arr;
			arr.type = type;
			arr.elements.resize(type->arraySize);

			int c = slot;
			for (int i = 0; i < arr.size(); ++i)
			{
				Value& elem = arr.elements[i];
				switch (type->elementType->kind)
				{
				case TypeKind::Bool  : elem = getBoolAt  (c); break;
				case TypeKind::Int   : elem = getIntAt   (c); break;
				case TypeKind::Double: elem = getDoubleAt(c); break;
				case TypeKind::Vec2  : elem = getVec2At  (c); break;
				case TypeKind::Vec3  : elem = getVec3At  (c); break;
				default:
					throw std::runtime_error("Unsupported type in getArrayAt");
				}
				c += (int)type->elementType->size();
			}
			return arr;
		}

		StructValue getStructAt(int slot, Type type)
		{
			assert(type->kind == TypeKind::Struct);
			StructValue obj;
			obj.type = type;
			obj.fields.resize(type->fields.size());
			int c = slot;
			for (int i = 0; i < type->fields.size(); ++i)
			{
				Value& field = obj.fields[i];
				switch (type->fields[i].first->kind)
				{
				case TypeKind::Bool  : field = getBoolAt  (c); break;
				case TypeKind::Int   : field = getIntAt   (c); break;
				case TypeKind::Double: field = getDoubleAt(c); break;
				case TypeKind::Vec2  : field = getVec2At  (c); break;
				case TypeKind::Vec3  : field = getVec3At  (c); break;
				case TypeKind::Array : field = std::make_shared<ArrayValue >(getArrayAt (c, type->fields[i].first)); break;
				case TypeKind::Struct: field = std::make_shared<StructValue>(getStructAt(c, type->fields[i].first)); break;
				default:
					throw std::runtime_error("Unsupported type in getStrucAt");
				}
				c += (int)type->elementType->size();
			}
			return obj;
		}

		void copy(size_t dest, size_t src, size_t size)
		{
			for (size_t i = 0; i < size; ++i)
				m_stack[dest + i] = m_stack[src + i];

			if (dest + size > stackTop)
				stackTop = dest + size;
		}

	private:
		const Program* m_program;
		size_t globalStackSize = 0;

		std::vector<double> m_stack;
		size_t stackTop = 0;
		Ref ref;

		CallFrame m_frames[MAX_CALL_DEPTH];
		size_t frameCount = 0;
		const uint8_t*  ip = nullptr; // current instruction pointer

		bool m_debug = false;
		std::ostream* m_debugOutput = nullptr;
	};

	Value runScript(const std::string& script);
}
