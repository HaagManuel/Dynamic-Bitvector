#pragma once

class TaggedPointer {
    private:
        uintptr_t val;
        static const uintptr_t mask = 1;
    public:
        TaggedPointer(void* ref = NULL, bool mark = false) {
            val = ((uintptr_t)ref & ~mask) | (mark ? 1 : 0);
        }
        void* getRef() const { return (void*)(val & ~mask); }
        bool getMark() const { return (val & mask); }
};
