#include "buddy.h"
#include <iostream>
#include <cassert>
#include <stdexcept>

struct TestObj {
    static int construct_count;
    static int destruct_count;
    int value;
    
    TestObj(int v = 0) : value(v) {
        ++construct_count;
        std::cout << "[TestObj] Constructed with value " << value << std::endl;
    }
    
    ~TestObj() {
        ++destruct_count;
        std::cout << "[TestObj] Destructed, value was " << value << std::endl;
    }
};

int TestObj::construct_count = 0;
int TestObj::destruct_count = 0;

void reset_counters() {
    TestObj::construct_count = 0;
    TestObj::destruct_count = 0;
}

int main() {
    bool all_passed = true;
    
    try {
        std::cout << "=== Test 1: Single object allocation and access ===" << std::endl;
        reset_counters();
        {
            auto handle = allocate<TestObj>(BuddyAllocator(1024, 16), 42);
            if (!handle.valid()) {
                std::cerr << "Failed: Invalid handle after allocation" << std::endl;
                all_passed = false;
            } else {
                assert(TestObj::construct_count == 1);
                assert(handle->value == 42);
                handle->value = 100;
                assert(handle->value == 100);
                std::cout << "Pass: Allocation, access, and modification successful" << std::endl;
            }
        }
        assert(TestObj::destruct_count == 1);
        std::cout << "Pass: Destructor called and memory freed" << std::endl;
    } catch (...) {
        std::cerr << "Test 1 failed" << std::endl;
        all_passed = false;
    }
    
    try {
        std::cout << "\n=== Test 2: Move semantics ===" << std::endl;
        reset_counters();
        {
            auto handle1 = allocate<TestObj>(BuddyAllocator(1024, 16), 10);
            assert(handle1.valid());
            assert(TestObj::construct_count == 1);
            
            auto handle2 = std::move(handle1);
            assert(!handle1.valid());
            assert(handle2.valid());
            assert(TestObj::construct_count == 1);  // No new construction on move
            
            handle2->value = 20;
            assert(handle2->value == 20);
        }
        assert(TestObj::destruct_count == 1);  // Only one destruction
        std::cout << "Pass: Move transfers ownership, original invalid, no extra ctor/dtor" << std::endl;
    } catch (...) {
        std::cerr << "Test 2 failed" << std::endl;
        all_passed = false;
    }
    
    try {
        std::cout << "\n=== Test 3: Invalid handle dereference throws ===" << std::endl;
        reset_counters();
        {
            auto handle1 = allocate<TestObj>(BuddyAllocator(1024, 16));
            auto handle2 = std::move(handle1);
            assert(!handle1.valid());
            try {
                *handle1;
                std::cerr << "Failed: No throw on invalid dereference" << std::endl;
                all_passed = false;
            } catch (const std::logic_error&) {
                std::cout << "Pass: Throws logic_error on invalid dereference" << std::endl;
            }
        }
    } catch (...) {
        std::cerr << "Test 3 failed" << std::endl;
        all_passed = false;
    }
    
    try {
        std::cout << "\n=== Test 4: Array allocation and access ===" << std::endl;
        reset_counters();
        {
            auto arr_handle = allocate<TestObj[3]>(BuddyAllocator(1024, 16));
            if (!arr_handle.valid()) {
                std::cerr << "Failed: Invalid array handle after allocation" << std::endl;
                all_passed = false;
            } else {
                assert(TestObj::construct_count == 3);  // Default ctor for each element
                TestObj (&arr)[3] = *arr_handle;
                arr[0].value = 1;
                arr[1].value = 2;
                arr[2].value = 3;
                assert(arr[0].value == 1);
                assert(arr[1].value == 2);
                assert(arr[2].value == 3);
                std::cout << "Pass: Array allocation, access, and modification successful" << std::endl;
            }
        }
        assert(TestObj::destruct_count == 3);  // Dtor for each element
        std::cout << "Pass: Array destructor called for each element" << std::endl;
    } catch (...) {
        std::cerr << "Test 4 failed" << std::endl;
        all_passed = false;
    }
    
    try {
        std::cout << "\n=== Test 5: Ownership and lifetime management ===" << std::endl;
        reset_counters();
        {
            auto handle = allocate<TestObj>(BuddyAllocator(512, 16), 50);
            assert(handle.valid());
            
            handle->value = 60;
            assert(handle->value == 60);
            std::cout << "Pass: Allocation, access, and modification with owned allocator successful" << std::endl;
        }
        assert(TestObj::destruct_count == 1);
        std::cout << "Pass: Handle destructor properly manages object and allocator lifetime" << std::endl;
    } catch (...) {
        std::cerr << "Test 5 failed" << std::endl;
        all_passed = false;
    }
    
    if (all_passed) {
        std::cout << "\nAll tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << "\nSome tests failed!" << std::endl;
        return 1;
    }
}