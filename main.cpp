#include <iostream>
#include "buddy.h"

int main() {
    try {
        // Setup: Total size 4096, min chunk size 64
        BuddyAllocator allocator(4096, 64);
        std::cout << "--- Initial State ---" << std::endl;
        allocator.print_report();

        // --- Allocation Phase ---
        std::cout << "\n--- Allocating p1 (100 -> 128 bytes) ---" << std::endl;
        void* p1 = allocator.allocate(100);
        allocator.print_report();

        std::cout << "\n--- Allocating p2 (30 -> 64 bytes) ---" << std::endl;
        void* p2 = allocator.allocate(30);
        allocator.print_report();

        std::cout << "\n--- Allocating p3 (250 -> 256 bytes) ---" << std::endl;
        void* p3 = allocator.allocate(250);
        allocator.print_report();

        std::cout << "\n--- Allocating p4 (500 -> 512 bytes) ---" << std::endl;
        void* p4 = allocator.allocate(500);
        allocator.print_report();

        // --- Deallocation and Merge Phase ---
        std::cout << "\n--- Deallocating p1 (128 bytes) ---" << std::endl;
        allocator.deallocate(p1);
        allocator.print_report();

        std::cout << "\n--- Deallocating p3 (256 bytes) ---" << std::endl;
        allocator.deallocate(p3);
        allocator.print_report();

        std::cout << "\n--- Deallocating p2 (64 bytes) ---" << std::endl;
        // Deallocating p2 should cause a merge with its now-free buddy (the first 64 bytes of p1's old block)
        // But p1 and p3 are not buddies, so this won't cause a large merge yet.
        allocator.deallocate(p2);
        allocator.print_report();

        std::cout << "\n--- Deallocating p4 (512 bytes) ---" << std::endl;
        // This should trigger a cascade of merges, returning the system to a single free block.
        allocator.deallocate(p4);
        std::cout << "\n--- Final State ---" << std::endl;
        allocator.print_report();

        // --- Edge Case Testing ---
        std::cout << "\n--- Testing Edge Cases ---" << std::endl;
        std::cout << "Allocating more than total size..." << std::endl;
        void* p5 = allocator.allocate(5000);
        assert(p5 == nullptr);
        std::cout << "Allocation returned nullptr as expected." << std::endl;

        std::cout << "\nDeallocating nullptr..." << std::endl;
        allocator.deallocate(nullptr);
        std::cout << "Deallocation of nullptr completed without error as expected." << std::endl;
        allocator.print_report();


    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}