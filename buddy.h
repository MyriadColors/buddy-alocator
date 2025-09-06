#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <memory>
#include <cassert>
#include <unordered_set>
#include <cmath>

// Helper to check if a number is a power of two
static inline bool is_power_of_two(size_t n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

// Helper to get the next power of two using bit twiddling
static inline size_t next_power_of_two(size_t n) {
    if (n == 0) return 1;
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    if (sizeof(size_t) == 8) n |= n >> 32; // For 64-bit systems
    return n + 1;
}

// The memory chunk metadata
struct MemChunk {
    size_t offset;
    size_t size;
    bool is_free;
    size_t level;

    MemChunk(size_t offset, size_t size, bool is_free, size_t level)
        : offset(offset), size(size), is_free(is_free), level(level) {}
};

class BuddyAllocator {
private:
    size_t total_pool_size;
    size_t min_chunk_size;
    size_t num_levels;

    // The actual memory pool
    std::vector<std::byte> pool;
    // Pointer to the start of the pool for address calculations
    std::byte* pool_start;

    // Canonical owner of all memory chunks
    std::unordered_map<size_t, std::shared_ptr<MemChunk>> chunks_by_offset;

    // Free lists, one for each size level
    std::vector<std::list<std::shared_ptr<MemChunk>>> free_lists;

    // --- Helper Methods ---

    // Calculates the free list index for a given chunk size
    ssize_t size_to_level(size_t size) const {
        size_t level = 0;
        size_t s = min_chunk_size;
        while (s < size) { s <<= 1; ++level; }
        return level;
    }

    // Calculates the chunk size for a given free list level
    size_t level_to_size(size_t level) const {
        return min_chunk_size << level;
    }

    void split_chunk(std::shared_ptr<MemChunk> chunk) {
        // 1. The parent chunk is being replaced, so remove it from the canonical map.
        chunks_by_offset.erase(chunk->offset);

        // 2. Calculate properties of the two new smaller chunks (buddies).
        size_t child_level = chunk->level - 1;
        size_t child_size = level_to_size(child_level);
        size_t offset1 = chunk->offset;
        size_t offset2 = chunk->offset + child_size;

        // 3. Create the two new buddy chunks.
        auto child1 = std::make_shared<MemChunk>(offset1, child_size, true, child_level);
        auto child2 = std::make_shared<MemChunk>(offset2, child_size, true, child_level);

        // 4. Add the new chunks to the canonical map.
        chunks_by_offset[offset1] = child1;
        chunks_by_offset[offset2] = child2;

        // 5. Add the new chunks to the front of the lower-level free list.
        // The allocate loop will pick up the first one to continue splitting.
        free_lists[child_level].push_front(child2);
        free_lists[child_level].push_front(child1);
    }

public:
    BuddyAllocator(size_t total_size, size_t min_size)
        : total_pool_size(total_size), min_chunk_size(min_size), pool_start(nullptr) {

        // 1. Validate inputs
        if (!is_power_of_two(total_pool_size) || !is_power_of_two(min_chunk_size)) {
            throw std::invalid_argument("Total size and min chunk size must be powers of two.");
        }
        if (min_chunk_size > total_pool_size) {
            throw std::invalid_argument("Min chunk size cannot be larger than total size.");
        }

        // 2. Initialize the memory pool
        pool.resize(total_pool_size);
        pool_start = pool.data();

        // 3. Calculate and resize the free lists
        num_levels = size_to_level(total_pool_size) + 1;
        free_lists.resize(num_levels);

        // 4. Create the initial root chunk that covers the entire pool
        size_t root_level = num_levels - 1;
        auto root_chunk = std::make_shared<MemChunk>(0, total_pool_size, true, root_level);

        // 5. Add the root chunk to the canonical map and the top-level free list
        chunks_by_offset[0] = root_chunk;
        free_lists[root_level].push_back(root_chunk);

        std::cout << "Buddy Allocator initialized." << std::endl;
        std::cout << "Total size: " << total_pool_size << " bytes" << std::endl;
        std::cout << "Min chunk size: " << min_chunk_size << " bytes" << std::endl;
        std::cout << "Number of levels: " << num_levels << std::endl;
    }

    // Destructor
    ~BuddyAllocator() = default;

    void* allocate(size_t request_size) {
        // 1. Determine the required chunk size
        size_t required_size = std::max(min_chunk_size, next_power_of_two(request_size));

        if (required_size > total_pool_size) {
            std::cerr << "Allocation failed: Requested size " << request_size << " is too large." << std::endl;
            return nullptr; // Request exceeds total memory
        }

        // 2. Find the level for this size
        size_t target_level = size_to_level(required_size);

        // 3. Find the first available chunk
        size_t current_level = 0;
        bool found = false;
        for (size_t i = target_level; i < num_levels; ++i) {
            if (!free_lists[i].empty()) {
                current_level = i;
                found = true;
                break;
            }
        }

        if (!found) {
            std::cerr << "Allocation failed: No free chunk large enough." << std::endl;
            return nullptr; // No suitable chunk found
        }

        // 4. Retrieve the chunk and remove it from its free list
        std::shared_ptr<MemChunk> chunk_to_allocate = free_lists[current_level].front();
        free_lists[current_level].pop_front();

        // 5. Split the chunk down to the required size
        while (current_level > target_level) {
            split_chunk(chunk_to_allocate);
            current_level--;
            chunk_to_allocate = free_lists[current_level].front();
            free_lists[current_level].pop_front();
        }

        // 6. Mark the final chunk as allocated and return the pointer
        chunk_to_allocate->is_free = false;

#ifndef NDEBUG
        check_invariants();
#endif

        return pool_start + chunk_to_allocate->offset;
    }

    void deallocate(void* ptr) {
        if (ptr == nullptr) {
            return; // Deallocating null is a no-op.
        }

        // 1. Convert pointer to offset and validate it.
        ptrdiff_t offset_diff = static_cast<std::byte*>(ptr) - pool_start;
        assert(offset_diff >= 0 && "Deallocation pointer is before pool start!");
        size_t offset = static_cast<size_t>(offset_diff);
        assert(offset < total_pool_size && "Deallocation pointer is out of bounds!");

        auto it = chunks_by_offset.find(offset);
        assert(it != chunks_by_offset.end() && "Invalid pointer passed to deallocate (chunk not found)!");

        auto chunk = it->second;
        assert(!chunk->is_free && "Double deallocation detected!");

        // 2. Mark the chunk as free.
        chunk->is_free = true;

        // 3. Iteratively merge with buddies.
        while (chunk->level < num_levels - 1) {
            size_t buddy_offset = chunk->offset ^ chunk->size;
            auto buddy_it = chunks_by_offset.find(buddy_offset);

            // Check if buddy is available and suitable for merging.
            if (buddy_it == chunks_by_offset.end() || !buddy_it->second->is_free || buddy_it->second->level != chunk->level) {
                break; // Buddy not available, stop merging.
            }

            auto buddy = buddy_it->second;

            // --- Merge is possible ---
            // a. The buddy is in a free list, remove it.
            free_lists[chunk->level].remove(buddy);

            // b. Remove the two smaller chunks from the canonical map.
            chunks_by_offset.erase(chunk->offset);
            chunks_by_offset.erase(buddy->offset);

            // c. Create the new, larger parent chunk.
            size_t parent_offset = std::min(chunk->offset, buddy->offset);
            size_t parent_level = chunk->level + 1;
            auto parent_chunk = std::make_shared<MemChunk>(parent_offset, level_to_size(parent_level), true, parent_level);

            // d. Add parent to the map and continue the merge check from the parent's level.
            chunks_by_offset[parent_offset] = parent_chunk;
            chunk = parent_chunk;
        }

        // 4. Add the final (potentially merged) chunk to its free list.
        free_lists[chunk->level].push_back(chunk);

#ifndef NDEBUG
        check_invariants();
#endif
    }

    void print_report() const {
        print_header();
        print_memory_usage();
        print_free_lists_summary();
        print_all_chunks_list();
    }

private:
    // --- Reporting Helpers ---

    void print_header() const {
        std::cout << "\n=================================================" << std::endl;
        std::cout << "        Buddy Allocator State Report" << std::endl;
        std::cout << "=================================================" << std::endl;
        std::cout << "Total Size: " << total_pool_size << " bytes" << std::endl;
        std::cout << "Min Chunk Size: " << min_chunk_size << " bytes" << std::endl;
        std::cout << "Pool Start Address: " << static_cast<const void*>(pool_start) << std::endl << std::endl;
    }

    void print_memory_usage() const {
        auto [memory_in_use, allocated_chunks_count] = get_memory_in_use();
        auto [memory_free, free_chunks_count] = get_free_memory();
        double external_fragmentation = calculate_external_fragmentation();

        std::cout << "--- Memory Usage ---" << std::endl;
        std::cout << "Used: " << memory_in_use << " bytes (" << allocated_chunks_count << " chunks)" << std::endl;
        std::cout << "Free: " << memory_free << " bytes (" << free_chunks_count << " chunks)" << std::endl;
        std::cout << "External Fragmentation (approx): " << external_fragmentation << "%" << std::endl << std::endl;
    }

    void print_free_lists_summary() const {
        std::cout << "--- Free Lists ---" << std::endl;
        for (size_t i = 0; i < num_levels; i++) {
            size_t level_size = level_to_size(i);
            std::cout << "Level " << i << " (" << level_size << " bytes): " << free_lists[i].size() << " free chunks" << std::endl;
        }
        std::cout << std::endl;
    }

    void print_all_chunks_list() const {
        auto sorted_chunks = get_sorted_chunks();
        std::cout << "--- All Chunks (" << chunks_by_offset.size() << " total) ---" << std::endl;
        for (const auto& chunk : sorted_chunks) {
            std::cout << "  "
                    << (chunk->is_free ? "[FREE]" : "[USED]")
                    << "\tAddr: " << static_cast<const void*>(pool_start + chunk->offset)
                    << "\tSize: " << chunk->size << " bytes"
                    << "\tLevel: " << chunk->level
                    << std::endl;
        }
        std::cout << "=================================================" << std::endl << std::endl;
    }

    // --- Calculation Helpers for Reporting ---

    std::pair<size_t, size_t> get_memory_in_use() const {
        size_t memory_in_use = 0;
        size_t allocated_chunks_count = 0;
        for (const auto& pair : chunks_by_offset) {
            if (!pair.second->is_free) {
                memory_in_use += pair.second->size;
                allocated_chunks_count++;
            }
        }
        return {memory_in_use, allocated_chunks_count};
    }

    std::pair<size_t, size_t> get_free_memory() const {
        size_t memory_free = 0;
        size_t free_chunks_count = 0;
#ifndef NDEBUG
        std::unordered_set<const MemChunk*> counted_chunks;
#endif
        for(const auto& list : free_lists) {
            free_chunks_count += list.size();
            for(const auto& chunk : list) {
                memory_free += chunk->size;
#ifndef NDEBUG
                assert(counted_chunks.find(chunk.get()) == counted_chunks.end() && "Double count detected in free list!");
                counted_chunks.insert(chunk.get());
#endif
            }
        }
        return {memory_free, free_chunks_count};
    }

    size_t get_largest_free_chunk_size() const {
        for (size_t i = num_levels; i-- > 0;) {
            if (!free_lists[i].empty()) {
                return level_to_size(i);
            }
        }
        return 0;
    }

    double calculate_external_fragmentation() const {
        auto [memory_free, free_chunks_count] = get_free_memory();
        if (memory_free == 0 || free_chunks_count <= 1) {
            return 0.0;
        }

        size_t largest_free_chunk = get_largest_free_chunk_size();
        if (memory_free > largest_free_chunk) {
            return 100.0 * (1.0 - static_cast<double>(largest_free_chunk) / memory_free);
        }
        return 0.0;
    }

    std::vector<std::shared_ptr<MemChunk>> get_sorted_chunks() const {
        std::vector<std::shared_ptr<MemChunk>> sorted_chunks;
        sorted_chunks.reserve(chunks_by_offset.size());
        for(const auto& pair : chunks_by_offset) {
            sorted_chunks.push_back(pair.second);
        }
        std::sort(sorted_chunks.begin(), sorted_chunks.end(), [](const auto& a, const auto& b){
            return a->offset < b->offset;
        });
        return sorted_chunks;
    }

#ifndef NDEBUG
    void check_invariants() const {
        // 1. Every chunk in a free_list must be in chunks_by_offset and marked free.
        std::unordered_set<std::shared_ptr<MemChunk>> all_free_chunks;
        for (size_t i = 0; i < num_levels; ++i) {
            for (const auto& chunk : free_lists[i]) {
                assert(chunk->level == i && "Invariant Fail: Chunk in wrong free list level.");
                
                auto it = chunks_by_offset.find(chunk->offset);
                assert(it != chunks_by_offset.end() && "Invariant Fail: Free chunk not in canonical map.");
                assert(it->second == chunk && "Invariant Fail: Free chunk in map is not the same object.");
                assert(it->second->is_free && "Invariant Fail: Chunk in free list is not marked as free.");

                all_free_chunks.insert(chunk);
            }
        }

        // 2. Every chunk in chunks_by_offset that is free must be in a free_list.
        for (const auto& pair : chunks_by_offset) {
            const auto& chunk = pair.second;
            if (chunk->is_free) {
                assert(all_free_chunks.count(chunk) && "Invariant Fail: Free chunk in map is not in any free list.");
            }
        }
    }
#endif
};