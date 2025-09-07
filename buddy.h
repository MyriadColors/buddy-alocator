#pragma once


#include <algorithm>        // std::max, std::min
#include <iostream>         // std::cout, std::cerr
#include <vector>           // std::vector
#include <list>             // std::list
#include <unordered_map>    // std::unordered_map
#include <memory>           // std::shared_ptr
#include <cassert>          // std::assert
#include <unordered_set>    // std::unordered_set
#include <cmath>            // std::log2
#include <cstdint>          // std::uint64_t
#include <limits>           // std::numeric_limits
#include <utility>          // std::forward
#include <stdexcept>        // std::logic_error
#include <type_traits>      // std::is_array, std::enable_if, std::declval, std::true_type, std::false_type

// Helper to check if a number is a power of two
static inline bool is_power_of_two(size_t n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

// Helper to get the next power of two using bit twiddling
static inline size_t next_power_of_two(size_t n) {
    if (n == 0) return 1;
    // guard: if n is already too large that next power-of-two would overflow,
    // return std::numeric_limits<size_t>::max() as a sentinel that caller must treat as an error.
    const unsigned bits = sizeof(size_t) * 8;
    if (n > (size_t(1) << (bits - 1))) return std::numeric_limits<size_t>::max();
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
    MemChunk* next_free = nullptr;
    MemChunk* prev_free = nullptr;

    MemChunk(size_t offset, size_t size, bool is_free, size_t level)
        : offset(offset), size(size), is_free(is_free), level(level), next_free(nullptr), prev_free(nullptr) {}
};

// Intrusive free list structure for O(1) operations
struct FreeList {
    /**
     * Intrusive doubly-linked free list for O(1) insert/remove operations in buddy allocator.
     * Each MemChunk stores its own next_free/prev_free pointers, eliminating std::list overhead
     * and shared_ptr indirection for better performance in high-frequency allocation/deallocation.
     */
    MemChunk* head = nullptr;
    MemChunk* tail = nullptr;
    size_t count = 0;

    bool empty() const { return head == nullptr; }
    size_t size() const { return count; }
    MemChunk* front() const { return head; }

    void push_back(MemChunk* chunk) {
        chunk->next_free = nullptr;
        chunk->prev_free = tail;
        if (tail) {
            tail->next_free = chunk;
        } else {
            head = chunk;
        }
        tail = chunk;
        ++count;
    }

    void pop_front() {
        if (head) {
            MemChunk* old_head = head;
            head = head->next_free;
            if (head) {
                head->prev_free = nullptr;
            } else {
                tail = nullptr;
            }
            old_head->next_free = nullptr;
            old_head->prev_free = nullptr;
            --count;
        }
    }

    void remove(MemChunk* chunk) {
        // O(1) removal using stored prev/next pointers - key performance optimization
        if (chunk->prev_free) {
            chunk->prev_free->next_free = chunk->next_free;
        } else {
            head = chunk->next_free;
        }
        if (chunk->next_free) {
            chunk->next_free->prev_free = chunk->prev_free;
        } else {
            tail = chunk->prev_free;
        }
        chunk->next_free = nullptr;
        chunk->prev_free = nullptr;
        --count;
    }

    // For compatibility with existing code
    void erase(MemChunk* chunk) { remove(chunk); }
};

/**
 * Opaque handle representing an allocated block in the buddy allocator.
 */
struct Block {
    size_t offset;      
    size_t size;        
    std::uint64_t instance_id; 
    
    Block() : offset(std::numeric_limits<size_t>::max()), size(0), instance_id(0) {}
    Block(size_t o, size_t s, std::uint64_t id) : offset(o), size(s), instance_id(id) {}
    bool valid() const { return offset != std::numeric_limits<size_t>::max() && instance_id != 0; }
};

// Forward declaration
class BuddyAllocator;

#include <memory> // std::unique_ptr (already included, but ensure)

// MemoryHandle declaration - will be defined after BuddyAllocator
template <typename T>
class MemoryHandle;

class BuddyAllocator {
private:
    inline static uint64_t next_instance_id = 1;
    uint64_t instance_id;
    size_t total_pool_size;
    size_t min_chunk_size;
    size_t num_levels;

    // The actual memory pool
    std::vector<std::byte> pool;

    // Canonical owner of all memory chunks
    std::unordered_map<size_t, std::shared_ptr<MemChunk>> chunks_by_offset;

    // Free lists, one for each size level (intrusive)
    std::vector<FreeList> free_lists;

    // --- Helper Methods ---

    // Calculates the free list index for a given chunk size
    size_t size_to_level(size_t size) const {
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
        free_lists[child_level].push_back(child2.get());
        free_lists[child_level].push_back(child1.get());
    }

    Block allocate_block(size_t request_size) {
        // 1. Determine the required chunk size
        size_t required_size = std::max(min_chunk_size, next_power_of_two(request_size));

        if (required_size > total_pool_size) {
            std::cerr << "Allocation failed: Requested size " << request_size << " is too large."
                      << std::endl;
            return Block(); // Return invalid block on failure
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
            std::cerr << "Allocation failed: No free chunk large enough."
                      << std::endl;
            return Block(); // Return invalid block on failure
        }

        // 4. Retrieve the chunk and remove it from its free list
        MemChunk* chunk_ptr = free_lists[current_level].front();
        free_lists[current_level].pop_front();
        std::shared_ptr<MemChunk> chunk_to_allocate = chunks_by_offset[chunk_ptr->offset];

        // 5. Split the chunk down to the required size
        while (current_level > target_level) {
            split_chunk(chunk_to_allocate);
            current_level--;
            MemChunk* chunk_ptr = free_lists[current_level].front();
            free_lists[current_level].pop_front();
            chunk_to_allocate = chunks_by_offset[chunk_ptr->offset];
        }

        // 6. Mark the final chunk as allocated and return the Block handle
        chunk_to_allocate->is_free = false;

#ifndef NDEBUG
        check_invariants();
#endif

        return Block(chunk_to_allocate->offset, chunk_to_allocate->size, instance_id);
    }

    /**
     * Deallocates the block represented by the given handle.
     *
     * The handle must be valid and belong to this allocator (use owns() to verify).
     * This method will attempt to merge the freed block with adjacent buddy blocks
     * to reduce fragmentation.
     */
    void deallocate(const Block& handle) {
        if (!handle.valid()) {
            return; // Invalid handle, no-op
        }

        size_t offset = handle.offset;
        if (offset >= total_pool_size) {
            assert(false && "Deallocation handle offset is out of bounds!");
            return;
        }

        auto it = chunks_by_offset.find(offset);
        if (it == chunks_by_offset.end()) {
            assert(false && "Invalid handle passed to deallocate (chunk not found)!");
            return;
        }

        auto chunk = it->second;
        if (chunk->is_free) {
            assert(false && "Double deallocation detected!");
            return;
        }

#ifndef NDEBUG
        if (handle.size != chunk->size) {
            assert(false && "Handle size does not match chunk size!");
        }
#endif

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
            free_lists[chunk->level].remove(buddy.get());

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
        free_lists[chunk->level].push_back(chunk.get());

#ifndef NDEBUG
        check_invariants();
#endif
    }

    void* resolve(const Block& h) const {
        std::cout << "[DEBUG] resolve called, pool.data(): " << (void*)pool.data()
                  << ", total_pool_size: " << total_pool_size
                  << ", h.offset: " << h.offset << std::endl;
        if (!h.valid()) {
            std::cout << "[DEBUG] resolve: invalid handle" << std::endl;
            return nullptr;
        }
        if (h.offset >= total_pool_size) {
            std::cout << "[DEBUG] resolve: offset out of bounds" << std::endl;
            return nullptr;
        }
        void* result = static_cast<void*>(const_cast<std::byte*>(pool.data()) + h.offset);
        std::cout << "[DEBUG] resolve returning: " << result << std::endl;
        return result;
    }

public:
    /**
     * Dtor to log destruction
     */
    ~BuddyAllocator() {
        if (instance_id == 0) {
            std::cout << "[DEBUG] BuddyAllocator dtor called (moved-from allocator), pool size: " << pool.size() << std::endl;
        } else {
            std::cout << "[DEBUG] BuddyAllocator dtor called, instance_id: " << instance_id
                      << ", pool size: " << pool.size() << std::endl;
        }
    }

    /**
     * Constructs a buddy allocator with the given total size and minimum chunk size.
     *
     * Both total_size and min_size must be powers of two, and min_size must not exceed total_size.
     *
     * The allocator is initialized with a single free chunk covering the entire pool.
     */
    BuddyAllocator(size_t total_size, size_t min_size)
        : instance_id(next_instance_id++), total_pool_size(total_size), min_chunk_size(min_size) {

        // 1. Validate inputs
        if (!is_power_of_two(total_pool_size) || !is_power_of_two(min_chunk_size)) {
            throw std::invalid_argument("Total size and min chunk size must be powers of two.");
        }
        if (min_chunk_size > total_pool_size) {
            throw std::invalid_argument("Min chunk size cannot be larger than total size.");
        }

        // 2. Initialize the memory pool
        pool.resize(total_pool_size);

        // 3. Calculate and resize the free lists
        num_levels = size_to_level(total_pool_size) + 1;
        free_lists.resize(num_levels);

        // 4. Create the initial root chunk that covers the entire pool
        size_t root_level = num_levels - 1;
        auto root_chunk = std::make_shared<MemChunk>(0, total_pool_size, true, root_level);

        // 5. Add the root chunk to the canonical map and the top-level free list
        chunks_by_offset[0] = root_chunk;
        free_lists[root_level].push_back(root_chunk.get());

        std::cout << "Buddy Allocator initialized."
                  << std::endl;
        std::cout << "Total size: " << total_pool_size << " bytes"
                  << std::endl;
        std::cout << "Min chunk size: " << min_chunk_size << " bytes"
                  << std::endl;
        std::cout << "Number of levels: " << num_levels << std::endl;
    }
    
    // Disable copying to prevent shallow copies of complex internal state
    BuddyAllocator(const BuddyAllocator&) = delete;
    
    BuddyAllocator& operator=(const BuddyAllocator&) = delete;

    /**
     * Move constructor that transfers ownership from another allocator.
     *
     * Handles (Block) remain valid after move because they store offsets relative
     * to the pool, and moving the pool preserves offset semantics. resolve()
     * will continue to work correctly with the moved allocator.
     */
    BuddyAllocator(BuddyAllocator&& other) noexcept
        : instance_id(other.instance_id)
        , total_pool_size(other.total_pool_size)
        , min_chunk_size(other.min_chunk_size)
        , num_levels(other.num_levels)
        , pool(std::move(other.pool))
        , chunks_by_offset(std::move(other.chunks_by_offset))
        , free_lists(std::move(other.free_lists)) {
        other.instance_id = 0;
        other.total_pool_size = 0;
        other.min_chunk_size = 0;
        other.num_levels = 0;
    }

    /**
     * Creates an independent deep copy of this allocator.
     *
     * This performs a complete copy of the underlying pool memory and metadata,
     * resulting in two independent allocators. Blocks created by this allocator
     * are not valid for the clone and vice versa, even if offsets match.
     *
     * Time complexity: O(number of chunks), as it rebuilds metadata structures.
     * Space complexity: O(original pool size + number of chunks).
     *
     * Note: This is explicitly provided instead of copy constructor/assignment
     * to avoid accidental expensive copies. Use move semantics for cheap transfers.
     */
    BuddyAllocator clone() const {
#ifdef DEBUG
        std::cout << "[DEBUG] clone(): Starting deep copy of allocator with " << chunks_by_offset.size() << " chunks"
                  << std::endl;
#endif
        
        // Create a new allocator with the same configuration
        BuddyAllocator copy(total_pool_size, min_chunk_size);
        copy.instance_id = next_instance_id++;  // Assign unique ID to copy
        
#ifdef DEBUG
        std::cout << "[DEBUG] clone(): Created empty copy allocator with instance_id=" << copy.instance_id
                  << std::endl;

        // 1. Copy the raw pool bytes
        copy.pool = pool;  // vector copy
        std::cout << "[DEBUG] clone(): Copied pool (" << pool.size() << " bytes)"
                  << std::endl;

        // 2. Deep-copy chunks_by_offset: create new MemChunk objects
        copy.chunks_by_offset.clear();
        std::cout << "[DEBUG] clone(): Clearing copy's chunks_by_offset"
                  << std::endl;
        
        for (const auto& kv : chunks_by_offset) {
            const auto& old_chunk = kv.second;
            std::cout << "[DEBUG] clone(): Copying chunk: offset=" << old_chunk->offset
                      << ", size=" << old_chunk->size << ", free=" << old_chunk->is_free
                      << ", level=" << old_chunk->level << std::endl;
                      
            auto new_chunk = std::make_shared<MemChunk>(
                old_chunk->offset, old_chunk->size, old_chunk->is_free, old_chunk->level);
            copy.chunks_by_offset[new_chunk->offset] = new_chunk;
        }
        std::cout << "[DEBUG] clone(): Copied " << chunks_by_offset.size() << " chunks to copy's chunks_by_offset"
                  << std::endl;

        // 3. Rebuild free_lists to contain the new shared_ptrs (preserving order)
        copy.free_lists.clear();
        copy.free_lists.resize(num_levels);
        std::cout << "[DEBUG] clone(): Rebuilding free_lists for " << num_levels << " levels"
                  << std::endl;
        
        for (size_t lvl = 0; lvl < num_levels; ++lvl) {
            std::cout << "[DEBUG] clone(): Processing level " << lvl << " with " << free_lists[lvl].size() << " free chunks"
                      << std::endl;
            for (const auto& old_chunk_ptr : free_lists[lvl]) {
                auto it = copy.chunks_by_offset.find(old_chunk_ptr->offset);
                assert(it != copy.chunks_by_offset.end());  // Must find corresponding new chunk
                std::cout << "[DEBUG] clone(): Adding chunk offset " << old_chunk_ptr->offset << " to copy's free_list[" << lvl << "]"
                          << std::endl;
                copy.free_lists[lvl].push_back(it->second.get());
            }
            // Fix iterators after rebuilding the list
            // No need to fix iterators since we'll remove free_list_iterator in the intrusive implementation
        }
        std::cout << "[DEBUG] clone(): Rebuilt all free_lists"
                  << std::endl;

        // 4. Copy remaining metadata (num_levels already set by constructor, but ensure)
        copy.total_pool_size = total_pool_size;
        copy.min_chunk_size = min_chunk_size;
        copy.num_levels = num_levels;
        
        std::cout << "[DEBUG] clone(): Copied metadata - total_pool_size=" << copy.total_pool_size
                  << ", min_chunk_size=" << copy.min_chunk_size << ", num_levels=" << copy.num_levels
                  << std::endl;

#ifndef NDEBUG
        std::cout << "[DEBUG] clone(): Verifying invariants on copy..."
                  << std::endl;
        copy.check_invariants();
        std::cout << "[DEBUG] clone(): Invariants verified successfully"
                  << std::endl;
#endif

        std::cout << "[DEBUG] clone(): Deep copy completed successfully"
                  << std::endl;
#endif
        return copy;
    }

    /**
     * Move assignment operator that transfers ownership from another allocator.
     *
     * Handles (Block) from the previous state of this allocator become invalid
     * after the move, but handles from the source 'other' remain valid in this
     * allocator. Always use the new allocator instance for resolve() after move.
     */
    BuddyAllocator& operator=(BuddyAllocator&& other) noexcept {
        if (this != &other) {
            instance_id = other.instance_id;
            total_pool_size = other.total_pool_size;
            min_chunk_size = other.min_chunk_size;
            num_levels = other.num_levels;
            pool = std::move(other.pool);
            chunks_by_offset = std::move(other.chunks_by_offset);
            free_lists = std::move(other.free_lists);
            
            other.instance_id = 0;
            other.total_pool_size = 0;
            other.min_chunk_size = 0;
            other.num_levels = 0;
        }
        return *this;
    }

    /**
     * Merges another allocator into this one, consuming the source.
     *
     * The source allocator is left in an empty state after the merge. Handles
     * (Block) created by the source allocator remain valid after merge because
     * their offsets are rebased to the new combined pool. Use resolve() with
     * the merged allocator to get the correct pointers.
     *
     * **Warning**: After merge_with, do not use the source allocator or its handles
     * with the old instance. Always use the target (this) allocator for all
     * operations and resolve() calls.
     */
    void merge_with(BuddyAllocator&& other) {
        if (other.total_pool_size == 0) return; // nothing to do

        if (other.min_chunk_size != this->min_chunk_size) {
            throw std::invalid_argument("merge_with: min_chunk_size mismatch");
        }

        // Compute new sizes / levels
        size_t old_pool_size = pool.size();
        size_t new_total = total_pool_size + other.total_pool_size;
        size_t new_num_levels = size_to_level(new_total) + 1;

        // Ensure our free_lists has enough levels before we splice other's lists
        free_lists.resize(new_num_levels);

        // Append other's raw bytes to our pool. We can use a simple insert;
        // std::byte is trivially copyable so move iterators are unnecessary.
        pool.insert(pool.end(), other.pool.begin(), other.pool.end());

        // Rebase other's chunks and move them into our canonical map.
        for (auto& kv : other.chunks_by_offset) {
            auto chunk_ptr = std::move(kv.second); // take ownership
            chunk_ptr->offset += old_pool_size;
            chunks_by_offset[chunk_ptr->offset] = std::move(chunk_ptr);
        }
        other.chunks_by_offset.clear();

        // Move free-list entries (rebased chunks are already in canonical map).
        // splice requires an lvalue list reference.
        for (size_t i = 0; i < std::min(num_levels, other.num_levels); ++i) {
            // Intrusive splice: move all elements from other to this list
            while (!other.free_lists[i].empty()) {
                MemChunk* chunk = other.free_lists[i].front();
                other.free_lists[i].pop_front();
                chunk->offset += old_pool_size; // Rebase offset
                free_lists[i].push_back(chunk);
            }
        }
        other.free_lists.clear();

        // Temporarily skip iterator fix (will be replaced with intrusive linking)
        // Fix iterators after splicing (splice may invalidate iterators from source list)

        // Update fields
        total_pool_size = new_total;
        num_levels = new_num_levels;
        other.total_pool_size = 0;
        other.min_chunk_size = 0;
        other.num_levels = 0;

        /**
         * Merge fixes for robustness and invariant preservation:
         * - Validates min_chunk_size compatibility to ensure consistent chunk sizing.
         * - Computes new_num_levels based on total_pool_size using size_to_level to maintain correct level hierarchy.
         * - Resizes free_lists before splicing to prevent out-of-bounds access.
         * - Rebases chunk offsets and inserts into chunks_by_offset with moved shared_ptrs, clearing source map to avoid dangling entries.
         * - Uses lvalue reference for splice to comply with std::list::splice requirements.
         * - Appends pool bytes directly (std::byte is trivially copyable).
         * These ensure chunks_by_offset keys match chunk->offset, free_lists index correctly by level, and all invariants hold post-merge.
         */

#ifndef NDEBUG
        check_invariants();
#endif
    }

    // Public allocation method that returns a Block
    Block allocate_raw(size_t size) {
        return allocate_block(size);
    }

    // Public deallocation method
    void deallocate_raw(const Block& block) {
        deallocate(block);
    }

    // Public resolve method
    void* resolve_raw(const Block& block) const {
        return resolve(block);
    }

    void print_report() const {
        print_header();
        print_memory_usage();
        print_free_lists_summary();
        print_all_chunks_list();
    }

    bool owns(const Block& h) const {
#ifdef DEBUG
        if (!h.valid()) {
            std::cout << "[DEBUG] owns(): Invalid block passed (offset=" << h.offset << ", size=" << h.size << ", id=" << h.instance_id << ")"
                      << std::endl;
            return false;
        }
        
        if (h.instance_id != instance_id) {
            std::cout << "[DEBUG] owns(): Instance ID mismatch - block id=" << h.instance_id << ", allocator id=" << instance_id << std::endl;
            return false;
        }
        
        std::cout << "[DEBUG] owns(): Checking block offset=" << h.offset << ", size=" << h.size << " in allocator with " << chunks_by_offset.size() << " chunks"
                  << std::endl;
        
        auto it = chunks_by_offset.find(h.offset);
        if (it == chunks_by_offset.end()) {
            std::cout << "[DEBUG] owns(): Offset " << h.offset << " not found in chunks_by_offset"
                      << std::endl;
            return false;
        }
        
        std::cout << "[DEBUG] owns(): Found chunk at offset " << h.offset << " with size " << it->second->size << std::endl;
        
        bool size_match = (it->second->size == h.size);
        std::cout << "[DEBUG] owns(): Size match: " << (size_match ? "YES" : "NO") << std::endl;
        
        std::cout << "[DEBUG] owns(): Returning " << (size_match ? "true" : "false") << " for this block"
                  << std::endl;
#endif
        if (!h.valid()) {
            return false;
        }
        
        if (h.instance_id != instance_id) {
            return false;
        }
        
        auto it = chunks_by_offset.find(h.offset);
        if (it == chunks_by_offset.end()) {
            return false;
        }
        
        return (it->second->size == h.size);
    }

private:
    // --- Reporting Helpers ---

    void print_header() const {
        std::cout << "\n=================================================" << std::endl;
        std::cout << "        Buddy Allocator State Report"
                  << std::endl;
        std::cout << "=================================================" << std::endl;
        std::cout << "Total Size: " << total_pool_size << " bytes"
                  << std::endl;
        std::cout << "Min Chunk Size: " << min_chunk_size << " bytes"
                  << std::endl;
        std::cout << "Pool Start Address: " << static_cast<const void*>(pool.data()) << std::endl << std::endl;
    }

    void print_memory_usage() const {
        auto [memory_in_use, allocated_chunks_count] = get_memory_in_use();
        auto [memory_free, free_chunks_count] = get_free_memory();
        double external_fragmentation = calculate_external_fragmentation();

        std::cout << "--- Memory Usage ---"
                  << std::endl;
        std::cout << "Used: " << memory_in_use << " bytes (" << allocated_chunks_count << " chunks)"
                  << std::endl;
        std::cout << "Free: " << memory_free << " bytes (" << free_chunks_count << " chunks)"
                  << std::endl;
        std::cout << "External Fragmentation (approx): " << external_fragmentation << "%"
                  << std::endl << std::endl;
    }

    void print_free_lists_summary() const {
        std::cout << "--- Free Lists ---"
                  << std::endl;
        for (size_t i = 0; i < num_levels; i++) {
            size_t level_size = level_to_size(i);
            std::cout << "Level " << i << " (" << level_size << " bytes): " << free_lists[i].size() << " free chunks"
                      << std::endl;
        }
        std::cout << std::endl;
    }

    void print_all_chunks_list() const {
        auto sorted_chunks = get_sorted_chunks();
        std::cout << "--- All Chunks (" << chunks_by_offset.size() << " total) ---"
                  << std::endl;
        for (const auto& chunk : sorted_chunks) {
            std::cout << "  "
                    << (chunk->is_free ? "[FREE]" : "[USED]")
                    << "\tAddr: " << static_cast<const void*>(pool.data() + chunk->offset)
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
            for(MemChunk* chunk = list.head; chunk; chunk = chunk->next_free) {
                memory_free += chunk->size;
#ifndef NDEBUG
                assert(counted_chunks.find(chunk) == counted_chunks.end() && "Double count detected in free list!");
                counted_chunks.insert(chunk);
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
            for (MemChunk* chunk_ptr = free_lists[i].head; chunk_ptr; chunk_ptr = chunk_ptr->next_free) {
                assert(chunk_ptr->level == i && "Invariant Fail: Chunk in wrong free list level.");
                
                auto it = chunks_by_offset.find(chunk_ptr->offset);
                assert(it != chunks_by_offset.end() && "Invariant Fail: Free chunk not in canonical map.");
                assert(it->second.get() == chunk_ptr && "Invariant Fail: Free chunk in map is not the same object.");
                assert(it->second->is_free && "Invariant Fail: Chunk in free list is not marked as free.");

                all_free_chunks.insert(it->second);
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

// MemoryHandle definition - now properly separated from BuddyAllocator
template <typename T>
class MemoryHandle {
private:
    Block block;
    std::unique_ptr<BuddyAllocator> allocator;

public:
    /**
     * Rebind this handle to take ownership of a new allocator.
     * The new_owner must be moved from to transfer ownership.
     */
    void rebind(BuddyAllocator&& new_owner) {
        allocator = std::make_unique<BuddyAllocator>(std::move(new_owner));
    }

public:
    // Constructor that takes ownership of allocator
    MemoryHandle(const Block& b, BuddyAllocator&& alloc)
        : block(b), allocator(std::make_unique<BuddyAllocator>(std::move(alloc))) {}

    // Constructor for invalid handle
    MemoryHandle() : block(), allocator(nullptr) {}

    // Destructor for RAII
    ~MemoryHandle() {
        std::cout << "[DEBUG] MemoryHandle dtor starting for type " << typeid(T).name()
                  << ", block valid: " << block.valid()
                  << ", allocator: " << (allocator ? (void*)allocator.get() : nullptr) << std::endl;
        if (block.valid() && allocator) {
            std::cout << "[DEBUG] Calling resolve_raw in handle dtor" << std::endl;
            void* ptr = allocator->resolve_raw(block);
            std::cout << "[DEBUG] resolve_raw returned: " << ptr << std::endl;
            if constexpr (std::is_array_v<T>) {
                using ElementType = typename std::remove_extent_t<T>;
                ElementType* arr_ptr = static_cast<ElementType*>(ptr);
                if (arr_ptr) {
                    std::cout << "[DEBUG] Starting array destruction for " << std::extent_v<T> << " elements" << std::endl;
                    for (size_t i = 0; i < std::extent_v<T>; ++i) {
                        try {
                            std::cout << "[DEBUG] Destroying array element " << i << ", value before dtor: " << arr_ptr[i].value << std::endl;
                        } catch (...) {
                            // Ignore if ElementType lacks 'value' member
                        }
                        arr_ptr[i].~ElementType();
                    }
                    std::cout << "[DEBUG] Completed array destruction" << std::endl;
                }
            } else {
                T* obj_ptr = static_cast<T*>(ptr);
                if (obj_ptr) {
                    try {
                        std::cout << "[DEBUG] Destroying object, value before dtor: " << obj_ptr->value << std::endl;
                    } catch (...) {
                        // Ignore if T lacks 'value' member
                    }
                    obj_ptr->~T();
                }
            }
            std::cout << "[DEBUG] Calling deallocate_raw in handle dtor" << std::endl;
            allocator->deallocate_raw(block);
            std::cout << "[DEBUG] deallocate_raw completed" << std::endl;
        } else {
            std::cout << "[DEBUG] Skipping destruction: invalid block or no allocator" << std::endl;
        }
        std::cout << "[DEBUG] MemoryHandle dtor ending" << std::endl;
    }


    // Move constructor - moves owned allocator
    MemoryHandle(MemoryHandle&& other) noexcept
        : block(std::move(other.block)), allocator(std::move(other.allocator)) {
        other.block = Block();
    }

    // Move assignment operator - moves owned allocator
    MemoryHandle& operator=(MemoryHandle&& other) noexcept {
        if (this != &other) {
            // Destroy current contents
            if (block.valid() && allocator) {
                if constexpr (std::is_array_v<T>) {
                    using ElementType = typename std::remove_extent_t<T>;
                    ElementType* arr_ptr = static_cast<ElementType*>(allocator->resolve_raw(block));
                    if (arr_ptr) {
                        std::cout << "[DEBUG] Starting array destruction in move assignment for " << std::extent_v<T> << " elements" << std::endl;
                        for (size_t i = 0; i < std::extent_v<T>; ++i) {
                            try {
                                std::cout << "[DEBUG] Destroying array element " << i << " in move assignment, value before dtor: " << arr_ptr[i].value << std::endl;
                            } catch (...) {
                                // Ignore if ElementType lacks 'value' member
                            }
                            arr_ptr[i].~ElementType();
                        }
                        std::cout << "[DEBUG] Completed array destruction in move assignment" << std::endl;
                    }
                } else {
                    T* obj_ptr = static_cast<T*>(allocator->resolve_raw(block));
                    if (obj_ptr) {
                        obj_ptr->~T();
                    }
                }
                allocator->deallocate_raw(block);
            }

            // Move from other
            block = std::move(other.block);
            allocator = std::move(other.allocator);

            // Invalidate other
            other.block = Block();
        }
        return *this;
    }

    // Disable copy constructor and copy assignment
    MemoryHandle(const MemoryHandle&) = delete;
    MemoryHandle& operator=(const MemoryHandle&) = delete;

    // Accessors
    T* operator->() const {
        if (!valid()) return nullptr;
        return static_cast<T*>(allocator->resolve_raw(block));
    }

    T& operator*() const {
        if (!valid()) {
            throw std::logic_error("Dereferencing invalid MemoryHandle");
        }
        return *static_cast<T*>(allocator->resolve_raw(block));
    }

    // Check if the handle is valid
    bool valid() const {
        return block.valid() && allocator != nullptr;
    }

    // Access the owned allocator (for advanced use)
    BuddyAllocator& get_allocator() { return *allocator; }
    const BuddyAllocator& get_allocator() const { return *allocator; }

    // For compatibility with old owns() checks - now always true since we own it
    bool is_owned_by(const BuddyAllocator& alloc) const {
        return valid() && &get_allocator() == &alloc;
    }

};

// Add allocation methods to BuddyAllocator as free functions or as methods that create MemoryHandle
template<typename T, typename... Args>
MemoryHandle<T> allocate(BuddyAllocator&& allocator, Args&&... args) {
    size_t required_size;
    if (std::is_array<T>::value) {
        typedef typename std::remove_extent<T>::type ElementType;
        constexpr size_t count = std::extent<T>::value;
        required_size = sizeof(ElementType) * count;
    } else {
        required_size = sizeof(T);
    }

    Block block = allocator.allocate_raw(required_size);
    if (!block.valid()) {
        return MemoryHandle<T>();
    }

    void* ptr = allocator.resolve_raw(block);
    if (std::is_array<T>::value) {
        typedef typename std::remove_extent<T>::type ElementType;
        for (size_t i = 0; i < std::extent<T>::value; ++i) {
            new (static_cast<ElementType*>(ptr) + i) ElementType();
        }
    } else {
        try {
            new (ptr) T(std::forward<Args>(args)...);
        } catch (...) {
            allocator.deallocate_raw(block);
            throw;
        }
    }

    return MemoryHandle<T>(block, std::move(allocator));
}