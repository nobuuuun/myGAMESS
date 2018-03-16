/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Private declarations for handling distributed memory.
 * Added May 2014 by A.G.
 */

/** @brief Structure describing "memory holes" (unallocated segments below the top of memory stack) */
typedef struct {
  size_t offset; /* offset in the memory pool */
  size_t size;   /* size of the hole */
} memhole_t;

/**
   The list of memholes is kept at the top of the memory pool, at
   offsets above *gv(mem_attic), and grows down, pushing
   *gv(mem_attic) down. The memholes are sorted by their offsets in
   the memory pool in ascending order; two adjacent memholes always
   coalesce, and zero-sized memholes are always eliminated.

   Two memholes never have the same offset.

   Global variables used:
   gv(mem_addr): address of the memory pool
   gv(mem_attic): pointer to the current offset of the attic in the memory pool
   gv(mem_nmholes): pointer to the current number of memholes
   gv(mem_top): pointer to the current offset of the top of the memory pool

   Public API (declared in ddi_base.h):
   void DDI_Memory_Alloc(size_t sz, void** addr_ptr, size_t* offset_ptr);
   void DDI_Memory_Free(size_t offset, size_t size);
**/

/* Memory global variables: */
static char   *gv(mem_addr)  = NULL; /* address of the memory pool */
static size_t *gv(mem_total) = NULL;
static size_t *gv(mem_top)  = NULL; /* current offset of the top of the memory pool */
static size_t *gv(mem_max)   = NULL;
/* offset of the service information at the top of the memory pool: */
static size_t *gv(mem_attic) = NULL;
/* current number of memory holes (free memory segments below the top mark): */
static size_t *gv(mem_nmholes) = NULL;

# if defined DDI_LAPI
static size_t gv(mem_heap_max) = 0;
static size_t gv(mem_heap_total) = 0;
# endif

  
/** @brief Service function.
    Insert a new memhole pointing to offset @param offs, of size @param sz;
    
    @returns 1 if successful, 0 if out of available memory for the attic.
    
    @author Alexander Gaenko
    @date May 15 2014
**/
static int DDI_Memory_mhinsert(size_t offs, size_t sz);

/** @brief Service function.
    Indicates that the memhole at index @param idx is filled with data of size @sz,
    and returns its offset at @param offset_ptr

    @author Alexander Gaenko
    @date May 15 2014
**/
static void DDI_Memory_mhfill(size_t idx, size_t sz, size_t* offset_ptr);

/** @brief Service function.
    Finds a memhole of size not less than and as close as possible to @param sz; if found,
    fills @param idx_ptr with the index of the memhole and @param hsize_ptr with the size.

    @return 1 (true) if search is successful, 0 (false) otherwise

    @author Alexander Gaenko
    @date May 15 2014
*/
static int DDI_Memory_mhfind(size_t sz, size_t* idx_ptr, size_t* hsize_ptr);

/** @brief Service function.
    Checks if the topmost memhole touches mem_top, and erases it if so;
    fills @param offset_ptr with the offset of the memhole, if any.
    @returns 1 (true) if there was a top memhole, 0 (false) if not.

    @author Alexander Gaenko
    @date May 15 2014
*/
static int DDI_Memory_mhchecktop(size_t* offset_ptr);


/** @brief Sevice function.
    Calculates the amount of memory to be reserved
    for service data structures at the bottom of the memory stack.
    @returns Number of bytes to be reserved.

    @note This function should be changed as new variables are added to the 
*/
static size_t DDI_Memory_reserved_bytes();
