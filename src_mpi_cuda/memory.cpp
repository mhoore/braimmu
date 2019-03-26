#include <mpi.h>
#include "math.h"

#include "memory.h"

using namespace brain_NS;

/* ----------------------------------------------------------------------*/
Memory::Memory() {
}

/* ----------------------------------------------------------------------*/
Memory::~Memory() {
}

/* ----------------------------------------------------------------------*/
void Memory::sfree(void *ptr) {
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------*/
void *Memory::smalloc(size_t nbytes, const char *name) {
  if (nbytes == 0) return NULL;

  void *ptr = malloc(nbytes);
  if (ptr == NULL) {
    printf("Failed to allocate %lld bytes for array %s. \n", nbytes,name);
    exit(1);
  }
  return ptr;
}
