#pragma once
#include <blosc2.h>
typedef uint8_t ZARRType_1;

typedef struct {
    int locationRootLength;
    char *location;
  
    unsigned char compressedData[sizeof(ZARRType_1)*7077888+BLOSC2_MAX_OVERHEAD];
    ZARRType_1 buffers[80][192][192][192];
    int bufferIndex[80][3];
    unsigned char written[80];
    uint64_t bufferUsed[80];
    ZARRType_1 (*buffer)[192][192][192];

    int index;
  
    uint64_t counter;

    // How many of the 80 buffers are really used. The array keeps its fixed size, but pages
    // only become resident when touched, so limiting this number limits the memory actually
    // occupied. Chosen at run time from ZARR_BUFFERS, which the caller derives from the memory
    // available and the number of readers it will keep open at once.
    int nbuf;

    /* Identifies which zarr this is, for the shared store of decompressed chunks. */
    unsigned long rootKey;
} ZARR_1;

ZARR_1 *ZARROpen_1(const char *location);
int ZARRFlushOne_1(ZARR_1 *z, int i);
int ZARRFlush_1(ZARR_1 *z);
int ZARRClose_1(ZARR_1 *z);
/* Make buffer point to the chunk and index contain the index*/
int ZARRCheckChunk_1(ZARR_1 *z, int c[3]);
ZARRType_1 ZARRRead_1(ZARR_1 *za,int x0,int x1,int x2);
/* Reads a block into a local array, if the whole block lies inside one chunk. Returns 1 when it
   filled the array, 0 when the block straddles chunks and the caller must do as before. */
int ZARRReadBlock_1(ZARR_1 *za,int x0,int x1,int x2,int n0,int n1,int n2,ZARRType_1 *outBlock);
/* READ-ONLY read: looks for the chunk among those already cached without touching any mutable
   state (not the buffer, not the index, not the counter), so several threads can call it on the
   same reader. When the chunk is absent it returns 0 and sets *found to 0: it is up to the caller
   to have loaded it beforehand in a serial pass. */
ZARRType_1 ZARRReadRO_1(const ZARR_1 *za,int x0,int x1,int x2,int *found,int *hint);
// Read several values from the last dimensions
void ZARRReadN_1(ZARR_1 *za,int x0,int x1,int x2,int n,ZARRType_1 *v);
int ZARRWrite_1(ZARR_1 *za,int x0,int x1,int x2,ZARRType_1 value);
void ZARRWriteN_1(ZARR_1 *za,int x0,int x1,int x2,int n, ZARRType_1 *v);
// Assumes that we have already written at least once to this chunk
void ZARRNoCheckWriteN_1(ZARR_1 *za,int x0,int x1,int x2,int n, ZARRType_1 *v);
