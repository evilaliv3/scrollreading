
#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <blosc2.h>
#include <omp.h>
long zarrMissingChunks_1 = 0; // chunks read as zeros because the file was missing


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

/* ---------------------------------------------------------------------------------------------
   A store of decompressed chunks, shared by every reader in the process.

   Why it exists. Each reader keeps its own buffers and cannot see what another reader has already
   decompressed, so a g 300 run loaded 1651 chunks that were only 316 distinct ones: a reload
   factor of 5.2, measured on 2026-09-14. Of the 32.5 ms a load costs, 99.7 % is
   blosc2_decompress and 0.2 % is fread, so the repeated work is decompression and nothing else.
   Serving a repeat from here is a 7 MB memcpy instead.

   Why it is safe. Readers keep their own buffers: the store hands out a copy, never a pointer, so
   nothing aliases and a writer cannot disturb another reader. And the moment anything is written
   to a zarr, that zarr is marked and the store stops serving and storing it, so a written chunk
   is always reread from the file it was flushed to. Growth never writes the surface field, which
   is where the whole cost is.

   SIMPAPER_SHARED_CHUNKS=0 turns it off and restores the previous behaviour exactly.
   --------------------------------------------------------------------------------------------- */
#define ZARR_STORE_SLOTS 8192
typedef struct { unsigned long key; int c[3]; ZARRType_1 *data; } ZarrStoreEntry;
static ZarrStoreEntry zarrStore[ZARR_STORE_SLOTS];
static long zarrStoreCount = 0;
static long zarrStoreCap = -1;          /* how many chunks may be held, -1 = not decided yet */
static unsigned long zarrDirty[64];
static int zarrDirtyCount = 0;
long zarrStoreHits_1 = 0, zarrStoreMisses_1 = 0, zarrStoreHeld_1 = 0;

static unsigned long zarrRootKeyOf(const char *s)
{
	unsigned long h = 1469598103934665603UL;
	while (*s) { h ^= (unsigned char)*s++; h *= 1099511628211UL; }
	return h ? h : 1;
}

static int zarrIsDirty(unsigned long rootKey)
{
	int i;
	for(i = 0; i < zarrDirtyCount; i++) if (zarrDirty[i] == rootKey) return 1;
	return 0;
}

void ZARRMarkWritten_1(ZARR_1 *z)
{
	#pragma omp critical(zarrStore)
	{
		if (!zarrIsDirty(z->rootKey) && zarrDirtyCount < 64)
			zarrDirty[zarrDirtyCount++] = z->rootKey;
	}
}

static long zarrStoreSlot(unsigned long rootKey, int c[3])
{
	unsigned long h = rootKey ^ ((unsigned long)c[0]*2654435761UL)
	                ^ ((unsigned long)c[1]*40503UL) ^ ((unsigned long)c[2]*2246822519UL);
	long i, slot = (long)(h % ZARR_STORE_SLOTS);
	for(i = 0; i < ZARR_STORE_SLOTS; i++)
	{
		long k = (slot+i) % ZARR_STORE_SLOTS;
		if (!zarrStore[k].data) return k;                       /* free */
		if (zarrStore[k].key == rootKey && zarrStore[k].c[0]==c[0]
		    && zarrStore[k].c[1]==c[1] && zarrStore[k].c[2]==c[2]) return k;   /* hit */
	}
	return -1;
}

/* How many chunks to hold. A quarter of the memory the machine reports as available, which on a
   small machine is a few and on this one is thousands: the run needs 316. */
static void zarrStoreDecideCap(void)
{
	const char *e = getenv("SIMPAPER_SHARED_CHUNKS");
	if (e && atoi(e) == 0) { zarrStoreCap = 0; return; }
	{
		long available_kb = 0;
		FILE *mi = fopen("/proc/meminfo","r");
		if (mi)
		{
			char line[256];
			while(fgets(line,sizeof(line),mi))
				if (sscanf(line,"MemAvailable: %ld kB",&available_kb)==1) break;
			fclose(mi);
		}
		zarrStoreCap = available_kb > 0 ? (long)((available_kb/1024.0) * 0.25 / 6.75) : 256;
		if (zarrStoreCap > ZARR_STORE_SLOTS/2) zarrStoreCap = ZARR_STORE_SLOTS/2;
		if (zarrStoreCap < 0) zarrStoreCap = 0;
		if (e && atoi(e) > 0 && atoi(e) < zarrStoreCap) zarrStoreCap = atoi(e);
	}
}

/* Returns 1 and fills dest when the chunk is already decompressed somewhere in the process. */
static int zarrStoreTake(ZARR_1 *z, int c[3], ZARRType_1 *dest)
{
	int got = 0;
	#pragma omp critical(zarrStore)
	{
		if (zarrStoreCap < 0) zarrStoreDecideCap();
		if (zarrStoreCap > 0 && !zarrIsDirty(z->rootKey))
		{
			long k = zarrStoreSlot(z->rootKey, c);
			if (k >= 0 && zarrStore[k].data)
			{
				memcpy(dest, zarrStore[k].data, sizeof(ZARRType_1)*7077888);
				zarrStoreHits_1++; got = 1;
			}
		}
	}
	return got;
}

static void zarrStoreKeep(ZARR_1 *z, int c[3], const ZARRType_1 *src)
{
	#pragma omp critical(zarrStore)
	{
		if (zarrStoreCap > 0 && zarrStoreCount < zarrStoreCap && !zarrIsDirty(z->rootKey))
		{
			long k = zarrStoreSlot(z->rootKey, c);
			if (k >= 0 && !zarrStore[k].data)
			{
				ZARRType_1 *p = (ZARRType_1 *)malloc(sizeof(ZARRType_1)*7077888);
				if (p)
				{
					memcpy(p, src, sizeof(ZARRType_1)*7077888);
					zarrStore[k].key = z->rootKey;
					zarrStore[k].c[0]=c[0]; zarrStore[k].c[1]=c[1]; zarrStore[k].c[2]=c[2];
					zarrStore[k].data = p;
					zarrStoreCount++; zarrStoreHeld_1 = zarrStoreCount;
				}
			}
		}
		zarrStoreMisses_1++;
	}
}

ZARR_1 *ZARROpen_1(const char *location)
{
	ZARR_1 *z = (ZARR_1 *)malloc(sizeof(ZARR_1));
	
	z->locationRootLength = strlen(location);
	z->rootKey = zarrRootKeyOf(location);
	
	z->location = (char *)malloc(z->locationRootLength + 1 + 100);
	
	z->buffer = NULL;
	z->index = -1;

	// How many buffers to keep. Default: all 80, as before. Code that opens several readers at
	// once can lower it through ZARR_BUFFERS to keep resident memory down: each buffer holds one
	// decompressed 192^3 chunk, that is 6.75 MiB.
	{
		const char *e = getenv("ZARR_BUFFERS");
		int n = e ? atoi(e) : 80;
		if (n < 1) n = 1;
		if (n > 80) n = 80;
		z->nbuf = n;
	}

    for(int i = 0; i<z->nbuf; i++)
	{
	  z->written[i] = 0;
      for(int j = 0; j<3; j++)
        z->bufferIndex[i][j] = -1;
	}
	
    strcpy(z->location,location);
	
	z->counter = 1;

    for(int i = 0; i<z->nbuf; i++)
      z->bufferUsed[i] = 0;
  
	return z;
}

int ZARRFlushOne_1(ZARR_1 *z, int i)
{
    if (z->written[i])
	{
      sprintf(z->location+z->locationRootLength,"/%d/%d/%d",z->bufferIndex[i][0],z->bufferIndex[i][1],z->bufferIndex[i][2]);

      blosc1_set_compressor("zstd");
	  int compressed_len = blosc2_compress(1,1,sizeof(ZARRType_1),z->buffers[i],sizeof(ZARRType_1)*7077888,z->compressedData,sizeof(ZARRType_1)*7077888+BLOSC2_MAX_OVERHEAD);

      if (compressed_len <= 0) {
        return -1;
      }

	  FILE *f = fopen(z->location,"wb");
	  fwrite(z->compressedData,1,compressed_len,f);
	  fclose(f);
	
	  z->written[i] = 0;
	}
	
	return 0;
}

int ZARRFlush_1(ZARR_1 *z)
{
	for(int i = 0; i<z->nbuf; i++)
	{
		if (z->bufferIndex[i][0] != -1)
		{
			ZARRFlushOne_1(z,i);
		}
	}
	
	return 0;
}

int ZARRClose_1(ZARR_1 *z)
{
    ZARRFlush_1(z);
    free(z->location);
    free(z);
	
	return 0;
}

/* Make buffer point to the chunk and index contain the index*/
int ZARRCheckChunk_1(ZARR_1 *z, int c[3])
{	
	if (z->buffer && c[0] == z->bufferIndex[z->index][0] && c[1] == z->bufferIndex[z->index][1] && c[2] == z->bufferIndex[z->index][2])
		return 0;

	for(z->index = 0; z->index < z->nbuf; z->index++)
	{
		if (1  && c[0] == z->bufferIndex[z->index][0] && c[1] == z->bufferIndex[z->index][1] && c[2] == z->bufferIndex[z->index][2])
		{
			z->buffer = &z->buffers[z->index];
			z->bufferUsed[z->index] = z->counter++;
			return 0;
		}
	}
		
	for(z->index = 0; z->index < z->nbuf; z->index++)
	{
		if (z->bufferIndex[z->index][0]==-1)
			break;
	}
  	
	if (z->index == z->nbuf)
	{
		/* Find the buffer that was least recently used and free it up */
		printf("Ran out of buffers - flushing oldest\n");
		int oldestIndex = 0;
		uint64_t oldestAge = z->bufferUsed[0];
		for(int i = 1; i<z->nbuf; i++)
		{
			if (z->bufferUsed[i]<oldestAge)
			{
				oldestAge = z->bufferUsed[i];
				oldestIndex = i;
			}
		}
		
		ZARRFlushOne_1(z,oldestIndex);
		z->index = oldestIndex;
	}
    z->bufferIndex[z->index][0] = c[0];
    z->bufferIndex[z->index][1] = c[1];
    z->bufferIndex[z->index][2] = c[2];

	z->buffer = &z->buffers[z->index];
	z->bufferUsed[z->index] = z->counter++;
	z->written[z->index] = 0;
    sprintf(z->location+z->locationRootLength,"/%d/%d/%d",z->bufferIndex[z->index][0],z->bufferIndex[z->index][1],z->bufferIndex[z->index][2]);

    if (zarrStoreTake(z,c,(ZARRType_1 *)z->buffer)) return 0;

    FILE *f = fopen(z->location,"rb");

    //printf("Opening:%s\n",z->location);	
	if (!f)
	{
		printf("Did not find file:%s\n",z->location); // Useful to display this message because it often indicates a file naming problem
		zarrMissingChunks_1++;
		// The path goes to a list the caller can check upstream: in zarr an absent chunk means
		// fill_value by definition, so a missing file is only an error when the chunk does exist
		// in the source the box was copied from.
		{
			FILE *ml = fopen(getenv("ZARR_MISSING_LIST") ? getenv("ZARR_MISSING_LIST") : "/dev/null","a");
			if (ml) { fprintf(ml,"%s\n",z->location); fclose(ml); }
		}

		memset(z->buffer,0,sizeof(ZARRType_1)*7077888);

		//No need to count it as written to yet - if it remains empty then just leave the file as non-existent
	    //z->written[z->index] = 1;
		
		return 0;
	}

	if (f)
	{
		fseek(f,0,SEEK_END);
		long fsize = ftell(f);
		fseek(f,0,SEEK_SET);
        fread(z->compressedData,1,fsize,f);
		fclose(f);
		

        blosc1_set_compressor("zstd");
        int decompressed_size = blosc2_decompress(z->compressedData, fsize, z->buffer, sizeof(ZARRType_1)*7077888);
        if (decompressed_size < 0) {
            return 0;
        }
        zarrStoreKeep(z,c,(const ZARRType_1 *)z->buffer);

	}	
	
	return 0;
}
ZARRType_1 ZARRReadRO_1(const ZARR_1 *za,int x0,int x1,int x2,int *found,int *hint)
{
	int c0 = x0/192, m0 = x0%192;
	int c1 = x1/192, m1 = x1%192;
	int c2 = x2/192, m2 = x2%192;

	/* Fast path: the block found last time by THIS caller. The hint lives on the caller's
	   stack, so it is private per thread while the cache stays shared and immutable. Without
	   it every read scanned all the buffers: measured on 2026-09-11, with about 124 reads per
	   tap all landing in the same block, that scan costs 24 times a single comparison, and it
	   was the reason the shared version stayed slow. */
	if (hint)
	{
		int h = *hint;
		if (h >= 0 && h < za->nbuf &&
		    c0 == za->bufferIndex[h][0] && c1 == za->bufferIndex[h][1] && c2 == za->bufferIndex[h][2])
		{
			*found = 1;
			return za->buffers[h][m0][m1][m2];
		}
	}

	for(int i = 0; i<za->nbuf; i++)
		if (c0 == za->bufferIndex[i][0] && c1 == za->bufferIndex[i][1] && c2 == za->bufferIndex[i][2])
		{
			*found = 1;
			if (hint) *hint = i;
			return za->buffers[i][m0][m1][m2];
		}
	*found = 0;
	return 0;
}

/* Reads an n0 x n1 x n2 block into a local array in one go.
   Why: computing the field performs 1.8 billion reads per run, and every one of them goes
   through ZARRRead_1, which does three divisions, three moduli and the chunk check. But all the
   taps needed by one smoothing evaluation lie inside a 9x9x9 cube, that is 729 bytes.
   The buffer's last index is contiguous, so a row of n2 voxels is n2 adjacent bytes: one chunk
   lookup and n0*n1 row copies, instead of n0*n1*n2 reads.
   Returns 1 when the block lay entirely inside a single chunk and was filled, 0 otherwise; in
   that case the caller takes the old path, so no voxel is ever read that the original code
   would not read. */
int ZARRReadBlock_1(ZARR_1 *za,int x0,int x1,int x2,int n0,int n1,int n2,ZARRType_1 *outBlock)
{
	int c[3];
	c[0] = x0/192; c[1] = x1/192; c[2] = x2/192;
	/* does the whole block sit inside one chunk? */
	if ((x0+n0-1)/192 != c[0] || (x1+n1-1)/192 != c[1] || (x2+n2-1)/192 != c[2])
		return 0;
	ZARRCheckChunk_1(za,c);
	{
		int m0 = x0%192, m1 = x1%192, m2 = x2%192;
		int i,j;
		for(i = 0; i<n0; i++)
			for(j = 0; j<n1; j++)
				memcpy(outBlock + ((size_t)i*n1 + j)*n2, &(*za->buffer)[m0+i][m1+j][m2], n2);
	}
	return 1;
}
ZARRType_1 ZARRRead_1(ZARR_1 *za,int x0,int x1,int x2)
{
	int c[3],m[3];
    c[0] = x0/192;
    m[0] = x0%192;
    c[1] = x1/192;
    m[1] = x1%192;
    c[2] = x2/192;
    m[2] = x2%192;
	
	ZARRCheckChunk_1(za,c);
	
    return (*za->buffer)[m[0]][m[1]][m[2]];	  	  	  	  
}


// Read several values from the last dimensions
void ZARRReadN_1(ZARR_1 *za,int x0,int x1,int x2,int n,ZARRType_1 *v)
{
	int c[3],m[3];
    c[0] = x0/192;
    m[0] = x0%192;
    c[1] = x1/192;
    m[1] = x1%192;
    c[2] = x2/192;
    m[2] = x2%192;
	
	ZARRCheckChunk_1(za,c);
			  
	memcpy(v,&(*za->buffer)[m[0]][m[1]][m[2]],n*sizeof(ZARRType_1));
}

int ZARRWrite_1(ZARR_1 *za,int x0,int x1,int x2,ZARRType_1 value)
{
	int c[3],m[3];
    c[0] = x0/192;
    m[0] = x0%192;
    c[1] = x1/192;
    m[1] = x1%192;
    c[2] = x2/192;
    m[2] = x2%192;
	
	ZARRCheckChunk_1(za,c);
			  
	(*za->buffer)[m[0]][m[1]][m[2]] = value;

	za->written[za->index] = 1; ZARRMarkWritten_1(za);
	
	return 0;
}


void ZARRWriteN_1(ZARR_1 *za,int x0,int x1,int x2,int n, ZARRType_1 *v)
{
	int c[3],m[3];
    c[0] = x0/192;
    m[0] = x0%192;
    c[1] = x1/192;
    m[1] = x1%192;
    c[2] = x2/192;
    m[2] = x2%192;
	
	ZARRCheckChunk_1(za,c);
			  
	memcpy(&(*za->buffer)[m[0]][m[1]][m[2]],v,n*sizeof(ZARRType_1));

	za->written[za->index] = 1; ZARRMarkWritten_1(za);  
}

// Assumes that we have already written at least once to this chunk
void ZARRNoCheckWriteN_1(ZARR_1 *za,int x0,int x1,int x2,int n, ZARRType_1 *v)
{
	int m[3];
    m[0] = x0%192;
    m[1] = x1%192;
    m[2] = x2%192;
	
	memcpy(&(*za->buffer)[m[0]][m[1]][m[2]],v,n*sizeof(ZARRType_1));

	za->written[za->index] = 1; ZARRMarkWritten_1(za);  
}
