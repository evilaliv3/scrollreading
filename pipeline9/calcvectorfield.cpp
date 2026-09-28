#include <omp.h>
#include "calcvectorfield.h"

float kernel[5][5][5] = {
  {
    {1.0546170681737936e-05,0.00010989362203613103,0.00024002973835478245,0.00010989362203613103,1.0546170681737936e-05},
    {0.00010989362203613103,0.0011451178374281624,0.0025011673089900184,0.0011451178374281624,0.00010989362203613103},
    {0.00024002973835478245,0.0025011673089900184,0.005463051664281514,0.0025011673089900184,0.00024002973835478245},
    {0.00010989362203613103,0.0011451178374281624,0.0025011673089900184,0.0011451178374281624,0.00010989362203613103},
    {1.0546170681737936e-05,0.00010989362203613103,0.00024002973835478245,0.00010989362203613103,1.0546170681737936e-05},
  },
  {
    {0.00010989362203613103,0.0011451178374281624,0.0025011673089900184,0.0011451178374281624,0.00010989362203613103},
    {0.0011451178374281624,0.011932401874651296,0.026062761849591155,0.011932401874651296,0.0011451178374281624},
    {0.0025011673089900184,0.026062761849591155,0.056926305563971345,0.026062761849591155,0.0025011673089900184},
    {0.0011451178374281624,0.011932401874651296,0.026062761849591155,0.011932401874651296,0.0011451178374281624},
    {0.00010989362203613103,0.0011451178374281624,0.0025011673089900184,0.0011451178374281624,0.00010989362203613103},
  },
  {
    {0.00024002973835478245,0.0025011673089900184,0.005463051664281514,0.0025011673089900184,0.00024002973835478245},
    {0.0025011673089900184,0.026062761849591155,0.056926305563971345,0.026062761849591155,0.0025011673089900184},
    {0.005463051664281514,0.056926305563971345,0.12433848276956383,0.056926305563971345,0.005463051664281514},
    {0.0025011673089900184,0.026062761849591155,0.056926305563971345,0.026062761849591155,0.0025011673089900184},
    {0.00024002973835478245,0.0025011673089900184,0.005463051664281514,0.0025011673089900184,0.00024002973835478245},
  },
  {
    {0.00010989362203613103,0.0011451178374281624,0.0025011673089900184,0.0011451178374281624,0.00010989362203613103},
    {0.0011451178374281624,0.011932401874651296,0.026062761849591155,0.011932401874651296,0.0011451178374281624},
    {0.0025011673089900184,0.026062761849591155,0.056926305563971345,0.026062761849591155,0.0025011673089900184},
    {0.0011451178374281624,0.011932401874651296,0.026062761849591155,0.011932401874651296,0.0011451178374281624},
    {0.00010989362203613103,0.0011451178374281624,0.0025011673089900184,0.0011451178374281624,0.00010989362203613103},
  },
  {
    {1.0546170681737936e-05,0.00010989362203613103,0.00024002973835478245,0.00010989362203613103,1.0546170681737936e-05},
    {0.00010989362203613103,0.0011451178374281624,0.0025011673089900184,0.0011451178374281624,0.00010989362203613103},
    {0.00024002973835478245,0.0025011673089900184,0.005463051664281514,0.0025011673089900184,0.00024002973835478245},
    {0.00010989362203613103,0.0011451178374281624,0.0025011673089900184,0.0011451178374281624,0.00010989362203613103},
    {1.0546170681737936e-05,0.00010989362203613103,0.00024002973835478245,0.00010989362203613103,1.0546170681737936e-05},
  },
};


VectorFieldCalculator::VectorFieldCalculator(ZARR_1 *sf)
{

    // Precompute nearest voxels and associated distance and direction vector
	int i = 0;
    for(int zo = -VF_RAD; zo<=VF_RAD; zo++)
    for(int yo = -VF_RAD; yo<=VF_RAD; yo++)
    for(int xo = -VF_RAD; xo<=VF_RAD; xo++)
      if (xo!=0 || yo!=0 || zo != 0)
	  {
		float dist = sqrt(zo*zo+yo*yo+xo*xo);
        sortedDistances[i][0] = dist;
        sortedDistances[i][1] = xo;
        sortedDistances[i][2] = yo;
        sortedDistances[i][3] = zo;
        sortedDistances[i][4] = xo/dist;
        sortedDistances[i][5] = yo/dist;
        sortedDistances[i][6] = zo/dist;
		i++;
	  }

	// Bubble sort to get them in order 
	int done = 0;
	while(!done)
    {
		done = 1;
		for(int i = 0; i<SORTED_DIST_SIZE-1; i++)
		{
			if (sortedDistances[i][0]>sortedDistances[i+1][0])
			{
				for(int j = 0; j<7; j++)
				{
					float tmp = sortedDistances[i][j];
					sortedDistances[i][j] = sortedDistances[i+1][j];
					sortedDistances[i+1][j] = tmp;
				}
				done = 0;
			}
		}
	}
/*
	for(int i = 0; i<SORTED_DIST_SIZE; i++)
	{	printf("%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",sortedDistances[i][0],sortedDistances[i][1],sortedDistances[i][2],sortedDistances[i][3],sortedDistances[i][4],sortedDistances[i][5],sortedDistances[i][6]);
    }
*/	
	
  surfaceZarr = sf;
}

VectorFieldCalculator::~VectorFieldCalculator()
{
}

void VectorFieldCalculator::GetVectorField(int x, int y, int z, Vec3 &v)
{
	uint64_t p = (((uint64_t)x)<<32)+(((uint64_t)y)<<16)+(uint64_t)z;
	
	extern long n_hit, n_miss;
	// This used to be two table lookups per query, count and then operator[], over seventy
	// million calls per run. find() does one and returns the value already.
	const Vec3 *found = vectorFieldLookup.find(p);
	if (found)
	{
		n_hit++;
		v = *found;
		return;
	}

    uint8_t foundValue = ZARRRead_1(surfaceZarr,z+VOL_OFFSET_Z,y+VOL_OFFSET_Y,x+VOL_OFFSET_X);
	Vec3 vf;
	
	{
        int minCount = 0;
        float lastDist = 0;
		int i=0;
        for(; i<SORTED_DIST_SIZE; i++)
		{ 
            if (ZARRRead_1(surfaceZarr,z+sortedDistances[i][3]+VOL_OFFSET_Z,y+sortedDistances[i][2]+VOL_OFFSET_Y,x+sortedDistances[i][1]+VOL_OFFSET_X) != foundValue)
			{

                if (sortedDistances[i][0] != lastDist)
				{
                    if (minCount>0)
					{
						/* If we're in a surface voxel, then point the vector away from the edge of the surface, but
						   if we're not in a surface voxel point it towards the surface */
                        v = (vf/minCount) * (foundValue==0?1.0f:-1.0f);
						break;
					}
        
					vf.x=vf.y=vf.z=0.0;
					lastDist = sortedDistances[i][0];
				}
                
				vf.x += sortedDistances[i][4];
                vf.y += sortedDistances[i][5];
                vf.z += sortedDistances[i][6];
                minCount += 1;
			}
		}
			
		if (i==SORTED_DIST_SIZE && minCount>0)
		{
			/* If we're in a surface voxel, then point the vector away from the edge of the surface, but
			   if we're not in a surface voxel point it towards the surface */
            v = (vf/minCount) * (foundValue==0?1.0f:-1.0f);
		}
	}
	
	n_miss++;
	vectorFieldLookup.insert(p,v);
}

void VectorFieldCalculator::ComputeField(int x, int y, int z, Vec3 &v, ZARR_1 *reader)
{
	

    uint8_t foundValue = ZARRRead_1(reader,z+VOL_OFFSET_Z,y+VOL_OFFSET_Y,x+VOL_OFFSET_X);
	Vec3 vf;
	
	{
        int minCount = 0;
        float lastDist = 0;
		int i=0;
        for(; i<SORTED_DIST_SIZE; i++)
		{ 
            if (ZARRRead_1(reader,z+sortedDistances[i][3]+VOL_OFFSET_Z,y+sortedDistances[i][2]+VOL_OFFSET_Y,x+sortedDistances[i][1]+VOL_OFFSET_X) != foundValue)
			{

                if (sortedDistances[i][0] != lastDist)
				{
                    if (minCount>0)
					{
						/* If we're in a surface voxel, then point the vector away from the edge of the surface, but
						   if we're not in a surface voxel point it towards the surface */
                        v = (vf/minCount) * (foundValue==0?1.0f:-1.0f);
						break;
					}
        
					vf.x=vf.y=vf.z=0.0;
					lastDist = sortedDistances[i][0];
				}
                
				vf.x += sortedDistances[i][4];
                vf.y += sortedDistances[i][5];
                vf.z += sortedDistances[i][6];
                minCount += 1;
			}
		}
			
		if (i==SORTED_DIST_SIZE && minCount>0)
		{
			/* If we're in a surface voxel, then point the vector away from the edge of the surface, but
			   if we're not in a surface voxel point it towards the surface */
            v = (vf/minCount) * (foundValue==0?1.0f:-1.0f);
		}
	}
	
}


void VectorFieldCalculator::GetSmoothedVectorField(int x, int y, int z, Vec3 &v)
{
	const int L = 2*SMOOTH_WINDOW+1;
	const int N = L*L*L;

	int cx[N], cy[N], cz[N];
	uint64_t key[N];
	float weight[N];
	Vec3 value[N];
	bool missing[N];

	int k = 0;
	for(int zo = z-SMOOTH_WINDOW; zo<=z+SMOOTH_WINDOW; zo++)
	{
	  int zod=zo;
	  if (zo<0) zod = -zo;
	  if (zo>VOL_SIZE_Z-1) zod = 2*VOL_SIZE_Z-zo-1;
	  for(int yo = y-SMOOTH_WINDOW; yo<=y+SMOOTH_WINDOW; yo++)
	  {
	    int yod=yo;
        if (yo<0) yod = -yo;
	    if (yo>VOL_SIZE_Y-1) yod = 2*VOL_SIZE_Y-yo-1;
	    for(int xo = x-SMOOTH_WINDOW; xo<=x+SMOOTH_WINDOW; xo++)
	    {
		  int xod=xo;
	      if (xo<0) xod = -xo;
		  if (xo>VOL_SIZE_X-1) xod = 2*VOL_SIZE_X-xo-1;

		  cx[k]=xo; cy[k]=yo; cz[k]=zo;
		  weight[k] = kernel[zod-z+SMOOTH_WINDOW][yod-y+SMOOTH_WINDOW][xod-x+SMOOTH_WINDOW];
		  key[k] = (((uint64_t)xo)<<32)+(((uint64_t)yo)<<16)+(uint64_t)zo;
		  // The hit counter belongs here, where the lookup actually happens: after the
		  // restructuring it no longer goes through GetVectorField, and leaving it there made
		  // the code report zero hits, which is a false diagnostic.
		  extern long n_hit;
		  const Vec3 *tr = vectorFieldLookup.find(key[k]);
		  if (tr) { n_hit++; value[k] = *tr; missing[k] = false; }
		  else missing[k] = true;
		  k++;
	    }
	  }
	}

	int toDo[N], nToDo = 0;
	for(int i = 0; i<N; i++) if (missing[i]) toDo[nToDo++] = i;

	if (nToDo > 0)
	{
		{
			// Dense block. Every tap the missing points need lies inside a cube of side
			// 2*(SMOOTH_WINDOW+VF_RAD)+1, that is 9 and 729 bytes, which fits in L1. Read in one
			// go it saves about 800 calls to ZARRRead_1 per evaluation, each with three
			// divisions, three moduli and the chunk check: 1.8 billion reads measured over a
			// 300-patch growth.
			// ComputeFieldDense performs exactly the same operations in the same order as
			// ComputeField: only the source of the bytes changes, so the output is identical
			// byte for byte. When the cube straddles chunks, or leaves the volume, the old path
			// is taken, so no voxel is ever read that the original code would not read.
			const int RADIUS = SMOOTH_WINDOW + VF_RAD;
			const int SIDE = 2*RADIUS + 1;
			unsigned char dense[(2*(SMOOTH_WINDOW+VF_RAD)+1)*(2*(SMOOTH_WINDOW+VF_RAD)+1)*(2*(SMOOTH_WINDOW+VF_RAD)+1)];
			int oz = z+VOL_OFFSET_Z-RADIUS, oy = y+VOL_OFFSET_Y-RADIUS, ox = x+VOL_OFFSET_X-RADIUS;
			int inside = (oz >= 0 && oy >= 0 && ox >= 0 &&
			              oz+SIDE <= VOL_SIZE_Z && oy+SIDE <= VOL_SIZE_Y && ox+SIDE <= VOL_SIZE_X);
			if (inside && ZARRReadBlock_1(surfaceZarr,oz,oy,ox,SIDE,SIDE,SIDE,dense))
			{
				for(int j = 0; j<nToDo; j++)
				{
					int i = toDo[j];
					ComputeFieldDense(cx[i],cy[i],cz[i],value[i],dense,oz,oy,ox,SIDE,SIDE,SIDE);
				}
			}
			else
			{
				for(int j = 0; j<nToDo; j++)
				{
					int i = toDo[j];
					ComputeField(cx[i],cy[i],cz[i],value[i],surfaceZarr);
				}
			}
		}
		extern long n_miss;
		for(int j = 0; j<nToDo; j++) { n_miss++; vectorFieldLookup.insert(key[toDo[j]],value[toDo[j]]); }
	}

	// Weighted sum in the original order: left untouched, because floating-point addition is
	// not associative and reordering it would change the result.
	v.x=v.y=v.z=0.0;
	for(int i = 0; i<N; i++)
		v += weight[i]*value[i];

	v = v*GAUSS_SCALE;
}

void VectorFieldCalculator::GetSmoothedVectorFieldInt8(int x, int y, int z, Vec3 &v)
{
	GetSmoothedVectorField(x,y,z,v);

    v = v * 127.0;
}


// Computes the vector field from a DENSE snapshot of the region instead of through the zarr
// reader. The buffer covers absolute coordinates [oz,oz+nz) x [oy,oy+ny) x [ox,ox+nx) and is
// immutable while the computation runs, so several threads may read it without locks and nothing
// is decompressed twice. Decompression is paid ONCE while filling the block; the computation
// itself divides cleanly.
//
// The body is transcribed from the original GetVectorField line by line, including the order of
// the sums, which is left untouched because floating-point addition is not associative.
void VectorFieldCalculator::ComputeFieldDense(int x, int y, int z, Vec3 &v,
                                              const unsigned char *dense, int oz, int oy, int ox,
                                              int nz, int ny, int nx)
{
	// The casts to int are explicit because sortedDistances is an array of floats and in the
	// original the conversion was done by the signature of ZARRRead_1.
	#define DENSE(ZZ,YY,XX) dense[((size_t)((int)(ZZ)-oz)*ny + ((int)(YY)-oy))*nx + ((int)(XX)-ox)]

	uint8_t foundValue = DENSE(z+VOL_OFFSET_Z, y+VOL_OFFSET_Y, x+VOL_OFFSET_X);
	Vec3 vf;

	{
		int minCount = 0;
		float lastDist = 0;
		int i=0;
		for(; i<SORTED_DIST_SIZE; i++)
		{
			if (DENSE(z+sortedDistances[i][3]+VOL_OFFSET_Z,
			          y+sortedDistances[i][2]+VOL_OFFSET_Y,
			          x+sortedDistances[i][1]+VOL_OFFSET_X) != foundValue)
			{
				if (sortedDistances[i][0] != lastDist)
				{
					if (minCount>0)
					{
						v = (vf/minCount) * (foundValue==0?1.0f:-1.0f);
						break;
					}
					vf.x=vf.y=vf.z=0.0;
					lastDist = sortedDistances[i][0];
				}
				vf.x += sortedDistances[i][4];
				vf.y += sortedDistances[i][5];
				vf.z += sortedDistances[i][6];
				minCount += 1;
			}
		}

		if (i==SORTED_DIST_SIZE && minCount>0)
			v = (vf/minCount) * (foundValue==0?1.0f:-1.0f);
	}
	#undef DENSE
}



// Flat table, open addressing, linear probing.
// Key 0 cannot occur, because the coordinates are positive, so it serves as the empty marker.
//
// The memory comes from posix_memalign on a 2 MiB boundary and is offered to the kernel as a huge
// page. Transparent huge pages are in "madvise" mode on this machine, so without the advice the
// table is paged in 4 KiB pieces and every probe that misses the TLB walks the page tables.
// A zeroed slot is an empty slot: the key is 0 and Vec3's default constructor zeroes its three
// floats, so zeroed memory is exactly what assign(FieldSlot{0,Vec3()}) produced.
#include <sys/mman.h>
#include <cstdlib>
#include <cstring>

static int fieldCacheHugePages(void)
{
	static int decided = -1;
	if (decided < 0)
	{
		const char *e = getenv("SIMPAPER_HUGE_PAGES");
		decided = (e && atoi(e) == 0) ? 0 : 1;
	}
	return decided;
}

VectorFieldCalculator::FieldSlot *VectorFieldCalculator::FieldCache::allocSlots(size_t n)
{
	size_t bytes = n * sizeof(FieldSlot);
	void *p = NULL;
	const size_t HUGE = 2u*1024u*1024u;
	if (fieldCacheHugePages() && bytes >= HUGE)
	{
		size_t rounded = ((bytes + HUGE - 1) / HUGE) * HUGE;
		if (posix_memalign(&p, HUGE, rounded) != 0) p = NULL;
		if (p)
		{
			madvise(p, rounded, MADV_HUGEPAGE);
			memset(p, 0, rounded);
			return (FieldSlot *)p;
		}
	}
	p = calloc(n, sizeof(FieldSlot));
	if (!p) { fprintf(stderr,"FieldCache: out of memory for %zu slots\n", n); exit(5); }
	return (FieldSlot *)p;
}

void VectorFieldCalculator::FieldCache::freeSlots(FieldSlot *p, size_t)
{
	free(p);
}

VectorFieldCalculator::FieldCache::FieldCache()
{
	cap = 1024;
	tab = allocSlots(cap);
	mask = cap - 1;
	used = 0;
}

VectorFieldCalculator::FieldCache::~FieldCache()
{
	freeSlots(tab, cap);
	tab = nullptr;
}

void VectorFieldCalculator::FieldCache::clear()
{
	freeSlots(tab, cap);
	cap = 1024;
	tab = allocSlots(cap);
	mask = cap - 1;
	used = 0;
}

void VectorFieldCalculator::FieldCache::grow()
{
	FieldSlot *old = tab;
	size_t oldCap = cap;
	cap = (oldCap ? oldCap : 1024) * 2;
	tab = allocSlots(cap);
	mask = cap - 1;
	// Same order as before: the old table is walked from slot 0 upwards, so two entries that
	// collide in the new table land in the same relative order they did with std::vector.
	for(size_t j = 0; j < oldCap; j++)
	{
		const FieldSlot &s = old[j];
		if (s.k)
		{
			size_t i = (size_t)(s.k * 0x9E3779B97F4A7C15ull >> 32) & mask;
			while (tab[i].k) i = (i+1) & mask;
			tab[i] = s;
		}
	}
	freeSlots(old, oldCap);
}

const Vec3 *VectorFieldCalculator::FieldCache::find(uint64_t k) const
{
	size_t i = (size_t)(k * 0x9E3779B97F4A7C15ull >> 32) & mask;
	while (tab[i].k)
	{
		if (tab[i].k == k) return &tab[i].v;
		i = (i+1) & mask;
	}
	return NULL;
}

void VectorFieldCalculator::FieldCache::insert(uint64_t k, const Vec3 &v)
{
	if ((used+1)*10 >= cap*7) grow();
	size_t i = (size_t)(k * 0x9E3779B97F4A7C15ull >> 32) & mask;
	while (tab[i].k)
	{
		if (tab[i].k == k) { tab[i].v = v; return; }
		i = (i+1) & mask;
	}
	tab[i].k = k; tab[i].v = v; used++;
}
