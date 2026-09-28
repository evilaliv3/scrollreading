#pragma once

#include <cstdio>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>

#include "zarr_1.h"
#include "vec3.h"
#include "parameters.h"

// Calculate a single vector field point

#define VF_RAD 2
#define SORTED_DIST_SIZE ((VF_RAD*2+1)*(VF_RAD*2+1)*(VF_RAD*2+1)-1)

#define SMOOTH_WINDOW 2
#define GAUSS_SCALE 1

class VectorFieldCalculator
{
	public:
		VectorFieldCalculator(ZARR_1 *sf);
		~VectorFieldCalculator();
		
		void GetVectorField(int x, int y, int z, Vec3 &v);
		// Pure version: computes and nothing else. It neither reads nor writes the cache, and
		// it reads through the reader it is given, so several voxels can be computed at once,
		// each with its own reader and no shared state.
		void ComputeField(int x, int y, int z, Vec3 &v, ZARR_1 *reader);
		// Computes the field from a dense, immutable snapshot of the region instead of through
		// the zarr reader: no state, no lock, and the same buffer can serve several callers.
		// The arithmetic is identical to ComputeField, operation for operation and in the same
		// order, so the output is identical byte for byte. All the taps one smoothing evaluation
		// needs fit in a 9x9x9 cube, 729 bytes, which sits in L1.
		void ComputeFieldDense(int x, int y, int z, Vec3 &v,
		                       const unsigned char *dense, int oz, int oy, int ox,
		                       int nz, int ny, int nx);

        void GetSmoothedVectorField(int x, int y, int z, Vec3 &v);
		void GetSmoothedVectorFieldInt8(int x, int y, int z, Vec3 &v);
		

		ZARR_1 *surfaceZarr;
		float sortedDistances[SORTED_DIST_SIZE][7];

		// Cache of the field taps. It used to be a std::unordered_map: O(1), but every lookup
		// reads the bucket vector and then CHASES a pointer to a separately allocated node, that
		// is two cold memory accesses, seventy million lookups per run over twenty-three million
		// entries. Here it is a flat open-addressed table: one cache line per lookup, no
		// per-entry allocation, and less memory.
		// Dense paging by coordinate would be worse: the voxels touched form a surface, not a
		// solid volume, so nine pages out of ten would stay empty.
		struct FieldSlot { uint64_t k; Vec3 v; };
		class FieldCache
		{
			std::vector<FieldSlot> tab;
			size_t mask = 0, used = 0;
			void grow();
		public:
			FieldCache();
			const Vec3 *find(uint64_t k) const;
			void insert(uint64_t k, const Vec3 &v);
			size_t size() const { return used; }
			void clear() { tab.assign(1024, FieldSlot{0,Vec3()}); mask = 1023; used = 0; }
		};
		FieldCache vectorFieldLookup;

};

