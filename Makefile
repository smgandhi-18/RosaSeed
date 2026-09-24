# /*************************************************************************************
#                            The MIT License

#    BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
#    Copyright (C) 2019  Intel Corporation, Heng Li.

#    Permission is hereby granted, free of charge, to any person obtaining
#    a copy of this software and associated documentation files (the
#    "Software"), to deal in the Software without restriction, including
#    without limitation the rights to use, copy, modify, merge, publish,
#    distribute, sublicense, and/or sell copies of the Software, and to
#    permit persons to whom the Software is furnished to do so, subject to
#    the following conditions:

#    The above copyright notice and this permission notice shall be
#    included in all copies or substantial portions of the Software.

#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

# Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
# Authors: Gandhi Shyama <smgandhi@ualberta.ca> [RosaSeed code addition]

# *****************************************************************************************/

ifneq ($(portable),)
	STATIC_GCC=-static-libgcc -static-libstdc++
endif

EXE=		bwa-mem2

ifeq ($(CXX), icpc)
	CC=icc
else ifeq ($(CXX), g++)
	CC="gcc -Wno-error=implicit-function-declaration"
else
	CC="gcc -Wno-error=implicit-function-declaration"
	CXX?=g++
endif

ARCH_FLAGS=	-msse -msse2 -msse3 -mssse3 -msse4.1
MEM_FLAGS=	-DSAIS=1

CPPFLAGS+=	-DENABLE_PREFETCH -DV17=1 -DMATE_SORT=0 $(MEM_FLAGS) -D_GNU_SOURCE
CPPFLAGS+=	$(CPPFLAGS_EXTRA)

INCLUDES=	-Isrc -Iext/safestringlib/include

LIBS=		-lpthread -lm -lz -L. -lbwa -Lext/safestringlib -lsafestring $(STATIC_GCC)

# ============================================================
# RosaSeed selection
#
# Normal BWA-MEM2:
#   make arch=native CXX=g++
#
# 2-step RosaSeed:
#   make arch=native CXX=g++ ROSASEED=1
#
# RosaSeed-Compact (1-step):
#   make arch=native CXX=g++ ROSASEED=1 ROSASEED_1STEP=1
#
# Extra flags:
#   CPPFLAGS_EXTRA="-DSA_COMPRESSION_FACTOR_POWER=2"
#   CPPFLAGS_EXTRA="-DSA_COMPRESSION_FACTOR_POWER=2 -DENABLE_PHASE_II"
#   CPPFLAGS_EXTRA="-DSA_COMPRESSION_FACTOR_POWER=2 -DGAPFILL_ALWAYS_RUN_BOTH_STRANDS -DGAPFILL_EARLY_EXIT"
# ============================================================

ifdef ROSASEED
	CPPFLAGS += -DROSASEED_INPROCESS
	ifdef ROSASEED_1STEP
		CPPFLAGS += -DROSASEED_1STEP
		ROSA_INC = -Isrc/rosaseed_compact
		ROSA_OBJS = \
			src/rosaseed_inprocess_compact.o \
			src/rosaseed_core_bridge_compact.o \
			src/rosaseed_compact/mem_alloc_compact.o \
			src/rosaseed_compact/load_data_compact_mmap.o \
			src/rosaseed_compact/helper_functions.o \
			src/rosaseed_compact/profiling.o \
			src/rosaseed_compact/bwa.o \
			src/rosaseed_compact/rosaseed_compact_phaseA.o \
			src/rosaseed_compact/rosaseed_compact_gapfill_PhaseC.o \
			src/rosaseed_compact/rosaseed_compact_phaseB.o
	else
		ROSA_INC = -Isrc/rosaseed
		ROSA_OBJS = \
			src/rosaseed_inprocess.o \
			src/rosaseed_core_bridge.o \
			src/rosaseed/read_init.o \
			src/rosaseed/mem_alloc.o \
			src/rosaseed/load_data_mmap.o \
			src/rosaseed/helper_functions.o \
			src/rosaseed/profiling.o \
			src/rosaseed/bwa.o \
			src/rosaseed/rosaseed_phaseA.o \
			src/rosaseed/rosaseed_gapfill_PhaseC.o \
			src/rosaseed/rosaseed_phaseB.o
	endif
else
	ROSA_INC =
	ROSA_OBJS =
endif

OBJS=		src/fastmap.o src/bwtindex.o src/utils.o src/memcpy_bwamem.o src/kthread.o \
			src/kstring.o src/ksw.o src/bntseq.o src/bwamem.o src/profiling.o src/bandedSWA.o \
			src/FMI_search.o src/read_index_ele.o src/bwamem_pair.o src/kswv.o src/bwa.o \
			src/bwamem_extra.o src/kopen.o \
			$(ROSA_OBJS)

BWA_LIB=	libbwa.a
SAFE_STR_LIB=	ext/safestringlib/libsafestring.a

ifeq ($(arch),sse41)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-msse4.1
	else
		ARCH_FLAGS=-msse -msse2 -msse3 -mssse3 -msse4.1
	endif
else ifeq ($(arch),sse42)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-msse4.2
	else
		ARCH_FLAGS=-msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2
	endif
else ifeq ($(arch),avx)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-mavx
	else
		ARCH_FLAGS=-mavx
	endif
else ifeq ($(arch),avx2)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-march=core-avx2
	else
		ARCH_FLAGS=-mavx2
	endif
else ifeq ($(arch),avx512)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512
	else
		ARCH_FLAGS=-mavx512bw
	endif
else ifeq ($(arch),native)
	ARCH_FLAGS=-march=native
else ifneq ($(arch),)
	ARCH_FLAGS=$(arch)
else
myall:multi
endif

CXXFLAGS+=	-g -O3 -fpermissive $(ARCH_FLAGS)

.PHONY: all clean depend multi
.SUFFIXES: .cpp .c .o

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $(ROSA_INC) $< -o $@

.c.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $(ROSA_INC) $< -o $@

all: $(EXE)

multi:
	rm -f src/*.o src/rosaseed/*.o src/rosaseed_compact/*.o  $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse41 EXE=bwa-mem2.sse41 CXX=$(CXX) ROSASEED=$(ROSASEED) ROSASEED_1STEP=$(ROSASEED_1STEP) CPPFLAGS_EXTRA="$(CPPFLAGS_EXTRA)" all
	rm -f src/*.o src/rosaseed/*.o src/rosaseed_compact/*.o  $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse42 EXE=bwa-mem2.sse42 CXX=$(CXX) ROSASEED=$(ROSASEED) ROSASEED_1STEP=$(ROSASEED_1STEP) CPPFLAGS_EXTRA="$(CPPFLAGS_EXTRA)" all
	rm -f src/*.o src/rosaseed/*.o src/rosaseed_compact/*.o  $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx EXE=bwa-mem2.avx CXX=$(CXX) ROSASEED=$(ROSASEED) ROSASEED_1STEP=$(ROSASEED_1STEP) CPPFLAGS_EXTRA="$(CPPFLAGS_EXTRA)" all
	rm -f src/*.o src/rosaseed/*.o src/rosaseed_compact/*.o  $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx2 EXE=bwa-mem2.avx2 CXX=$(CXX) ROSASEED=$(ROSASEED) ROSASEED_1STEP=$(ROSASEED_1STEP) CPPFLAGS_EXTRA="$(CPPFLAGS_EXTRA)" all
	rm -f src/*.o src/rosaseed/*.o src/rosaseed_compact/*.o  $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx512 EXE=bwa-mem2.avx512bw CXX=$(CXX) ROSASEED=$(ROSASEED) ROSASEED_1STEP=$(ROSASEED_1STEP) CPPFLAGS_EXTRA="$(CPPFLAGS_EXTRA)" all
	$(CXX) -Wall -O3 src/runsimd.cpp -Iext/safestringlib/include -Lext/safestringlib/ -lsafestring $(STATIC_GCC) -o bwa-mem2

$(EXE): $(BWA_LIB) $(SAFE_STR_LIB) src/main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/main.o $(BWA_LIB) $(LIBS) -o $@

$(BWA_LIB): $(OBJS)
	ar rcs $(BWA_LIB) $(OBJS)

$(SAFE_STR_LIB):
	cd ext/safestringlib/ && $(MAKE) clean && $(MAKE) CC=$(CC) directories libsafestring.a

clean:
	rm -fr src/*.o src/rosaseed/*.o src/rosaseed_compact/*.o $(BWA_LIB) $(EXE) \
		bwa-mem2.sse41 bwa-mem2.sse42 bwa-mem2.avx bwa-mem2.avx2 bwa-mem2.avx512bw
	cd ext/safestringlib/ && $(MAKE) clean

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) -I. -- src/*.cpp)
