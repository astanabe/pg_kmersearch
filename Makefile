MODULE_big = pg_kmersearch
OBJS = kmersearch.o kmersearch_gin.o kmersearch_datatype.o kmersearch_kmer.o kmersearch_cache.o kmersearch_freq.o kmersearch_partition.o

EXTENSION = pg_kmersearch
DATA = pg_kmersearch--1.0.sql
PGFILEDESC = "pg_kmersearch - k-mer search for DNA sequences"

REGRESS = 01_basic_types 02_configuration 03_tables_indexes 04_search_operators 05_scoring_functions 06_advanced_search 07_length_functions 08_cache_management 09_highfreq_filter 10_parallel_cache 11_cache_hierarchy 12_management_views 13_partition_functions

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

override CPPFLAGS += -Wno-unused-variable -Wno-unused-function -std=c99
override CFLAGS := $(filter-out -Werror=vla, $(CFLAGS))

# SIMD Optimization Support
# To enable SIMD optimizations for DNA comparison functions, uncomment the
# appropriate lines below based on your target CPU architecture:
#
# For x86_64 processors:
#   AVX2 support:         override CPPFLAGS += -mavx2
#   AVX512F support:      override CPPFLAGS += -mavx512f
#   AVX512BW support:     override CPPFLAGS += -mavx512f -mavx512bw
#   All optimizations:    override CPPFLAGS += -mavx2 -mavx512f -mavx512bw
#
# For ARM64 processors:
#   NEON support:      override CPPFLAGS += -mfpu=neon (automatically enabled on ARM64)
#   SVE support:       override CPPFLAGS += -march=armv8-a+sve
#
# Example for enabling AVX2 on x86_64:
# override CPPFLAGS += -mavx2
#
# Example for enabling AVX512F (basic set) on x86_64:
# override CPPFLAGS += -mavx512f
#
# Example for enabling full AVX512BW support on x86_64:
# override CPPFLAGS += -mavx2 -mavx512f -mavx512bw
#
# Note: Ensure your target CPU supports these instruction sets before enabling