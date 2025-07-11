MODULE_big = pg_kmersearch
OBJS = pg_kmersearch.o kmersearch_gin.o kmersearch_datatype.o kmersearch_kmer.o

EXTENSION = pg_kmersearch
DATA = pg_kmersearch--1.0.sql
PGFILEDESC = "pg_kmersearch - k-mer search for DNA sequences"

REGRESS = 01_basic_types 02_configuration 03_tables_indexes 04_search_operators 05_scoring_functions 06_advanced_search 07_length_functions 08_cache_management 09_highfreq_filter 10_parallel_cache 11_cache_hierarchy

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

override CPPFLAGS += -Wno-unused-variable -Wno-unused-function