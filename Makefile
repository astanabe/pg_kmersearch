# pg_kmersearch Makefile

EXTENSION = pg_kmersearch
DATA = pg_kmersearch--1.0.sql
MODULES = pg_kmersearch
REGRESS = 01_basic_types 02_configuration 03_tables_indexes 04_search_operators 05_scoring_functions 06_advanced_search

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

# Additional compiler flags
PG_CFLAGS = $(shell $(PG_CONFIG) --cflags)
PG_CPPFLAGS = -I$(shell $(PG_CONFIG) --includedir)

SHLIB_LINK += $(shell $(PG_CONFIG) --libs)

# Add our custom flags
override CFLAGS += $(PG_CFLAGS) $(PG_CPPFLAGS) -std=c99 -Wall -Wextra