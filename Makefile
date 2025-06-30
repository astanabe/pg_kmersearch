# pg_kmersearch Makefile

EXTENSION = pg_kmersearch
DATA = pg_kmersearch--1.0.sql
MODULES = pg_kmersearch

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

# Additional compiler flags
PG_CFLAGS = $(shell $(PG_CONFIG) --cflags)
PG_CPPFLAGS = -I$(shell $(PG_CONFIG) --includedir)

SHLIB_LINK += $(shell $(PG_CONFIG) --libs)

# Add our custom flags
override CFLAGS += $(PG_CFLAGS) $(PG_CPPFLAGS) -std=c99 -Wall -Wextra