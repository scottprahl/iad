#ifndef CACHE_H
#define CACHE_H

#include <stddef.h>

typedef struct {
    double a, b, g;
    double ur1, ut1, uru, utu;
} CacheRecord;

typedef struct Cache Cache;

int cache_init(Cache * c);

void cache_free(Cache * c);

int cache_get(const Cache * c, double a, double b, double g, double *ur1, double *ut1, double *uru, double *utu);

void cache_put(Cache * c, double a, double b, double g, double ur1, double ut1, double uru, double utu);

size_t cache_size(const Cache * c);

int cache_get_by_index(const Cache * c, size_t i,
    double *a, double *b, double *g, double *ur1, double *ut1, double *uru, double *utu);

#endif
