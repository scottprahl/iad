#include "cache.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---------- low-level hash utilities ------------------------------------ */
static inline uint64_t dbl_bits(double x)
{
    uint64_t u;
    memcpy(&u, &x, sizeof u);        /* exact bit pattern */
    return u;
}

static uint64_t hash_key(double a, double b, double g)
{
    const uint64_t FNV_OFFSET = 0xcbf29ce484222325ULL;
    const uint64_t FNV_PRIME  = 0x100000001b3ULL;

    uint64_t h = FNV_OFFSET;
    #define STEP(v) do { h ^= (v); h *= FNV_PRIME; } while (0)
    STEP(dbl_bits(a)); STEP(dbl_bits(b)); STEP(dbl_bits(g));
    return h;
}

/* ---------- table layout ------------------------------------------------- */
typedef enum { SLOT_EMPTY, SLOT_FULL, SLOT_TOMB } slot_state;

typedef struct {
    slot_state   state;
    uint64_t     hash;
    CacheRecord  rec;
} Slot;

struct Cache {
    Slot  *tbl;
    size_t cap;        /* power of two */
    size_t sz;
    size_t tombs;
};

#define LOAD_NUM 7
#define LOAD_DEN 10

/* forward decls ---------------------------------------------------- */
static size_t lookup_slot_const(const Cache *c,
                                double a, double b, double g,
                                uint64_t h, int *found);

static size_t lookup_slot(Cache *c,
                          double a, double b, double g,
                          uint64_t h, int *found);

static void   resize(Cache *c, size_t new_cap);

/* ---------- public functions --------------------------------------------- */
int cache_init(Cache *c)
{
    if (!c) return -1;                    /* invalid pointer */

    /* If the caller re-uses a handle without calling cache_free first,
       release the old table to prevent a leak.                         */
    if (c->tbl) {
        free(c->tbl);
        c->tbl = NULL;
    }

    c->cap = 8;                           /* power-of-two ≥ 4 */
    c->sz = c->tombs = 0;

    c->tbl = calloc(c->cap, sizeof(Slot));
    if (!c->tbl) {                        /* allocation failed */
        c->cap = 0;
        return -1;                        /* propagate OOM to caller */
    }
    return 0;                             /* success */
}

void cache_free(Cache *c)
{
    if (!c) return;                       /* graceful on NULL handle */

    if (c->tbl) {
        free(c->tbl);
        c->tbl = NULL;
    }
    c->cap = c->sz = c->tombs = 0;        /* leave in "empty" state */
}

int cache_get(const Cache *c,
              double a, double b, double g,
              double *ur1, double *ut1, double *uru, double *utu)
{
    uint64_t h = hash_key(a,b,g);
    int found;
    size_t pos = lookup_slot_const(c, a, b, g, h, &found);

    if (!found) return 0;

    const CacheRecord *r = &c->tbl[pos].rec;
    if (ur1) *ur1 = r->ur1; if (ut1) *ut1 = r->ut1;
    if (uru) *uru = r->uru; if (utu) *utu = r->utu;
    return 1;
}

void cache_put(Cache *c,
               double a, double b, double g,
               double ur1, double ut1, double uru, double utu)
{
    uint64_t h = hash_key(a,b,g);
    if ((c->sz + c->tombs) * LOAD_DEN > c->cap * LOAD_NUM)
        resize(c, c->cap << 1);

    int found;
    size_t pos = lookup_slot(c, a, b, g, h, &found);
    Slot *s = &c->tbl[pos];

    if (found) {                      /* overwrite */
        s->rec.ur1 = ur1; s->rec.ut1 = ut1;
        s->rec.uru = uru; s->rec.utu = utu;
        return;
    }

    if (s->state == SLOT_TOMB) c->tombs--;
    s->state = SLOT_FULL;
    s->hash  = h;
    s->rec.a = a; s->rec.b = b; s->rec.g = g;
    s->rec.ur1 = ur1; s->rec.ut1 = ut1;
    s->rec.uru = uru; s->rec.utu = utu;
    c->sz++;
}

/* ---------- iteration ---------------------------------------------------- */
size_t cache_size(const Cache *c) { return c->sz; }

int cache_get_by_index(const Cache *c, size_t i,
                       double *a, double *b, double *g,
                       double *ur1, double *ut1, double *uru, double *utu)
{
    if (i >= c->sz) return 0;
    size_t seen = 0;
    for (size_t p = 0; p < c->cap; ++p)
        if (c->tbl[p].state == SLOT_FULL && seen++ == i) {
            const CacheRecord *r = &c->tbl[p].rec;
            if (a) *a = r->a; if (b) *b = r->b; if (g) *g = r->g;
            if (ur1) *ur1 = r->ur1; if (ut1) *ut1 = r->ut1;
            if (uru) *uru = r->uru; if (utu) *utu = r->utu;
            return 1;
        }
    return 0;
}

/* ---------- internals ---------------------------------------------------- */
static int key_equal(const CacheRecord *r,
                     double a, double b, double g)
{ return r->a == a && r->b == b && r->g == g; }

/* read-only probing – never modifies the table -------------------- */
static size_t lookup_slot_const(const Cache *c,
                                double a, double b, double g,
                                uint64_t h, int *found)
{
    size_t mask = c->cap - 1, tomb = SIZE_MAX, pos = h & mask;

    for (;;) {
        const Slot *s = &c->tbl[pos];                /* note the const */
        if (s->state == SLOT_EMPTY) {
            *found = 0;
            return (tomb != SIZE_MAX) ? tomb : pos;
        }
        if (s->state == SLOT_TOMB) {
            if (tomb == SIZE_MAX) tomb = pos;
        } else if (s->hash == h &&      /* key match */
                   s->rec.a == a && s->rec.b == b && s->rec.g == g) {
            *found = 1;
            return pos;
        }
        pos = (pos + 1) & mask;                    /* linear probe */
    }
}

static size_t lookup_slot(Cache *c,
                          double a, double b, double g,
                          uint64_t h, int *found)
{
    return lookup_slot_const((const Cache *)c, a, b, g, h, found);
}

static void resize(Cache *c, size_t new_cap)
{
    Slot *old = c->tbl; size_t old_cap = c->cap;
    c->tbl = calloc(new_cap, sizeof(Slot));
    if (!c->tbl) abort();
    c->cap = new_cap; c->sz = c->tombs = 0;
    for (size_t i=0;i<old_cap;++i)
        if (old[i].state == SLOT_FULL) {
            CacheRecord *r = &old[i].rec;
            cache_put(c, r->a,r->b,r->g, r->ur1,r->ut1,r->uru,r->utu);
        }
    free(old);
}
