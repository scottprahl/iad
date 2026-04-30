@** IAD Adaptive Grid.

This module implements an adaptive quadtree grid for the two-parameter
inverse searches |FIND_AB|, |FIND_AG|, and |FIND_BG|.  It replaces the
dense $101\times101$ grid with a smaller, targeted sample that still
seeds the Nelder-Mead simplex with a good warm start.

The adaptive grid works by recursively subdividing a
$[u_0,u_1]\times[v_0,v_1]$ rectangle in a nonlinear coordinate space
until the bilinear interpolation error at the cell centre falls below a
tolerance |AGRID_TOL|.  Subdivision is always performed to at least
|AGRID_MIN_DEPTH| and never beyond |AGRID_MAX_DEPTH|.

The coordinate transforms match those used by the dense grid:
$$a(u) = 1-(1-u)^2(1+2u), \qquad u\in[0,1]$$
$$g(v) = \bigl(1-2(1-v)^2(1+2v)\bigr)\cdot G_{\max}, \qquad v\in[0,1]$$
$$b = e^\ell, \qquad \ell\in[-8,\,\pm10]$$

Two helper functions, |abg_eval| and |abg_stored_distance|, live in
\.{iad\_calc.c} (where the global |MM|/|RR| state is accessible) and
are declared in \.{iad\_calc.h}.  Everything else is local to this file.

@ All output goes into \.{iad\_agrid.c}.

@(iad_agrid.c@>=

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ad_globl.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"
#include "iad_agrid.h"

@<AGrid constants@>@;
@<AGrid cache storage@>@;
@<AGrid stale-context storage@>@;

    @<Definition for |agrid_nonlinear_a|@>@;
    @<Definition for |agrid_nonlinear_g|@>@;
    @<Definition for |agrid_make_abg|@>@;
    @<Definition for |agrid_ensure_capacity|@>@;
    @<Definition for |agrid_find_entry|@>@;
    @<Definition for |agrid_add_or_get|@>@;
    @<Definition for |agrid_raw_interp_error|@>@;
    @<Definition for |agrid_subdivide|@>@;
    @<Definition for |agrid_save_context|@>@;
    @<Definition for |AGrid_Valid|@>@;
    @<Definition for |AGrid_Build|@>@;
    @<Definition for |AGrid_Fill_Guesses|@>@;

@ The header file exposes only the three public entry points.

@(iad_agrid.h@>=
    @<Prototype for |AGrid_Build|@>;
    @<Prototype for |AGrid_Valid|@>;
    @<Prototype for |AGrid_Fill_Guesses|@>;

@*1 Constants.

@<AGrid constants@>=
#define AGRID_TOL           0.05    /* subdivision threshold (sphere-corrected M_R+M_T space) */
#define AGRID_MIN_DEPTH     2       /* always subdivide to this depth */
#define AGRID_MAX_DEPTH     4       /* never go deeper */
#define AGRID_INIT_CAP      512     /* initial cache capacity */
#define AGRID_CLOSURE_N     4       /* top-N entries used in local closure */
#define AGRID_MIN_LOG_B     (-8.0)
#define AGRID_MAX_LOG_B_AB  8.0     /* find_ab, find_ag use this upper b limit */
#define AGRID_MAX_LOG_B_BG  10.0    /* find_bg uses a wider b range */

@*1 Cache.

Each entry stores the physical parameters and the raw adding-doubling
results.  The distance is NOT cached---it is recomputed at query time
from the stored RT values so that sphere-parameter changes (which happen
between Monte Carlo iterations) are reflected without rebuilding the RT
table.

@<AGrid cache storage@>=

typedef struct {
    double a, b, g;
    double ur1, ut1, uru, utu;
} agrid_entry_t;

static agrid_entry_t *AGrid_Cache = NULL;
static int AGrid_N = 0;     /* number of valid entries */
static int AGrid_Cap = 0;   /* allocated capacity */
static int AGrid_Search = -1;   /* search mode at last build */

@*1 Stale-context.

We store the key fields that determine whether a rebuild is needed.
These mirror the checks in |Valid_Grid| (|iad_calc.c|) plus the
search-specific fixed-axis value.

@<AGrid stale-context storage@>=

static int    AGrid_Initialized    = 0;
static double AGrid_slab_n         = 0;
static double AGrid_slab_n_top     = 0;
static double AGrid_slab_n_bottom  = 0;
static double AGrid_slab_b_top     = 0;
static double AGrid_slab_b_bottom  = 0;
static double AGrid_slab_cos_angle = 0;
static int    AGrid_num_spheres    = -1;
static int    AGrid_num_measures   = -1;
static double AGrid_m_u            = 0;
/* Lost-light variables removed: subdivision uses include_lost=False */
static double AGrid_fixed_g        = 0;  /* FIND_AB: fixed anisotropy */
static double AGrid_fixed_b        = 0;  /* FIND_AG: fixed thickness */
static double AGrid_fixed_a        = 0;  /* FIND_BG: fixed albedo */

@*1 Coordinate transforms.

The same nonlinear map used by the dense grid concentrates sample
points near the boundaries of the albedo and anisotropy ranges where
the inverse problem is most sensitive.

@<Definition for |agrid_nonlinear_a|@>=
static double agrid_nonlinear_a(double u)
{
    /* u in [0,1] -> a in [0,1], denser near 0 and 1 */
    return 1.0 - (1.0 - u) * (1.0 - u) * (1.0 + 2.0 * u);
}

@ @<Definition for |agrid_nonlinear_g|@>=
static double agrid_nonlinear_g(double v)
{
    /* v in [0,1] -> g in [-MAX_ABS_G, +MAX_ABS_G] */
    double vv = (1.0 - v) * (1.0 - v) * (1.0 + 2.0 * v);
    return (1.0 - 2.0 * vv) * MAX_ABS_G;
}

@ Convert the 2D adaptive coordinates $(u,v)$ to physical $(a,b,g)$.
The mapping depends on which axis is held fixed (i.e.\ which search is
active).

@<Definition for |agrid_make_abg|@>=
static void agrid_make_abg(double u, double v, int search,
                           double *a, double *b, double *g)
{
    switch (search) {
    case FIND_AB:
        /* u = a-coordinate [0,1], v = log_b in [MIN_LOG_B, MAX_LOG_B_AB] */
        *a = agrid_nonlinear_a(u);
        *b = exp(v);
        *g = AGrid_fixed_g;
        break;
    case FIND_AG:
        /* u = a-coordinate [0,1], v = g-coordinate [0,1] */
        *a = agrid_nonlinear_a(u);
        *b = AGrid_fixed_b;
        *g = agrid_nonlinear_g(v);
        break;
    case FIND_BG:
        /* u = log_b in [MIN_LOG_B, MAX_LOG_B_BG], v = g-coordinate [0,1] */
        *a = AGrid_fixed_a;
        *b = exp(u);
        *g = agrid_nonlinear_g(v);
        break;
    default:
        *a = 0.5; *b = 1.0; *g = 0.0;
        break;
    }

    /* clamp to physical range */
    if (*b < 1e-8)       *b = 1e-8;
    if (*g < -MAX_ABS_G) *g = -MAX_ABS_G;
    if (*g >  MAX_ABS_G) *g =  MAX_ABS_G;
}

@*1 Cache management.

@<Definition for |agrid_ensure_capacity|@>=
static void agrid_ensure_capacity(void)
{
    if (AGrid_N < AGrid_Cap) return;
    AGrid_Cap = (AGrid_Cap == 0) ? AGRID_INIT_CAP : AGrid_Cap * 2;
    AGrid_Cache = (agrid_entry_t *) realloc(AGrid_Cache,
                      (size_t) AGrid_Cap * sizeof(agrid_entry_t));
    if (AGrid_Cache == NULL) {
        fprintf(stderr, "AGrid: out of memory\n");
        exit(EXIT_FAILURE);
    }
}

@ Linear scan for an existing cache entry.  Returns the index if found,
$-1$ otherwise.  Exact floating-point comparison is intentional: we
only generate canonical $(u,v)$ pairs from the subdivision recursion,
so duplicates arise from shared cell corners.

@<Definition for |agrid_find_entry|@>=
static int agrid_find_entry(double a, double b, double g)
{
    int i;
    for (i = 0; i < AGrid_N; i++) {
        if (AGrid_Cache[i].a == a &&
            AGrid_Cache[i].b == b &&
            AGrid_Cache[i].g == g)
            return i;
    }
    return -1;
}

@ Evaluate the adding-doubling RT at $(a,b,g)$ if not already cached,
then return the index.  Uses |abg_eval| from \.{iad\_calc.c} which has
access to the global |MM|/|RR| state.

@<Definition for |agrid_add_or_get|@>=
static int agrid_add_or_get(double a, double b, double g)
{
    int idx = agrid_find_entry(a, b, g);
    if (idx >= 0) return idx;

    agrid_ensure_capacity();
    idx = AGrid_N;
    AGrid_Cache[idx].a = a;
    AGrid_Cache[idx].b = b;
    AGrid_Cache[idx].g = g;
    abg_eval(a, b, g,
             &AGrid_Cache[idx].ur1,
             &AGrid_Cache[idx].ut1,
             &AGrid_Cache[idx].uru,
             &AGrid_Cache[idx].utu);
    AGrid_N++;
    return idx;
}

@*1 Subdivision error metric.

The subdivision uses sphere-corrected predicted $(M_R, M_T)$ values (without
lost-light corrections, matching the Python |include\_lost=False| convention)
as the interpolation-error estimate.  Sphere corrections can amplify small
$(ur1, ut1)$ differences into large $(M_R, M_T)$ differences in the
high-albedo region, so a criterion in corrected space naturally places finer
grid points where they are needed.

The corrected sphere distance to the actual measurements is used when ranking
candidates at query time (see |AGrid_Fill_Guesses|).

@<Definition for |agrid_raw_interp_error|@>=
static double agrid_raw_interp_error(int c00, int c10, int c01, int c11, int cc)
{
    double mr00, mt00, mr10, mt10, mr01, mt01, mr11, mt11, mrcc, mtcc;
    double interp_mr, interp_mt;
    agrid_entry_t *e;

    e = &AGrid_Cache[c00];
    abg_sphere_mr_mt(e->a, e->b, e->g, e->ur1, e->ut1, e->uru, e->utu, &mr00, &mt00);
    e = &AGrid_Cache[c10];
    abg_sphere_mr_mt(e->a, e->b, e->g, e->ur1, e->ut1, e->uru, e->utu, &mr10, &mt10);
    e = &AGrid_Cache[c01];
    abg_sphere_mr_mt(e->a, e->b, e->g, e->ur1, e->ut1, e->uru, e->utu, &mr01, &mt01);
    e = &AGrid_Cache[c11];
    abg_sphere_mr_mt(e->a, e->b, e->g, e->ur1, e->ut1, e->uru, e->utu, &mr11, &mt11);
    e = &AGrid_Cache[cc];
    abg_sphere_mr_mt(e->a, e->b, e->g, e->ur1, e->ut1, e->uru, e->utu, &mrcc, &mtcc);

    interp_mr = 0.25 * (mr00 + mr10 + mr01 + mr11);
    interp_mt = 0.25 * (mt00 + mt10 + mt01 + mt11);
    return fabs(mrcc - interp_mr) + fabs(mtcc - interp_mt);
}

@*1 Recursive subdivision.

@<Definition for |agrid_subdivide|@>=
static void agrid_subdivide(double u0, double u1, double v0, double v1,
                            int depth, int search)
{
    double um = 0.5 * (u0 + u1);
    double vm = 0.5 * (v0 + v1);
    double a00, b00, g00;
    double a10, b10, g10;
    double a01, b01, g01;
    double a11, b11, g11;
    double amm, bmm, gmm;
    int c00, c10, c01, c11, cc;
    double err;
    int need_split;

    /* evaluate the four corners and the centre */
    agrid_make_abg(u0, v0, search, &a00, &b00, &g00);
    agrid_make_abg(u1, v0, search, &a10, &b10, &g10);
    agrid_make_abg(u0, v1, search, &a01, &b01, &g01);
    agrid_make_abg(u1, v1, search, &a11, &b11, &g11);
    agrid_make_abg(um, vm, search, &amm, &bmm, &gmm);

    c00 = agrid_add_or_get(a00, b00, g00);
    c10 = agrid_add_or_get(a10, b10, g10);
    c01 = agrid_add_or_get(a01, b01, g01);
    c11 = agrid_add_or_get(a11, b11, g11);
    cc  = agrid_add_or_get(amm, bmm, gmm);

    err = agrid_raw_interp_error(c00, c10, c01, c11, cc);
    need_split = (depth < AGRID_MIN_DEPTH) ||
                 (err > AGRID_TOL && depth < AGRID_MAX_DEPTH);

    if (need_split) {
        agrid_subdivide(u0, um, v0, vm, depth + 1, search);
        agrid_subdivide(um, u1, v0, vm, depth + 1, search);
        agrid_subdivide(u0, um, vm, v1, depth + 1, search);
        agrid_subdivide(um, u1, vm, v1, depth + 1, search);
    }
    /* leaf cell: corners and centre are already in the cache */
}

@*1 Stale-context helpers.

@<Definition for |agrid_save_context|@>=
static void agrid_save_context(struct measure_type m, struct invert_type r)
{
    AGrid_slab_n         = m.slab_index;
    AGrid_slab_n_top     = m.slab_top_slide_index;
    AGrid_slab_n_bottom  = m.slab_bottom_slide_index;
    AGrid_slab_b_top     = r.slab.b_top_slide;
    AGrid_slab_b_bottom  = r.slab.b_bottom_slide;
    AGrid_slab_cos_angle = m.slab_cos_angle;
    AGrid_num_spheres    = m.num_spheres;
    AGrid_num_measures   = m.num_measures;
    AGrid_m_u            = m.m_u;
    /* Lost-light values are NOT saved: subdivision uses include_lost=False,
       so the grid is valid across MC iterations for the same wavelength. */

    AGrid_fixed_g = r.slab.g;
    AGrid_fixed_b = r.slab.b;
    AGrid_fixed_a = (r.default_a == UNINITIALIZED) ? 0.0 : r.default_a;

    AGrid_Search      = r.search;
    AGrid_Initialized = 1;
}

@*1 Public interface.

|AGrid_Valid| returns non-zero if the current cache was built for the
same physical state as $(m, r)$.  The check mirrors |Valid_Grid| in
\.{iad\_calc.c} but is independent of the dense-grid allocation.

@<Prototype for |AGrid_Valid|@>=
int AGrid_Valid(struct measure_type m, struct invert_type r)

@ @<Definition for |AGrid_Valid|@>=
    @<Prototype for |AGrid_Valid|@>
{
    if (!AGrid_Initialized || AGrid_N == 0) return 0;
    if (AGrid_Search != r.search) return 0;

    if (m.slab_index             != AGrid_slab_n)        return 0;
    if (m.slab_cos_angle         != AGrid_slab_cos_angle) return 0;
    if (m.slab_top_slide_index   != AGrid_slab_n_top)    return 0;
    if (m.slab_bottom_slide_index!= AGrid_slab_n_bottom) return 0;
    if (r.slab.b_top_slide       != AGrid_slab_b_top)    return 0;
    if (r.slab.b_bottom_slide    != AGrid_slab_b_bottom) return 0;

    if (m.num_spheres            != AGrid_num_spheres)   return 0;
    if (m.num_measures           != AGrid_num_measures)  return 0;
    if (m.num_measures == 3 && m.m_u != AGrid_m_u)      return 0;

    /* Subdivision uses include_lost=False so lost-light terms do NOT
       invalidate the grid (matches Python AGrid stale-check behaviour).
       The distance at query time recomputes with current lost-light values. */

    /* fixed-axis value for the current search */
    if (r.search == FIND_AB && r.slab.g != AGrid_fixed_g) return 0;
    if (r.search == FIND_AG && r.slab.b != AGrid_fixed_b) return 0;
    if (r.search == FIND_BG) {
        double fa = (r.default_a == UNINITIALIZED) ? 0.0 : r.default_a;
        if (fa != AGrid_fixed_a) return 0;
    }

    return 1;
}

@ |AGrid_Build| resets the cache and fills it via adaptive quadtree
subdivision.  Call |AGrid_Valid| first to avoid unnecessary rebuilds.

@<Prototype for |AGrid_Build|@>=
void AGrid_Build(struct measure_type m, struct invert_type r)

@ @<Definition for |AGrid_Build|@>=
    @<Prototype for |AGrid_Build|@>
{
    double u0, u1, v0, v1;
    int search = r.search;

    /* reset cache */
    AGrid_N = 0;

    /* save fixed-axis values before setting calc state */
    if (search == FIND_AB) AGrid_fixed_g = r.slab.g;
    if (search == FIND_AG) AGrid_fixed_b = r.slab.b;
    if (search == FIND_BG)
        AGrid_fixed_a = (r.default_a == UNINITIALIZED) ? 0.0 : r.default_a;

    /* set global MM/RR state so abg_eval uses the right sphere parameters */
    Set_Calc_State(m, r);

    if (Debug(DEBUG_GRID)) {
        switch (search) {
        case FIND_AB:
            fprintf(stderr, "AGRID: Filling AB grid (g=%.5f)\n", AGrid_fixed_g);
            break;
        case FIND_AG:
            fprintf(stderr, "AGRID: Filling AG grid (b=%.5f)\n", AGrid_fixed_b);
            break;
        case FIND_BG:
            fprintf(stderr, "AGRID: Filling BG grid (a=%.5f)\n", AGrid_fixed_a);
            break;
        default:
            fprintf(stderr, "AGRID: Filling grid for search=%d\n", search);
            break;
        }
    }

    /* determine coordinate ranges */
    switch (search) {
    case FIND_AB:
        u0 = 0.0; u1 = 1.0;
        v0 = AGRID_MIN_LOG_B; v1 = AGRID_MAX_LOG_B_AB;
        break;
    case FIND_AG:
        u0 = 0.0; u1 = 1.0;
        v0 = 0.0; v1 = 1.0;
        break;
    case FIND_BG:
        u0 = AGRID_MIN_LOG_B; u1 = AGRID_MAX_LOG_B_BG;
        v0 = 0.0; v1 = 1.0;
        break;
    default:
        return;
    }

    agrid_subdivide(u0, u1, v0, v1, 0, search);

    agrid_save_context(m, r);

    if (Debug(DEBUG_GRID))
        fprintf(stderr, "AGRID: Built %d entries for search=%d\n", AGrid_N, search);
    if (Debug(DEBUG_LOST_LIGHT))
        fprintf(stderr, "GRID: %d AD evaluations to fill the grid\n", AGrid_N);
}

@ |AGrid_Fill_Guesses| finds the best candidates from the adaptive cache
and fills at most |max_n| entries in the caller's |guesses| array.  The
entries are sorted by sphere-corrected distance (smallest first).  Returns
the number of entries actually filled.

The function also performs a local Cartesian closure: the top-|AGRID_CLOSURE_N|
candidates contribute their axis values to a small cross-product set, and all
cross-combinations are evaluated.  This catches the case where the true minimum
lies between diagonal samples.

@<Prototype for |AGrid_Fill_Guesses|@>=
int AGrid_Fill_Guesses(double m_r, double m_t,
                       guess_type *guesses, int max_n)

@ @<Definition for |AGrid_Fill_Guesses|@>=
    @<Prototype for |AGrid_Fill_Guesses|@>
{
    int i, j, n_out;
    double *all_dist;
    double best_dist[AGRID_CLOSURE_N];
    int    best_idx [AGRID_CLOSURE_N];
    int    n_best;

    if (AGrid_N == 0) return 0;

    /* -- step 1: compute distances for all entries; keep top AGRID_CLOSURE_N -- */
    all_dist = (double *) malloc((size_t) AGrid_N * sizeof(double));
    if (all_dist == NULL) {
        fprintf(stderr, "AGrid: out of memory in AGrid_Fill_Guesses\n");
        exit(EXIT_FAILURE);
    }

    n_best = 0;
    for (i = 0; i < AGrid_N; i++) {
        double d = abg_stored_distance(AGrid_Cache[i].a,
                                       AGrid_Cache[i].b,
                                       AGrid_Cache[i].g,
                                       AGrid_Cache[i].ur1,
                                       AGrid_Cache[i].ut1,
                                       AGrid_Cache[i].uru,
                                       AGrid_Cache[i].utu);
        all_dist[i] = d;
        if (Debug(DEBUG_GRID) && d < 0.15) {
            fprintf(stderr, "AGRID_ALL: a=%8.5f b=%9.5f g=%7.4f ur1=%8.5f ut1=%8.5f d=%9.6f\n",
                    AGrid_Cache[i].a, AGrid_Cache[i].b, AGrid_Cache[i].g,
                    AGrid_Cache[i].ur1, AGrid_Cache[i].ut1, d);
        }
        if (n_best < AGRID_CLOSURE_N) {
            best_dist[n_best] = d;
            best_idx [n_best] = i;
            n_best++;
        } else {
            /* replace worst */
            int worst = 0;
            for (j = 1; j < n_best; j++)
                if (best_dist[j] > best_dist[worst]) worst = j;
            if (d < best_dist[worst]) {
                best_dist[worst] = d;
                best_idx [worst] = i;
            }
        }
    }

    /* -- step 2: local Cartesian closure over unique axis values --
       Scan ALL cache entries in order of increasing distance (matching
       Python AGrid._candidate_axis_values) and collect up to
       AGRID_CLOSURE_N distinct values per axis independently.  Python
       scans the full sorted ranked list for each axis separately, which
       finds good axis values even when they appear in low-ranked entries.
       Then enrich with midpoints between consecutive sorted axis values
       (geometric for b, arithmetic for a and g), matching Python's
       _enriched_axis_values. */
    {
        double axis0[AGRID_CLOSURE_N * 2];
        double axis1[AGRID_CLOSURE_N * 2];
        int n0 = 0, n1 = 0;
        int k, l, found;
        int use_geo0, use_geo1; /* use geometric midpoint for this axis? */

        /* axis0 uses geometric mean only for the b axis in find_bg */
        use_geo0 = (AGrid_Search == FIND_BG);
        /* axis1 uses geometric mean only for the b axis in find_ab */
        use_geo1 = (AGrid_Search == FIND_AB);

        /* collect top-AGRID_CLOSURE_N distinct axis0 values in order of
           increasing distance across the full cache */
        while (n0 < AGRID_CLOSURE_N) {
            double best_d = 1e30;
            int    best_i = -1;
            double vc;
            for (i = 0; i < AGrid_N; i++) {
                switch (AGrid_Search) {
                case FIND_AB: vc = AGrid_Cache[i].a; break;
                case FIND_AG: vc = AGrid_Cache[i].a; break;
                case FIND_BG: vc = AGrid_Cache[i].b; break;
                default:      vc = AGrid_Cache[i].a; break;
                }
                found = 0;
                for (l = 0; l < n0; l++) if (axis0[l] == vc) { found = 1; break; }
                if (!found && all_dist[i] < best_d) {
                    best_d = all_dist[i];
                    best_i = i;
                }
            }
            if (best_i < 0) break; /* no more distinct values */
            switch (AGrid_Search) {
            case FIND_AB: axis0[n0++] = AGrid_Cache[best_i].a; break;
            case FIND_AG: axis0[n0++] = AGrid_Cache[best_i].a; break;
            case FIND_BG: axis0[n0++] = AGrid_Cache[best_i].b; break;
            default:      axis0[n0++] = AGrid_Cache[best_i].a; break;
            }
        }

        /* collect top-AGRID_CLOSURE_N distinct axis1 values similarly */
        while (n1 < AGRID_CLOSURE_N) {
            double best_d = 1e30;
            int    best_i = -1;
            double vc;
            for (i = 0; i < AGrid_N; i++) {
                switch (AGrid_Search) {
                case FIND_AB: vc = AGrid_Cache[i].b; break;
                case FIND_AG: vc = AGrid_Cache[i].g; break;
                case FIND_BG: vc = AGrid_Cache[i].g; break;
                default:      vc = AGrid_Cache[i].b; break;
                }
                found = 0;
                for (l = 0; l < n1; l++) if (axis1[l] == vc) { found = 1; break; }
                if (!found && all_dist[i] < best_d) {
                    best_d = all_dist[i];
                    best_i = i;
                }
            }
            if (best_i < 0) break;
            switch (AGrid_Search) {
            case FIND_AB: axis1[n1++] = AGrid_Cache[best_i].b; break;
            case FIND_AG: axis1[n1++] = AGrid_Cache[best_i].g; break;
            case FIND_BG: axis1[n1++] = AGrid_Cache[best_i].g; break;
            default:      axis1[n1++] = AGrid_Cache[best_i].b; break;
            }
        }

        /* enrich axis0 with midpoints between consecutive ORIGINAL sorted
           values.  Mirrors Python _enriched_axis_values: sort the originals,
           then append (not insert) midpoints so further iterations don't
           cascade midpoints-of-midpoints. */
        if (n0 >= 2) {
            double sorted0[AGRID_CLOSURE_N];
            int n_sorted = n0;
            for (k = 0; k < n0; k++) sorted0[k] = axis0[k];
            for (k = 1; k < n_sorted; k++) {
                double tmp = sorted0[k];
                l = k - 1;
                while (l >= 0 && sorted0[l] > tmp) { sorted0[l+1] = sorted0[l]; l--; }
                sorted0[l+1] = tmp;
            }
            for (k = 0; k < n_sorted - 1 && n0 < AGRID_CLOSURE_N * 2; k++) {
                double lo = sorted0[k], hi = sorted0[k+1];
                double mid = (use_geo0 && lo > 0 && hi > 0) ? sqrt(lo * hi) : 0.5*(lo+hi);
                found = 0;
                for (l = 0; l < n0; l++) if (axis0[l] == mid) { found = 1; break; }
                if (!found) axis0[n0++] = mid;
            }
        }

        /* enrich axis1 similarly */
        if (n1 >= 2) {
            double sorted1[AGRID_CLOSURE_N];
            int n_sorted = n1;
            for (k = 0; k < n1; k++) sorted1[k] = axis1[k];
            for (k = 1; k < n_sorted; k++) {
                double tmp = sorted1[k];
                l = k - 1;
                while (l >= 0 && sorted1[l] > tmp) { sorted1[l+1] = sorted1[l]; l--; }
                sorted1[l+1] = tmp;
            }
            for (k = 0; k < n_sorted - 1 && n1 < AGRID_CLOSURE_N * 2; k++) {
                double lo = sorted1[k], hi = sorted1[k+1];
                double mid = (use_geo1 && lo > 0 && hi > 0) ? sqrt(lo * hi) : 0.5*(lo+hi);
                found = 0;
                for (l = 0; l < n1; l++) if (axis1[l] == mid) { found = 1; break; }
                if (!found) axis1[n1++] = mid;
            }
        }

        /* evaluate all cross-combinations and update top candidates */
        for (k = 0; k < n0; k++) {
            for (l = 0; l < n1; l++) {
                double a, b, g;
                double d;
                int idx;

                switch (AGrid_Search) {
                case FIND_AB: a = axis0[k]; b = axis1[l]; g = AGrid_fixed_g; break;
                case FIND_AG: a = axis0[k]; b = AGrid_fixed_b; g = axis1[l]; break;
                case FIND_BG: a = AGrid_fixed_a; b = axis0[k]; g = axis1[l]; break;
                default:      a = axis0[k]; b = axis1[l]; g = AGrid_fixed_g; break;
                }

                if (b < 1e-8) b = 1e-8;
                if (g < -MAX_ABS_G) g = -MAX_ABS_G;
                if (g >  MAX_ABS_G) g =  MAX_ABS_G;

                idx = agrid_add_or_get(a, b, g);
                d   = abg_stored_distance(AGrid_Cache[idx].a,
                                          AGrid_Cache[idx].b,
                                          AGrid_Cache[idx].g,
                                          AGrid_Cache[idx].ur1,
                                          AGrid_Cache[idx].ut1,
                                          AGrid_Cache[idx].uru,
                                          AGrid_Cache[idx].utu);

                /* update top-AGRID_CLOSURE_N if better */
                if (n_best < AGRID_CLOSURE_N) {
                    best_dist[n_best] = d;
                    best_idx [n_best] = idx;
                    n_best++;
                } else {
                    int worst = 0;
                    for (j = 1; j < n_best; j++)
                        if (best_dist[j] > best_dist[worst]) worst = j;
                    if (d < best_dist[worst]) {
                        best_dist[worst] = d;
                        best_idx [worst] = idx;
                    }
                }
            }
        }
    }

    free(all_dist);

    /* -- step 3: sort best_idx by distance (insertion sort, tiny N) -- */
    for (i = 1; i < n_best; i++) {
        int    ti = best_idx[i];
        double td = best_dist[i];
        j = i - 1;
        while (j >= 0 && best_dist[j] > td) {
            best_dist[j+1] = best_dist[j];
            best_idx [j+1] = best_idx[j];
            j--;
        }
        best_dist[j+1] = td;
        best_idx [j+1] = ti;
    }

    /* -- step 4: copy into caller's guesses array -- */
    n_out = (n_best < max_n) ? n_best : max_n;
    for (i = 0; i < n_out; i++) {
        agrid_entry_t *e = &AGrid_Cache[best_idx[i]];
        guesses[i].a        = e->a;
        guesses[i].b        = e->b;
        guesses[i].g        = e->g;
        guesses[i].distance = best_dist[i];
        guesses[i].ur1_lost = 0.0;
        guesses[i].ut1_lost = 0.0;
        guesses[i].uru_lost = 0.0;
        guesses[i].utu_lost = 0.0;
    }

    if (Debug(DEBUG_BEST_GUESS)) {
        fprintf(stderr, "BEST: AGRID GUESSES\n");
        fprintf(stderr, "BEST:  k      albedo          b          g   distance\n");
        for (i = 0; i < n_out && i < 7; i++) {
            fprintf(stderr, "BEST:%3d  ", i);
            fprintf(stderr, "%10.5f ", guesses[i].a);
            fprintf(stderr, "%10.5f ", guesses[i].b);
            fprintf(stderr, "%10.5f ", guesses[i].g);
            fprintf(stderr, "%10.5f\n", guesses[i].distance);
        }
    }

    return n_out;
}
