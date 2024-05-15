/* Autogenerated v3-16-3 from https://github.com/scottprahl/iad */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "nr_util.h"
#include "nr_zbrent.h"
#include "ad_globl.h"
#include "ad_frsnl.h"
#include "ad_prime.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"

#define ABIT 1e-6
#define A_COLUMN 1
#define B_COLUMN 2
#define G_COLUMN 3
#define URU_COLUMN 4
#define UTU_COLUMN 5
#define UR1_COLUMN 6
#define UT1_COLUMN 7
#define REFLECTION_SPHERE 1
#define TRANSMISSION_SPHERE 0
#define T_TRUST_FACTOR 1

#define SWAP(a,b) {double swap= (a);(a)= (b);(b)= swap;}

static int CALCULATING_GRID = 0;
static struct measure_type MM;
static struct invert_type RR;
static struct measure_type MGRID;
static struct invert_type RGRID;
static double **The_Grid = NULL;
static double GG_a;
static double GG_b;
static double GG_g;
static double GG_bs;
static double GG_ba;
static boolean_type The_Grid_Initialized = FALSE;
static boolean_type The_Grid_Search = -1;

void Set_Calc_State(struct measure_type m, struct invert_type r)
{
    memcpy(&MM, &m, sizeof(struct measure_type));
    memcpy(&RR, &r, sizeof(struct invert_type));
}

void Get_Calc_State(struct measure_type *m, struct invert_type *r)
{
    memcpy(m, &MM, sizeof(struct measure_type));
    memcpy(r, &RR, sizeof(struct invert_type));
}

boolean_type Same_Calc_State(struct measure_type m, struct invert_type r)
{
    if (The_Grid == NULL)
        return FALSE;
    if (!The_Grid_Initialized)
        return FALSE;

    if (r.search != RR.search)
        return FALSE;
    if (r.method.quad_pts != RR.method.quad_pts)
        return FALSE;
    if (r.slab.a != RR.slab.a)
        return FALSE;
    if (r.slab.b != RR.slab.b)
        return FALSE;
    if (r.slab.g != RR.slab.g)
        return FALSE;
    if (r.slab.phase_function != RR.slab.phase_function)
        return FALSE;
    if (r.slab.n_slab != RR.slab.n_slab)
        return FALSE;
    if (r.slab.n_top_slide != RR.slab.n_top_slide)
        return FALSE;
    if (r.slab.n_bottom_slide != RR.slab.n_bottom_slide)
        return FALSE;
    if (r.slab.b_top_slide != RR.slab.b_top_slide)
        return FALSE;
    if (r.slab.b_bottom_slide != RR.slab.b_bottom_slide)
        return FALSE;
    if (r.slab.cos_angle != RR.slab.cos_angle)
        return FALSE;
    if ((m.num_measures == 3) && (m.m_u != MGRID.m_u))
        return (FALSE);
    return TRUE;
}

void Fill_AB_Grid(struct measure_type m, struct invert_type r);

void Fill_AG_Grid(struct measure_type m, struct invert_type r);

void RT_Flip(int flip, int n, struct AD_slab_type *slab, double *UR1, double *UT1, double *URU, double *UTU)
{
    double correct_UR1, correct_URU;
    RT(n, slab, UR1, UT1, URU, UTU);
    if (flip) {
        correct_UR1 = *UR1;
        correct_URU = *URU;
        SWAP(slab->n_top_slide, slab->n_bottom_slide)
            SWAP(slab->b_top_slide, slab->b_bottom_slide)
            RT(n, slab, UR1, UT1, URU, UTU);
        SWAP(slab->n_top_slide, slab->n_bottom_slide)
            SWAP(slab->b_top_slide, slab->b_bottom_slide)
            * UR1 = correct_UR1;
        *URU = correct_URU;
    }
}

void Allocate_Grid(search_type s)
{
    (void) s;
    The_Grid = dmatrix(0, GRID_SIZE * GRID_SIZE, 1, 7);
    if (The_Grid == NULL)
        AD_error("unable to allocate the grid matrix");
    The_Grid_Initialized = FALSE;
}

boolean_type Valid_Grid(struct measure_type m, struct invert_type r)
{
    int s = r.search;

    if (The_Grid == NULL) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because NULL\n");
        return (FALSE);
    }
    if (!The_Grid_Initialized) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because not initialized\n");
        return (FALSE);
    }

    if (The_Grid_Search != s) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because search type changed\n");
        return (FALSE);
    }

    if ((m.num_measures == 3) && (m.m_u != MGRID.m_u)) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because unscattered light changed\n");
        return (FALSE);
    }

    if (m.slab_index != MGRID.slab_index) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because slab refractive index changed\n");
        return (FALSE);
    }
    if (m.slab_cos_angle != MGRID.slab_cos_angle) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because light angle changed\n");
        return (FALSE);
    }

    if (m.slab_top_slide_index != MGRID.slab_top_slide_index) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because top slide index changed\n");
        return (FALSE);
    }

    if (m.slab_bottom_slide_index != MGRID.slab_bottom_slide_index) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because bottom slide index changed\n");
        return (FALSE);
    }

    if (s == FIND_AB && r.slab.g != RGRID.slab.g) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because anisotropy changed\n");
        return (FALSE);
    }

    if (s == FIND_AG && r.slab.b != RGRID.slab.b) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because optical depth changed\n");
        return (FALSE);
    }

    if (s == FIND_BsG && r.default_ba != RGRID.default_ba) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because mu_a changed\n");
        return (FALSE);
    }

    if (s == FIND_BaG && r.default_bs != RGRID.default_bs) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr, "GRID: Fill because mu_s changed\n");
        return (FALSE);
    }

    return (TRUE);
}

static void fill_grid_entry(int i, int j)
{
    double ur1, ut1, uru, utu;

    if (RR.slab.b <= 1e-6)
        RR.slab.b = 1e-6;

    if (Debug(DEBUG_GRID_CALC) && i == 0 && j == 0) {
        fprintf(stderr, "+   i   j ");
        fprintf(stderr, "      a         b          g     |");
        fprintf(stderr, "     M_R        grid  |");
        fprintf(stderr, "     M_T        grid\n");
    }

    if (Debug(DEBUG_EVERY_CALC)) {
        if (!CALCULATING_GRID)
            fprintf(stderr, "a=%8.5f b=%10.5f g=%8.5f ", RR.slab.a, RR.slab.b, RR.slab.g);
        else {
            if (j == 0)
                fprintf(stderr, ".");
            if (i + 1 == GRID_SIZE && j == 0)
                fprintf(stderr, "\n");
        }
    }

    RT_Flip(MM.flip_sample, RR.method.quad_pts, &RR.slab, &ur1, &ut1, &uru, &utu);

    if (Debug(DEBUG_EVERY_CALC) && !CALCULATING_GRID)
        fprintf(stderr, "ur1=%8.5f ut1=%8.5f\n", ur1, ut1);

    The_Grid[GRID_SIZE * i + j][A_COLUMN] = RR.slab.a;
    The_Grid[GRID_SIZE * i + j][B_COLUMN] = RR.slab.b;
    The_Grid[GRID_SIZE * i + j][G_COLUMN] = RR.slab.g;
    The_Grid[GRID_SIZE * i + j][UR1_COLUMN] = ur1;
    The_Grid[GRID_SIZE * i + j][UT1_COLUMN] = ut1;
    The_Grid[GRID_SIZE * i + j][URU_COLUMN] = uru;
    The_Grid[GRID_SIZE * i + j][UTU_COLUMN] = utu;

    if (Debug(DEBUG_GRID_CALC)) {
        fprintf(stderr, "+ %3d %3d ", i, j);
        fprintf(stderr, "%10.5f %10.5f %10.5f |", RR.slab.a, RR.slab.b, RR.slab.g);
        fprintf(stderr, "%10.5f %10.5f |", MM.m_r, uru);
        fprintf(stderr, "%10.5f %10.5f \n", MM.m_t, utu);
    }
}

void Fill_Grid(struct measure_type m, struct invert_type r, int force_new)
{
    if (force_new || !Same_Calc_State(m, r)) {
        switch (r.search) {
        case FIND_AB:
            Fill_AB_Grid(m, r);
            break;
        case FIND_AG:
            Fill_AG_Grid(m, r);
            break;
        case FIND_BG:
            Fill_BG_Grid(m, r);
            break;
        case FIND_BaG:
            Fill_BaG_Grid(m, r);
            break;
        case FIND_BsG:
            Fill_BsG_Grid(m, r);
            break;
        default:
            AD_error("Attempt to fill grid for unknown search case.");
        }
    }

    Get_Calc_State(&MGRID, &RGRID);
}

void Near_Grid_Points(double r, double t, search_type s, int *i_min, int *j_min)
{
    int i, j;
    double fval;
    double smallest = 10.0;
    struct measure_type old_mm;
    struct invert_type old_rr;
    (void) r;
    (void) t;
    (void) s;

    if (Debug(DEBUG_GRID))
        fprintf(stderr, "GRID: Finding best grid points\n");
    Get_Calc_State(&old_mm, &old_rr);

    *i_min = 0;
    *j_min = 0;
    for (i = 0; i < GRID_SIZE; i++) {
        for (j = 0; j < GRID_SIZE; j++) {

            CALCULATING_GRID = 1;
            fval = Calculate_Grid_Distance(i, j);
            CALCULATING_GRID = 0;

            if (fval < smallest) {
                *i_min = i;
                *j_min = j;
                smallest = fval;
            }
        }
    }

    Set_Calc_State(old_mm, old_rr);
}

void Fill_AB_Grid(struct measure_type m, struct invert_type r)
{
    int i, j;
    double min_log_b = -8;
    double max_log_b = +8;

    if (Debug(DEBUG_GRID))
        fprintf(stderr, "GRID: Filling AB grid (g=%.5f)\n", RR.slab.g);

    if (The_Grid == NULL)
        Allocate_Grid(r.search);

    GG_a = 0.0;
    GG_b = 0.0;
    GG_g = 0.0;
    GG_bs = 0.0;
    GG_ba = 0.0;

    Set_Calc_State(m, r);

    GG_g = RR.slab.g;
    for (i = 0; i < GRID_SIZE; i++) {
        double x = (double) i / (GRID_SIZE - 1.0);
        RR.slab.b = exp(min_log_b + (max_log_b - min_log_b) * x);
        for (j = 0; j < GRID_SIZE; j++) {

            {
                double x = (double) j / (GRID_SIZE - 1.0);
                RR.slab.a = 1.0 - (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
            }

            fill_grid_entry(i, j);
        }
    }

    The_Grid_Initialized = TRUE;
    The_Grid_Search = FIND_AB;
}

void Fill_AG_Grid(struct measure_type m, struct invert_type r)
{
    int i, j;
    double max_a = -10;
    double min_a = 10;

    if (Debug(DEBUG_GRID))
        fprintf(stderr, "GRID: Filling AG grid\n");

    if (The_Grid == NULL)
        Allocate_Grid(r.search);

    GG_a = 0.0;
    GG_b = 0.0;
    GG_g = 0.0;
    GG_bs = 0.0;
    GG_ba = 0.0;

    Set_Calc_State(m, r);
    GG_b = r.slab.b;
    for (i = 0; i < GRID_SIZE; i++) {

        {
            double x = (double) i / (GRID_SIZE - 1.0);
            double xx = (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
            RR.slab.g = (1.0 - 2.0 * xx) * MAX_ABS_G;
        }

        for (j = 0; j < GRID_SIZE; j++) {

            {
                double x = (double) j / (GRID_SIZE - 1.0);
                RR.slab.a = 1.0 - (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
            }

            fill_grid_entry(i, j);
            if (RR.slab.a < min_a)
                min_a = RR.slab.a;
            if (RR.slab.a > max_a)
                max_a = RR.slab.a;
        }
    }

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: a        = %9.7f to %9.7f \n", min_a, max_a);
        fprintf(stderr, "GRID: b        = %9.5f \n", r.slab.b);
        fprintf(stderr, "GRID: g  range = %9.6f to %9.6f \n", -MAX_ABS_G, MAX_ABS_G);
    }

    The_Grid_Initialized = TRUE;
    The_Grid_Search = FIND_AG;
}

void Fill_BG_Grid(struct measure_type m, struct invert_type r)
{
    int i, j;
    double min_log_b = -8;
    double max_log_b = +10;

    if (The_Grid == NULL)
        Allocate_Grid(r.search);

    GG_a = 0.0;
    GG_b = 0.0;
    GG_g = 0.0;
    GG_bs = 0.0;
    GG_ba = 0.0;

    if (Debug(DEBUG_GRID))
        fprintf(stderr, "GRID: Filling BG grid\n");

    Set_Calc_State(m, r);
    RR.slab.a = RR.default_a;
    GG_a = RR.slab.a;

    for (i = 0; i < GRID_SIZE; i++) {
        double x = (double) i / (GRID_SIZE - 1.0);
        RR.slab.b = exp(min_log_b + (max_log_b - min_log_b) * x);
        for (j = 0; j < GRID_SIZE; j++) {

            {
                double x = (double) j / (GRID_SIZE - 1.0);
                double xx = (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
                RR.slab.g = (1.0 - 2.0 * xx) * MAX_ABS_G;
            }

            fill_grid_entry(i, j);
        }
    }

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: a        = %9.7f \n", RR.default_a);
        fprintf(stderr, "GRID: b  range = %9.5f to %9.3f \n", exp(min_log_b), exp(max_log_b));
        fprintf(stderr, "GRID: g  range = %9.6f to %9.6f \n", -MAX_ABS_G, MAX_ABS_G);
    }

    The_Grid_Initialized = TRUE;
    The_Grid_Search = FIND_BG;
}

void Fill_BaG_Grid(struct measure_type m, struct invert_type r)
{
    int i, j;
    double max_a = -10;
    double min_a = 10;
    double bs = r.default_bs;
    double min_log_ba = -8;
    double max_log_ba = +10;

    if (The_Grid == NULL)
        Allocate_Grid(r.search);

    GG_a = 0.0;
    GG_b = 0.0;
    GG_g = 0.0;
    GG_bs = 0.0;
    GG_ba = 0.0;

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: Filling BaG grid\n");
        fprintf(stderr, "GRID:       bs = %9.5f\n", bs);
        fprintf(stderr, "GRID: ba range = %9.6f to %9.3f \n", exp(min_log_ba), exp(max_log_ba));
    }

    Set_Calc_State(m, r);
    GG_bs = bs;
    for (i = 0; i < GRID_SIZE; i++) {
        double x = (double) i / (GRID_SIZE - 1.0);
        double ba = exp(min_log_ba + (max_log_ba - min_log_ba) * x);
        RR.slab.b = ba + bs;
        if (RR.slab.b > 0)
            RR.slab.a = bs / RR.slab.b;
        else
            RR.slab.a = 0;
        if (RR.slab.a < 0)
            RR.slab.a = 0;
        if (RR.slab.a < min_a)
            min_a = RR.slab.a;
        if (RR.slab.a > max_a)
            max_a = RR.slab.a;

        for (j = 0; j < GRID_SIZE; j++) {

            {
                double x = (double) j / (GRID_SIZE - 1.0);
                double xx = (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
                RR.slab.g = (1.0 - 2.0 * xx) * MAX_ABS_G;
            }

            fill_grid_entry(i, j);
        }
    }

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: a        = %9.7f to %9.7f \n", min_a, max_a);
        fprintf(stderr, "GRID: b  range = %9.5f to %9.3f \n", exp(min_log_ba) + bs, exp(max_log_ba) + bs);
        fprintf(stderr, "GRID: g  range = %9.6f to %9.6f \n", -MAX_ABS_G, MAX_ABS_G);
    }

    The_Grid_Initialized = TRUE;
    The_Grid_Search = FIND_BaG;
}

void Fill_BsG_Grid(struct measure_type m, struct invert_type r)
{
    int i, j;
    double max_a = -10;
    double min_a = 10;
    double ba = r.default_ba;
    double min_log_bs = -8;
    double max_log_bs = +10;

    if (The_Grid == NULL)
        Allocate_Grid(r.search);

    GG_a = 0.0;
    GG_b = 0.0;
    GG_g = 0.0;
    GG_bs = 0.0;
    GG_ba = 0.0;

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: Filling BsG grid\n");
        fprintf(stderr, "GRID:       ba = %9.5f\n", ba);
        fprintf(stderr, "GRID: bs range = %9.6f to %9.3f \n", exp(min_log_bs), exp(max_log_bs));
    }

    Set_Calc_State(m, r);
    GG_ba = RR.default_ba;
    for (i = 0; i < GRID_SIZE; i++) {
        double x = (double) i / (GRID_SIZE - 1.0);
        double bs = exp(min_log_bs + (max_log_bs - min_log_bs) * x);
        RR.slab.b = ba + bs;
        if (RR.slab.b > 0)
            RR.slab.a = 1 - RR.default_ba / RR.slab.b;
        else
            RR.slab.a = 0;
        if (RR.slab.a < 0)
            RR.slab.a = 0;
        if (RR.slab.a < min_a)
            min_a = RR.slab.a;
        if (RR.slab.a > max_a)
            max_a = RR.slab.a;

        for (j = 0; j < GRID_SIZE; j++) {

            {
                double x = (double) j / (GRID_SIZE - 1.0);
                double xx = (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
                RR.slab.g = (1.0 - 2.0 * xx) * MAX_ABS_G;
            }

            fill_grid_entry(i, j);
        }
    }

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: a  range = %9.7f to %9.7f \n", min_a, max_a);
        fprintf(stderr, "GRID: b  range = %9.5f to %9.3f \n", exp(min_log_bs) + ba, exp(max_log_bs) + ba);
        fprintf(stderr, "GRID: g  range = %9.6f to %9.6f \n", -MAX_ABS_G, MAX_ABS_G);
    }

    The_Grid_Initialized = TRUE;
    The_Grid_Search = FIND_BsG;
}

void Grid_ABG(int i, int j, guess_type *guess)
{
    if (0 <= i && i < GRID_SIZE && 0 <= j && j < GRID_SIZE) {
        guess->a = The_Grid[GRID_SIZE * i + j][A_COLUMN];
        guess->b = The_Grid[GRID_SIZE * i + j][B_COLUMN];
        guess->g = The_Grid[GRID_SIZE * i + j][G_COLUMN];
        guess->distance = Calculate_Grid_Distance(i, j);
    }
    else {
        guess->a = 0.5;
        guess->b = 0.5;
        guess->g = 0.5;
        guess->distance = 999;
    }
}

double Gain(int sphere, struct measure_type m, double uru_sample, double uru_third)
{
    double inv_gain;

    if (sphere == REFLECTION_SPHERE) {
        if (m.baffle_r) {
            inv_gain = m.rw_r + (m.at_r / m.aw_r) * uru_third;
            inv_gain *= m.aw_r + (1 - m.at_r) * (m.ad_r * m.rd_r + m.as_r * uru_sample);
            inv_gain = 1.0 - inv_gain;
        }
        else {
            inv_gain = 1.0 - m.aw_r * m.rw_r - m.ad_r * m.rd_r - m.as_r * uru_sample - m.at_r * uru_third;
        }
    }
    else if (m.baffle_t) {
        inv_gain = m.rw_t + (m.at_t / m.aw_t) * uru_third;
        inv_gain *= m.aw_t + (1 - m.at_t) * (m.ad_t * m.rd_t + m.as_t * uru_sample);
        inv_gain = 1.0 - inv_gain;
    }
    else {
        inv_gain = 1.0 - m.aw_t * m.rw_t - m.ad_t * m.rd_t - m.as_t * uru_sample - m.at_t * uru_third;
    }
    return 1.0 / inv_gain;
}

double Gain_11(struct measure_type m, double URU, double tdiffuse)
{
    double G, GP, G11;

    G = Gain(REFLECTION_SPHERE, m, URU, 0);
    GP = Gain(TRANSMISSION_SPHERE, m, URU, 0);

    G11 = G / (1 - m.as_r * m.as_t * m.aw_r * m.aw_t * (1 - m.at_r) * (1 - m.at_t)
        * G * GP * tdiffuse * tdiffuse);

    return G11;
}

double Gain_22(struct measure_type m, double URU, double tdiffuse)
{
    double G, GP, G22;

    G = Gain(REFLECTION_SPHERE, m, URU, 0);
    GP = Gain(TRANSMISSION_SPHERE, m, URU, 0);

    G22 = GP / (1 - m.as_r * m.as_t * m.aw_r * m.aw_t * (1 - m.at_r) * (1 - m.at_t)
        * G * GP * tdiffuse * tdiffuse);

    return G22;
}

double Two_Sphere_R(struct measure_type m, double UR1, double URU, double UT1, double UTU)
{
    double x, GP;
    GP = Gain(TRANSMISSION_SPHERE, m, URU, 0);

    x = m.ad_r * (1 - m.at_r) * m.rw_r * Gain_11(m, URU, UTU);
    x *= (1 - m.f_r) * UR1 + m.rw_r * m.f_r + (1 - m.f_r) * m.as_t * (1 - m.at_t) * m.rw_t * UT1 * UTU * GP;
    return x;
}

double Two_Sphere_T(struct measure_type m, double UR1, double URU, double UT1, double UTU)
{
    double x, G;
    G = Gain(REFLECTION_SPHERE, m, URU, 0);
    x = m.ad_t * (1 - m.at_t) * m.rw_t * Gain_22(m, URU, UTU);
    x *= (1 - m.f_r) * UT1 + (1 - m.at_r) * m.rw_r * m.as_r * UTU * (m.f_r * m.rw_r + (1 - m.f_r) * UR1) * G;
    return x;
}

void Calculate_Distance_With_Corrections(double UR1, double UT1,
    double Ru, double Tu, double URU, double UTU, double *M_R, double *M_T, double *dev)
{

    double UR1_calc, UT1_calc, URU_calc, UTU_calc;

    URU_calc = URU - MM.uru_lost;
    if (URU_calc < 0)
        URU_calc = 0;

    UTU_calc = UTU - MM.utu_lost;
    if (UTU_calc < 0)
        UTU_calc = 0;

    UR1_calc = UR1 - (1.0 - MM.fraction_of_ru_in_mr) * Ru - MM.ur1_lost;
    if (UR1_calc < 0)
        UR1_calc = 0;

    UT1_calc = UT1 - (1.0 - MM.fraction_of_tu_in_mt) * Tu - MM.ut1_lost;
    if (UT1_calc < 0)
        UT1_calc = 0;

    switch (MM.num_spheres) {
    case 0:

        {
            *M_R = UR1_calc;
            *M_T = UT1_calc;
        }

        break;

    case 1:
        if (MM.method == COMPARISON) {

            {
                *M_R = UR1_calc;
                *M_T = UT1_calc;
            }

        }
        else {

            double P_std, P, P_0, G, G_0, G_std, r_first, P_ss, P_su;
            r_first = 1;

            if (MM.baffle_r)
                r_first = MM.rw_r * (1 - MM.at_r);

            UR1_calc = UR1 - Ru - MM.ur1_lost;
            if (UR1_calc < 0)
                UR1_calc = 0;

            G_0 = Gain(REFLECTION_SPHERE, MM, 0.0, 0.0);
            G = Gain(REFLECTION_SPHERE, MM, URU_calc, 0.0);
            G_std = Gain(REFLECTION_SPHERE, MM, MM.rstd_r, 0.0);

            P_std = G_std * (MM.rstd_r * (1 - MM.f_r) + MM.f_r * MM.rw_r);
            P_0 = G_0 * (MM.f_r * MM.rw_r);
            P_ss = r_first * (UR1_calc * (1 - MM.f_r) + MM.f_r * MM.rw_r);
            P_su = MM.rw_r * (1 - MM.f_r) * MM.fraction_of_ru_in_mr * Ru;
            P = G * (P_ss + P_su);

            *M_R = MM.rstd_r * (P - P_0) / (P_std - P_0);

            if (Debug(DEBUG_SPHERE_GAIN) && !CALCULATING_GRID) {
                fprintf(stderr, "SPHERE: REFLECTION\n");
                fprintf(stderr, "SPHERE:      G0 = %7.3f      G  = %7.3f G_std = %7.3f\n", G_0, G, G_std);
                fprintf(stderr, "SPHERE:      P0 = %7.3f      P  = %7.3f P_std = %7.3f\n", P_0, P, P_std);
                fprintf(stderr, "SPHERE:     UR1 = %7.3f UR1calc = %7.3f   M_R = %7.3f\n", UR1, UR1_calc, *M_R);
            }

            {
                double r_first = 1;
                double r_third = MM.rstd_t;

                if (MM.fraction_of_tu_in_mt == 0)
                    r_third = 0;

                if (MM.baffle_t)
                    r_first = MM.rw_t * (1 - MM.at_t) + MM.rstd_t * MM.at_t;

                UT1_calc = UT1 - Tu - MM.ut1_lost;
                if (UT1_calc < 0)
                    UT1_calc = 0;

                G = Gain(TRANSMISSION_SPHERE, MM, URU_calc, r_third);
                G_std = Gain(TRANSMISSION_SPHERE, MM, 0, MM.rstd_t);

                *M_T = (r_third * Tu * MM.fraction_of_tu_in_mt + r_first * UT1_calc) * G / G_std;

                if (Debug(DEBUG_SPHERE_GAIN) && !CALCULATING_GRID) {
                    fprintf(stderr, "SPHERE: TRANSMISSION\n");
                    fprintf(stderr, "SPHERE:      G  = %7.3f   G_std = %7.3f\n", G, G_std);
                    fprintf(stderr, "SPHERE:     UT1 = %7.3f UT1calc = %7.3f T_c = %7.3f\n", UT1, UT1_calc, Tu);
                    fprintf(stderr, "SPHERE:     M_T = %7.3f\n", *M_T);
                    fprintf(stderr, "\n");
                }
            }

        }
        break;

    case 2:

        {
            double R_0, T_0;
            R_0 = Two_Sphere_R(MM, 0, 0, 0, 0);
            T_0 = Two_Sphere_T(MM, 0, 0, 0, 0);

            *M_R = MM.rstd_r * (Two_Sphere_R(MM, UR1_calc, URU_calc, UT1_calc, UTU_calc) - R_0) /
                (Two_Sphere_R(MM, MM.rstd_r, MM.rstd_r, 0, 0) - R_0);
            *M_T = (Two_Sphere_T(MM, UR1_calc, URU_calc, UT1_calc, UTU_calc) - T_0) /
                (Two_Sphere_T(MM, 0, 0, 1, 1) - T_0);
        }

        break;

    default:
        fprintf(stderr, "Bad number of spheres = %d\n", MM.num_spheres);
        exit(EXIT_FAILURE);
    }

    if (RR.search == FIND_A || RR.search == FIND_G || RR.search == FIND_B ||
        RR.search == FIND_Bs || RR.search == FIND_Ba) {

        if (MM.m_t > 0) {
            if (RR.metric == RELATIVE)
                *dev = fabs(MM.m_t - *M_T) / (MM.m_t + ABIT);
            else
                *dev = fabs(MM.m_t - *M_T);
        }
        else {
            if (RR.metric == RELATIVE)
                *dev = fabs(MM.m_r - *M_R) / (MM.m_r + ABIT);
            else
                *dev = fabs(MM.m_r - *M_R);
        }

    }
    else {

        if (RR.metric == RELATIVE) {
            if (MM.m_t > ABIT)
                *dev = T_TRUST_FACTOR * fabs(MM.m_t - *M_T) / (UTU_calc + ABIT);
            if (RR.default_a != 0) {
                *dev += fabs(MM.m_r - *M_R) / (URU_calc + ABIT);
            }
        }
        else {
            *dev = T_TRUST_FACTOR * fabs(MM.m_t - *M_T);
            if (RR.default_a != 0)
                *dev += fabs(MM.m_r - *M_R);
        }

    }

    if ((Debug(DEBUG_ITERATIONS) && !CALCULATING_GRID) || (Debug(DEBUG_GRID_CALC) && CALCULATING_GRID)) {
        fprintf(stderr, "%10.5f %10.4f %10.5f |", RR.slab.a, RR.slab.b, RR.slab.g);
        fprintf(stderr, " %10.5f %10.5f |", MM.m_r, *M_R);
        fprintf(stderr, " %10.5f %10.5f |", MM.m_t, *M_T);
        fprintf(stderr, "%10.3f\n", *dev);
    }

}

double Calculate_Grid_Distance(int i, int j)
{
    double ur1, ut1, uru, utu, Ru, Tu, b, dev, LR, LT;

    if (Debug(DEBUG_GRID_CALC) && i == 0 && j == 0) {
        fprintf(stderr, "+   i   j ");
        fprintf(stderr, "      a         b          g     |");
        fprintf(stderr, "     M_R        grid   |");
        fprintf(stderr, "     M_T        grid   |  distance\n");
    }

    if (Debug(DEBUG_GRID_CALC))
        fprintf(stderr, "g %3d %3d ", i, j);

    b = The_Grid[GRID_SIZE * i + j][B_COLUMN];
    ur1 = The_Grid[GRID_SIZE * i + j][UR1_COLUMN];
    ut1 = The_Grid[GRID_SIZE * i + j][UT1_COLUMN];
    uru = The_Grid[GRID_SIZE * i + j][URU_COLUMN];
    utu = The_Grid[GRID_SIZE * i + j][UTU_COLUMN];
    RR.slab.a = The_Grid[GRID_SIZE * i + j][A_COLUMN];
    RR.slab.b = The_Grid[GRID_SIZE * i + j][B_COLUMN];
    RR.slab.g = The_Grid[GRID_SIZE * i + j][G_COLUMN];

    Sp_mu_RT_Flip(MM.flip_sample,
        RR.slab.n_top_slide, RR.slab.n_slab, RR.slab.n_bottom_slide,
        RR.slab.b_top_slide, b, RR.slab.b_bottom_slide, RR.slab.cos_angle, &Ru, &Tu);

    CALCULATING_GRID = 1;
    Calculate_Distance_With_Corrections(ur1, ut1, Ru, Tu, uru, utu, &LR, &LT, &dev);
    CALCULATING_GRID = 0;

    return dev;
}

void Calculate_Distance(double *M_R, double *M_T, double *deviation)
{
    double Ru, Tu, ur1, ut1, uru, utu;

    if (RR.slab.b <= 1e-6)
        RR.slab.b = 1e-6;

    RT_Flip(MM.flip_sample, RR.method.quad_pts, &RR.slab, &ur1, &ut1, &uru, &utu);

    Sp_mu_RT_Flip(MM.flip_sample,
        RR.slab.n_top_slide, RR.slab.n_slab, RR.slab.n_bottom_slide,
        RR.slab.b_top_slide, RR.slab.b, RR.slab.b_bottom_slide, RR.slab.cos_angle, &Ru, &Tu);

    if ((!CALCULATING_GRID && Debug(DEBUG_ITERATIONS)) || (CALCULATING_GRID && Debug(DEBUG_GRID_CALC)))
        fprintf(stderr, "        ");

    Calculate_Distance_With_Corrections(ur1, ut1, Ru, Tu, uru, utu, M_R, M_T, deviation);
}

void abg_distance(double a, double b, double g, guess_type *guess)
{
    double m_r, m_t, distance;
    struct measure_type old_mm;
    struct invert_type old_rr;

    Get_Calc_State(&old_mm, &old_rr);

    RR.slab.a = a;
    RR.slab.b = b;
    RR.slab.g = g;

    Calculate_Distance(&m_r, &m_t, &distance);

    Set_Calc_State(old_mm, old_rr);

    guess->a = a;
    guess->b = b;
    guess->g = g;
    guess->distance = distance;
}

double Find_AG_fn(double x[])
{
    double m_r, m_t, deviation;
    RR.slab.a = acalc2a(x[1]);
    RR.slab.g = gcalc2g(x[2]);
    Calculate_Distance(&m_r, &m_t, &deviation);
    return deviation;
}

double Find_AB_fn(double x[])
{
    double m_r, m_t, deviation;
    RR.slab.a = acalc2a(x[1]);
    RR.slab.b = bcalc2b(x[2]);
    Calculate_Distance(&m_r, &m_t, &deviation);
    return deviation;
}

double Find_Ba_fn(double x)
{
    double m_r, m_t, deviation, ba, bs;

    bs = RR.slab.b;
    ba = bcalc2b(x);
    RR.slab.b = ba + bs;
    RR.slab.a = bs / (ba + bs);

    Calculate_Distance(&m_r, &m_t, &deviation);

    RR.slab.b = bs;
    return deviation;
}

double Find_Bs_fn(double x)
{
    double m_r, m_t, deviation, ba, bs;

    ba = RR.slab.b;
    bs = bcalc2b(x);
    RR.slab.b = ba + bs;
    RR.slab.a = bs / (ba + bs);

    Calculate_Distance(&m_r, &m_t, &deviation);

    RR.slab.b = ba;
    return deviation;
}

double Find_A_fn(double x)
{
    double m_r, m_t, deviation;
    RR.slab.a = acalc2a(x);
    Calculate_Distance(&m_r, &m_t, &deviation);
    return deviation;
}

double Find_B_fn(double x)
{
    double m_r, m_t, deviation;
    RR.slab.b = bcalc2b(x);
    Calculate_Distance(&m_r, &m_t, &deviation);
    return deviation;
}

double Find_G_fn(double x)
{
    double m_r, m_t, deviation;
    RR.slab.g = gcalc2g(x);
    Calculate_Distance(&m_r, &m_t, &deviation);
    return deviation;
}

double Find_BG_fn(double x[])
{
    double m_r, m_t, deviation;
    RR.slab.b = bcalc2b(x[1]);
    RR.slab.g = gcalc2g(x[2]);
    RR.slab.a = RR.default_a;
    Calculate_Distance(&m_r, &m_t, &deviation);
    return deviation;
}

double Find_BaG_fn(double x[])
{
    double m_r, m_t, deviation;

    RR.slab.b = bcalc2b(x[1]) + RR.default_bs;
    if (RR.slab.b <= 0)
        RR.slab.a = 0;
    else
        RR.slab.a = RR.default_bs / RR.slab.b;

    RR.slab.g = gcalc2g(x[2]);

    Calculate_Distance(&m_r, &m_t, &deviation);
    return deviation;
}

double Find_BsG_fn(double x[])
{
    double m_r, m_t, deviation;

    RR.slab.b = bcalc2b(x[1]) + RR.default_ba;
    if (RR.slab.b <= 0)
        RR.slab.a = 0;
    else
        RR.slab.a = 1.0 - RR.default_ba / RR.slab.b;

    RR.slab.g = gcalc2g(x[2]);
    Calculate_Distance(&m_r, &m_t, &deviation);
    return deviation;
}

double maxloss(double f)
{
    struct measure_type m_old;
    struct invert_type r_old;
    double m_r, m_t, deviation;

    Get_Calc_State(&m_old, &r_old);

    RR.slab.a = 1.0;
    MM.ur1_lost *= f;
    MM.ut1_lost *= f;

    Calculate_Distance(&m_r, &m_t, &deviation);

    Set_Calc_State(m_old, r_old);
    deviation = ((MM.m_r + MM.m_t) - (m_r + m_t));

    return deviation;
}

void Max_Light_Loss(struct measure_type m, struct invert_type r, double *ur1_loss, double *ut1_loss)
{
    struct measure_type m_old;
    struct invert_type r_old;

    *ur1_loss = m.ur1_lost;
    *ut1_loss = m.ut1_lost;

    if (Debug(DEBUG_LOST_LIGHT))
        fprintf(stderr, "\nlost before ur1=%7.5f, ut1=%7.5f\n", *ur1_loss, *ut1_loss);

    Get_Calc_State(&m_old, &r_old);

    Set_Calc_State(m, r);

    if (maxloss(1.0) * maxloss(0.0) < 0) {
        double frac;
        frac = zbrent(maxloss, 0.00, 1.0, 0.001);

        *ur1_loss = m.ur1_lost * frac;
        *ut1_loss = m.ut1_lost * frac;
    }

    Set_Calc_State(m_old, r_old);
    if (Debug(DEBUG_LOST_LIGHT))
        fprintf(stderr, "lost after  ur1=%7.5f, ut1=%7.5f\n", *ur1_loss, *ut1_loss);
}
