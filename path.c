/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2023 Chenxi Zhou <chnx.zhou@gmail.com>                          *
 *                                                                               *
 * Permission is hereby granted, free of charge, to any person obtaining a copy  *
 * of this software and associated documentation files (the "Software"), to deal *
 * in the Software without restriction, including without limitation the rights  *
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
 * copies of the Software, and to permit persons to whom the Software is         *
 * furnished to do so, subject to the following conditions:                      *
 *                                                                               *
 * The above copyright notice and this permission notice shall be included in    *
 * all copies or substantial portions of the Software.                           *
 *                                                                               *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE *
 * SOFTWARE.                                                                     *
 *********************************************************************************/

/********************************** Revision History *****************************
 *                                                                               *
 * 03/02/23 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <zlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "khashl.h"
#include "kstring.h"
#include "kseq.h"
#include "kvec.h"
#include "kdq.h"

#include "path.h"
#include "graph.h"
#include "misc.h"

#undef DEBUG_SRCC
#undef DEBUG_SEG_COV_EST
#undef DEBUG_SEG_COPY
#undef DEBUG_SEG_COPY_EST
#undef DEBUG_SEG_COV_BOUND
#undef DEBUG_SEG_COV_ADJUST
#undef DEBUG_BRUTE_FORCE_OPTIM
#undef DEBUG_SIM_ANNEAL_OPTIM
#undef DEBUG_PATH_FINDER

// Minimum kmer coverage to copy_number 1, based on prior knowledge.
extern int global_avg_cov;

void path_destroy(path_t *path)
{
    if (!path) return;
    if (path->sid)
        free(path->sid);
    if (path->v)
        free(path->v);
}

void path_v_destroy(path_v *path)
{
    if (!path) return;
    size_t i;
    for (i = 0; i < path->n; ++i)
        path_destroy(&path->a[i]);
    free(path->a);
}

static int u64_cmpfunc(const void *a, const void *b)
{
    return (*(uint64_t *) a > *(uint64_t *) b) - (*(uint64_t *) a < *(uint64_t *) b);
}

uint64_t asg_seg_len(asg_t *asg)
{
    asmg_t *g = asg->asmg;
    uint64_t i, seg_len = 0;
    for (i = 0; i < g->n_vtx; ++i) {
        if (g->vtx[i].del) continue;
        seg_len += g->vtx[i].len;
    }
    return seg_len;
}

static double graph_sequence_coverage_lower_bound(asg_t *asg, double cov_nq)
{
    uint32_t i, n_seg, len, cov;
    uint64_t *seqcovs, tot_seg_len, tot_cov, len_thresh;
    double cov_bound;
    asmg_t *g;

    g = asg->asmg;
    n_seg = 0;
    tot_seg_len = 0;
    MYMALLOC(seqcovs, g->n_vtx);
    MYBONE(seqcovs, g->n_vtx);
    for (i = 0; i < g->n_vtx; ++i) {
        if (g->vtx[i].del) continue;
        len = g->vtx[i].len;
        cov = g->vtx[i].cov;
        
        ++n_seg;
        tot_seg_len += len;
        seqcovs[i] = (uint64_t) cov << 32 | len;
    }

    qsort(seqcovs, g->n_vtx, sizeof(uint64_t), u64_cmpfunc);
    len_thresh = tot_seg_len * cov_nq;
    i = 0;
    len = (uint32_t) seqcovs[i];
    tot_seg_len = tot_cov = 0;
    while (tot_seg_len + len <= len_thresh) {
        tot_cov += (seqcovs[i] >> 32) * len;
        tot_seg_len += len;
        len = (uint32_t) seqcovs[++i];
    }
    if (tot_seg_len < len_thresh)
        tot_cov += (seqcovs[i] >> 32) * (len_thresh - tot_seg_len);
    cov_bound = (double) tot_cov / len_thresh;
    free(seqcovs);

    cov_bound *= 1 - cov_nq;

#ifdef DEBUG_SEG_COV_BOUND
    fprintf(stderr, "[DEBUG_SEG_COV_BOUND::%s] estimated sequence coverage lower boundary: %.3f\n", __func__, cov_bound);
#endif

    return cov_bound;
}

/*** old implementation
 *
static double graph_sequence_coverage_rough(asg_t *asg)
{
    uint32_t i, len, cov, deg_in, deg_out;
    double tot_seg_len, tot_cov, avg_cov;
    asmg_t *g;

    g = asg->asmg;
    tot_seg_len = tot_cov = 0;
    for (i = 0; i < g->n_vtx; ++i) {
        if (g->vtx[i].del) continue;
        len = g->vtx[i].len;
        cov = g->vtx[i].cov;
        avg_cov = cov;
        deg_out = MAX(1, asmg_arc_n1(g, i<<1));
        deg_in  = MAX(1, asmg_arc_n1(g, i<<1|1));
        avg_cov = avg_cov * 2 / (deg_out + deg_in);
#ifdef DEBUG_SEG_COV_EST
        // seg and vtx indices are always interchangeable as we never do graph clean
        fprintf(stderr, "[DEBUG_SEG_COV_EST::%s] %s %u %u [adj out deg: %u] [adj in deg: %u] [normalized: %.3f]\n",
                __func__, asg->seg[i].name, len, cov, deg_out, deg_in, avg_cov);
#endif
        tot_seg_len += len;
        tot_cov += avg_cov * len;
    }
    
    if (tot_seg_len == 0) return 0;

    avg_cov = tot_cov / tot_seg_len;

#ifdef DEBUG_SEG_COV_EST
    fprintf(stderr, "[DEBUG_SEG_COV_EST::%s] estimated sequence coverage: %.3f\n", __func__, avg_cov);
#endif

    return avg_cov;
}
**/

static double graph_sequence_coverage_rough(asg_t *asg, double min_cf)
{
    uint32_t i, j, len, cov, best1;
    double tot_len, tot_len_c, tot_rm, avg_cov, near1, diff1;
    asmg_t *g;
    kvec_t(uint64_t) lc_p;

    g = asg->asmg;
    kv_init(lc_p);
    for (i = 0; i < g->n_vtx; ++i) {
        if (g->vtx[i].del) continue;
        kv_push(uint64_t, lc_p, (uint64_t) (g->vtx[i].cov) << 32 | (g->vtx[i].len));
    }
    if (lc_p.n == 0) return .0;

    // Lowest coverage first.
    qsort(lc_p.a, lc_p.n, sizeof(uint64_t), u64_cmpfunc);

    best1 = 0;
    near1 = DBL_MAX;
    for (i = 0; i < lc_p.n; ++i) {
        avg_cov = (double) (lc_p.a[i] >> 32);
        if (avg_cov == 0)
            continue;

        // Compute total length of all remaining nodes from here
        tot_len = tot_len_c = tot_rm = 0;
        for (j = 0; j < lc_p.n; ++j) {
            len = (uint32_t) lc_p.a[j];
            cov = lc_p.a[j] >> 32;
            if (cov / avg_cov >= min_cf) {
                tot_len += len;
                tot_len_c += (double) len * cov / avg_cov;
            } else {
                tot_rm += len;
            }
        }

        // TODO is this good?
        if (tot_rm / (tot_rm + tot_len) > 0.7) break;

        if (tot_len > 0) {
            diff1 = fabs(tot_len_c / tot_len - 1.0);
            if (diff1 < near1) {
                near1 = diff1;
                best1 = i;
            }
#ifdef DEBUG_SEG_COV_EST
            fprintf(stderr, "[DEBUG_SEG_COV_EST::%s] seg %u [%u %lu] - len: %.0f; len_rm: %.0f; len_c: %.3f; diff1: %.3f; best: [%u %.3f]\n",
                    __func__, i, (uint32_t) lc_p.a[i], lc_p.a[i]>>32, tot_len, tot_rm, tot_len_c, diff1, best1, near1);
#endif
        }
    }
    if (near1 == DBL_MAX) {
        kv_destroy(lc_p);
        return .0;
    }

    avg_cov = (double) (lc_p.a[best1] >> 32);

#ifdef DEBUG_SEG_COV_EST
    fprintf(stderr, "[DEBUG_SEG_COV_EST::%s] estimated sequence coverage: %.3f\n", __func__, avg_cov);
#endif

    kv_destroy(lc_p);
    return avg_cov;
}

static void make_seg_dups(asg_t *asg, kh_u32_t *seg_dups, uint32_t s, uint32_t copy)
{
    // copy number include sequence itself
    // i.e., make 'copy' copies and remove the original
    // FIXME arcs involving self cycles are not copied since it will potentially cause problems
    // for example
    // in a graph containing path (X+)<->(a+,b+,b+,c+)<->(Y+)
    // the self cycle b+,b+ will increase the ambiguities of the path selection
    // valid path including
    // (X+)->(a+,b+,b+,c+)->(Y+)->(a-,b-,b-,c-)
    // (X+)->(a+,b+,b+,b+,c+)->(Y+)->(a-,b-,c-)
    // TODO solve the copy number of self cycles in postprocessing
    // now the self arcs is processed in this way
    // (X+)->(a+)[4]->(Y+) => (X+)->(a0+) ->(a1+) ->(a2+) ->(a3+) ->(Y+)
    //                            ->(a1+) ->(a2+) ->(a0+) ->(a1+) ->(a0+)
    //                            ->(a2+) ->(a3+) ->(a3+) ->(a2+) ->(a1+)
    //                            ->(a3+) ->(Y+)  ->(Y+)  ->(Y+)  ->(a2+)
    // (a+)->(a-) arcs are not copied
    uint32_t sid;
    uint64_t i, j, v, nv, pv, self_arc;
    kvec_t(uint64_t) arcs_diff;
    asmg_t *g;
    asmg_arc_t *av;
    asmg_vtx_t *vt;
    char *name;
    asg_seg_t *seg, *seg_c;
    khint32_t k;
    int absent;

    g = asg->asmg;
    // mark arcs from vtx
    kv_init(arcs_diff);
    self_arc = UINT64_MAX;

    for (i = 0; i < 2; ++i) {
        v = s<<1 | i;
        pv = g->idx_p[v];
        nv = asmg_arc_n(g, v);
        av = asmg_arc_a(g, v);
        for (j = 0; j < nv; ++j) {
            if (av[j].del) continue;
            if ((av[j].v >> 1) != (av[j].w >> 1))
                kv_push(uint64_t, arcs_diff, pv+j);
            else if (av[j].v == av[j].w && i == 0)
                // (a+)->(a+)
                // avoid adding v+->v+ and v-->v- twice
                self_arc = pv + j;
        }
    }

    for (i = 0; i < copy; ++i) {
        // make a copy of the segment
        seg = &asg->seg[s];
        MYMALLOC(name, strlen(seg->name) + floor(log10(abs(i + 1))) + 7);
        sprintf(name, "%s_copy%lu", seg->name, i);
        sid = asg_add_seg(asg, name, 0);
        free(name);
        // add to seg dups map
        k = kh_u32_put(seg_dups, sid, &absent);
        kh_val(seg_dups, k) = s;
        // only copy essential fields
        seg = &asg->seg[s]; // need to redo this as the seg might be reallocated
        seg_c = &asg->seg[sid];
        seg_c->len = seg->len;
        seg_c->cov = seg->cov;

        // add vtx
        asmg_vtx_addp(g, &vt);
        MYBZERO(vt, 1);
        vt->len = seg->len;
        // using vtx coverage instead of seg coverage
        //vt->cov = seg->cov / copy;
        vt->cov = g->vtx[s].cov / copy;

        // make copies of links
        // arc cov already changed
        for (j = 0; j < arcs_diff.n; ++j) {
            av = &g->arc[arcs_diff.a[j]];
            asmg_arc_add2(g, sid << 1 | (av->v&1), av->w, av->ln, av->ls, UINT64_MAX, av->cov / copy, av->comp);
        }
        
        if (self_arc != UINT64_MAX) {
            // tandem repeat
            // self arcs are added between copies
            av = &g->arc[self_arc];
            for (j = 0; j < i; ++j) {
                asmg_arc_add2(g, (sid - i + j) << 1, sid << 1, av->ln, av->ls, UINT64_MAX, av->cov / copy, 0);
                asmg_arc_add2(g, sid << 1, (sid - i + j) << 1, av->ln, av->ls, UINT64_MAX, av->cov / copy, 0);
            }
        }
    }

    kv_destroy(arcs_diff);

    // need to redo sort and index but not graph clean
    asmg_finalize(g, 0);

    // delete the original vtx
    asmg_vtx_del(g, s, 1);

    return;
}

#define EM_MAX_ITER 1000

// Any copy_number 0 nodes that have used arc data has their copy_number
// set to 1 anyway.  This can happen where the kmer data wasn't strong
// enough, but we've observed some transitions.  Permit minor usage count.
void graph_sequence_coverage_edge_rescue(asg_t *asg, int *copy_number) {
    asmg_t *g = asg->asmg;
    uint32_t n_seg = asg->n_seg;

    for (int dir = 0; dir < 2; dir++) {
        for (int i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            int v = (i<<1) | dir;
            int nv = asmg_arc_n(g, v);
            asmg_arc_t *av = asmg_arc_a(g, v);
            for (int j = 0; j < nv; j++) {
                if (av[j].del) continue;
                if (!av[j].cov) continue;
                int from = av[j].v>>1;
                int to   = av[j].w>>1;
                //fprintf(stderr, "%s->%s ec %d\n", asg->seg[from].name, asg->seg[to].name, av[j].cov);
                if (copy_number[from] && !copy_number[to]) {
                    fprintf(stderr, "Incr copy in %s\n", asg->seg[to].name);
                    copy_number[to]=1;
                }
                if (copy_number[to] && !copy_number[from]) {
                    fprintf(stderr, "Incr copy in %s\n", asg->seg[from].name);
                    copy_number[from]=1;
                }
            }
        }
    }
}

// Similarly look for trivial bubbles A->{P,Q}->B where both paths are
// unknown and A and B are in use.  Our choice is to break the path up
// and give a sub-string, or guess (which is what this code encourages).
//
// However, if one of P and Q *is* known, then delete the other node
// and force it to never to be used.
int graph_sequence_coverage_pick_bubbles(asg_t *asg, int *copy_number) {
    asmg_t *g = asg->asmg;
    uint32_t n_seg = asg->n_seg;
 
    // Similarly look for trivial bubbles A->{P,Q}->B where both paths are
    // unknown and A and B are in use.  Our choice is to break the path up
    // and give a sub-string, or guess (which is what this code encourages).
    //
    // However, if one of P and Q *is* known, then delete the other node
    // and force it to never to be used.

    // Step 1: allocate and fill out an adjacency matrix
    int **link_to   = malloc(n_seg * sizeof(int *));
    int **link_from = malloc(n_seg * sizeof(int *));
    if (!link_to || !link_from)
        return -1; // mem leak, but exiting
    int err = 0;
    for (int i = 0; i < n_seg; i++) {
        link_to[i]   = calloc(n_seg, sizeof(int));
        link_from[i] = calloc(n_seg, sizeof(int));
        err |= (!link_to[i] || !link_from[i]);
    }
    if (err)
        return -1; // mem leak, but exiting

    for (int dir = 0; dir < 2; dir++) {
        for (int i = 0; i < n_seg; i++) {
            if (g->vtx[i].del) continue;
            int vtx_cov = g->vtx[i].cov;
            int v = (i<<1) | dir;
            int nv = asmg_arc_n(g, v);
            asmg_arc_t *av = asmg_arc_a(g, v);
            for (int j = 0; j < nv; j++) {
                if (av[j].del) continue;

                //if (!av[j].cov) continue;
            
                int from = av[j].v>>1;
                int to   = av[j].w>>1;

                printf("Edge %d -> %d\n", from, to);
 
                if (vtx_cov)
                    link_from[from][to]++;
                
                if (!g->vtx[to].del && g->vtx[to].cov)
                    link_to[from][to]++;
            }
        }

        break; // one dir only to remove A>B>C and C>B>A dups
    }

    // Step 2: Identify trivial bubble routes to rescue.
    for (int l = 0; l < n_seg; l++) {
        for (int r = 0; r < n_seg; r++) {
            if (l == r) continue;
            int found = 0, routes = 0, bypass;
            //bypass = (link_from[l][r] || link_to[l][r]);
            bypass = 0; // the above isn't helping.
            for (int i = 0; i < n_seg; i++) {
                if (i == l || i == r) continue;
                if (!link_from[l][i] || !link_to[i][r]) continue;
                printf("Link %d -> %d -> %d, covered = %d,%d\n",
                       l, i, r, found, (int)copy_number[i]);
                found += copy_number[i]>0;
                routes++;
            }

            if (routes && !found && !bypass) {
                for (int i = 0; i < n_seg; i++) {
                    if (i == l || i == r) continue;
                    if (!link_from[l][i] || !link_to[i][r]) continue;
                    printf("Add copy to middle %d -> %d -> %d\n", l, i, r);
                    copy_number[i] = 1;
                }
            } else if (routes && (found || bypass)) {
                // Explicitly mark some nodes as being excluded as we have
                // an observed path.
                for (int i = 0; i < n_seg; i++) {
                    if (i == l || i == r) continue;
                    if (!link_from[l][i] || !link_to[i][r]) continue;
                    //copy_number[i] = -1;
                    g->vtx[i].del = 1;
                }
            }
        }
    }

    for (int i = 0; i < n_seg; i++) {
        free(link_to[i]);
        free(link_from[i]);
    }
    free(link_to);
    free(link_from);

    return 0;
}

void graph_sequence_coverage_add_neighbours(asg_t *asg, int *copy_number,
                                            int min_length) {
    asmg_t *g = asg->asmg;
    uint32_t n_seg = asg->n_seg;

    for (int dir = 0; dir < 2; dir++) {
        for (int i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            int v = (i<<1) | dir;
            int nv = asmg_arc_n(g, v);
            asmg_arc_t *av = asmg_arc_a(g, v);
            for (int j = 0; j < nv; j++) {
                if (av[j].del) continue;
                int from = av[j].v>>1;
                int to   = av[j].w>>1;
                //fprintf(stderr, "%s->%s ec %d\n", asg->seg[from].name, asg->seg[to].name, av[j].cov);
                if (copy_number[from] && !copy_number[to] && g->vtx[to].len < min_length) {
                    fprintf(stderr, "Incr copy in %s\n", asg->seg[to].name);
                    copy_number[to]=1;
                }
                if (copy_number[to] && !copy_number[from] && g->vtx[from].len < min_length) {
                    fprintf(stderr, "Incr copy in %s\n", asg->seg[from].name);
                    copy_number[from]=1;
                }
            }
        }
    }
}

double graph_sequence_coverage_precise(asg_t *asg, double min_cf, int min_copy, int max_copy, int edge_to_seq, int bub_check, int neighbour_steps, int min_neighbour_len, int **_copy_number)
{
    uint32_t i, n_seg, iter;
    int *copy_number;
    double total_covs, total_lens, avg_cov, new_avg_cov, min_avg_cov;
    asmg_t *g;
    
    n_seg = asg->n_seg;
    g = asg->asmg;
    min_avg_cov = graph_sequence_coverage_lower_bound(asg, 0.3);
    avg_cov = graph_sequence_coverage_rough(asg, min_cf);
    if (avg_cov < global_avg_cov)
        avg_cov = global_avg_cov;
#ifdef DEBUG_SEG_COPY_EST
    fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] min coverage: %.3f; rough avg coverage: %.3f\n",
            __func__, min_avg_cov, avg_cov);
#endif
    avg_cov = MAX(avg_cov, min_avg_cov);
    MYCALLOC(copy_number, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        if (global_avg_cov)
            copy_number[i] = MIN(MAX(min_copy, ceil((double) (g->vtx[i].cov) / avg_cov)), max_copy);
        else
            copy_number[i] = MIN(MAX(min_copy, lround((double) (g->vtx[i].cov) / avg_cov)), max_copy);
    }

    iter = 0;
    while (iter++ < EM_MAX_ITER) {
#ifdef DEBUG_SEG_COPY_EST
        fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] iteration %u avg coverage: %.3f\n", __func__, iter, avg_cov);
#endif
        total_covs = total_lens = 0;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            total_lens += (double) g->vtx[i].len * copy_number[i];
            total_covs += (double) g->vtx[i].len * g->vtx[i].cov;
        }
        // FIXME total_lens could be zero
        new_avg_cov = total_lens < FLT_EPSILON? DBL_MAX : (total_covs / total_lens);
        new_avg_cov = MAX(new_avg_cov, min_avg_cov);
        if (fabs(new_avg_cov - avg_cov) < FLT_EPSILON) 
            break; // converged
        if (avg_cov < global_avg_cov)
            avg_cov = global_avg_cov;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            if (global_avg_cov)
                copy_number[i] = MIN(MAX(min_copy, ceil((double) (g->vtx[i].cov) / avg_cov)), max_copy);
            else
                copy_number[i] = MIN(MAX(min_copy, lround((double) (g->vtx[i].cov) / avg_cov)), max_copy);
        }
    }

    // Any copy_number 0 nodes that have used arc data has their copy_number
    // set to 1 anyway.  This can happen where the kmer data wasn't strong
    // enough, but we've observed some transitions.  Permit minor usage count.
    if (edge_to_seq)
        graph_sequence_coverage_edge_rescue(asg, copy_number);

    if (bub_check)
        if (graph_sequence_coverage_pick_bubbles(asg, copy_number) < 0)
            abort();

    // Final method.  Neighbour boosting.  The benefit to doing this now
    // is we can still recognise the deleted node assessment from the previous
    // step.  Add coverage to x in A->x and x->A if A has coverage.
    for (int i = 0; i < neighbour_steps; i++)
        graph_sequence_coverage_add_neighbours(asg, copy_number, min_neighbour_len);

#ifdef DEBUG_SEG_COPY_EST
    fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] sequence copy number estimation finished in %u iterations with an average sequence coverage %.3f\n",
            __func__, iter, avg_cov);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        fprintf(stderr, "[DEBUG_SEG_COPY_EST::%s] %s %u %d\n", __func__, asg->seg[i].name, g->vtx[i].cov, copy_number[i]);
    }
#endif
    
    if (_copy_number)
        *_copy_number = copy_number;
    else
        free(copy_number);

    return avg_cov;
}

#define DVAL(var) ((var)->D)
// with BALANCE_IN_OUT defined
// the objective function considers balanced indegree and outdegree
#define BALANCE_IN_OUT
#ifdef BALANCE_IN_OUT
#define FVAL(fun) do { \
    int __i, __n = (fun)->N; \
    double __val[2] = {.0, .0}; \
    for (__i = 0; __i < __n; ++__i) \
        __val[(fun)->V[__i]&1] += DVAL((fun)->VAR[(fun)->V[__i]>>1]); \
    (fun)->VAL = (fun)->weight * (fabs((fun)->v_exp - __val[0]) / 2 + \
            fabs((fun)->v_exp - __val[1]) / 2 + \
            fabs(__val[0] - __val[1])); \
} while (0)
#else
#define FVAL(fun) do { \
    int __i, __n = (fun)->N; \
    double __val = (fun)->v_exp; \
    for (__i = 0; __i < __n; ++__i) \
        __val -= DVAL((fun)->VAR[(fun)->V[__i]]); \
    (fun)->VAL = (fun)->weight * __val * __val; \
} while (0)
#endif

typedef struct var {
    int B, D;
    struct var *prev, *next;
} var_t;

typedef struct {
    double weight; // weight
    double v_exp; // expected value
    var_t **VAR;
    int N, *V;
    double VAL;
} func_t;

#define BRUTE_FORCE_N_LIM 100000000

static inline double fvals(func_t *funcs, int n_func)
{
    int i;
    double fval = 0.;
    for (i = 0; i < n_func; ++i) {
        FVAL(&funcs[i]);
        fval += funcs[i].VAL;
    }
    return fval;
}

static inline void copy_results(var_t **vars, int n_var, int *res)
{
    int i;
    for (i = 0; i < n_var; ++i)
        res[i] = vars[i]->D;
}

// brute force optimization
static double estimate_arc_copy_number_brute_force_impl(func_t *funcs, int n_func, var_t **vars, int n_var, int *res, int64_t sol_space_size)
{
    double fval, m_fval;
    int64_t sol;
    int a, v;
    fval = fvals(funcs, n_func);
    m_fval = fval;
    copy_results(vars, n_var, res);
    sol = 0;
    while (++sol < sol_space_size) {
        a = 1, v = 0;
        while (a) {
            vars[v] = vars[v]->next;
            a = !vars[v]->B;
            ++v;
        }
        fval = fvals(funcs, n_func);
        if (fval < m_fval) {
            m_fval = fval;
            copy_results(vars, n_var, res);
        }
        if (fabs(m_fval) < FLT_EPSILON)
            break;
    }
#ifdef DEBUG_BRUTE_FORCE_OPTIM
    fprintf(stderr, "[DEBUG_BRUTE_FORCE_OPTIM::%s] brute force search finished after %ld/%ld attempts with a minimum fval: %.6f\n",
            __func__, sol, sol_space_size, m_fval);
#endif
    return m_fval;
}

#define SA_TEMPERATURE  1000
#define SA_COOLING_RATE .999
#define SA_MAX_ATTEMPTS 100
#define SA_RESTART_TEMP .99

static inline void random_walk_to_neighbour(var_t **var)
{
    if (rand() < RAND_MAX>>1) {
        // move to prev
        *var = (*var)->B == 0? (*var)->next : (*var)->prev;
    } else {
        // move to next
        *var = (*var)->next->B == 0? (*var)->prev : (*var)->next;
    }
}

static void set_vars(var_t **vars, int n_var, int *res)
{
    int i;
    for (i = 0; i < n_var; ++i) {
        while (vars[i]->D != res[i])
            vars[i] = vars[i]->next;
    }
}

// simulated annealing optimization
static double estimate_arc_copy_number_siman_impl(func_t *funcs, int n_func, var_t **vars, int n_var, int *res)
{
    int i, iter, n_iter;
    double optim_cost, current_cost, new_cost, temp0, temp, p;
    var_t *var;

    current_cost = fvals(funcs, n_func);
    optim_cost = current_cost;
    copy_results(vars, n_var, res);
    
    srand(1234);
    temp0 = SA_TEMPERATURE;
    n_iter = 0;
    for (iter = 0; iter < SA_MAX_ATTEMPTS; ++iter) {
        temp = temp0;
        while (temp > 1e-6) {
            // random select a var to update 
            i = rand() % n_var;
            // record the old var
            var = vars[i];
            // take a random walk to either prev or next
            random_walk_to_neighbour(&vars[i]);
            // calculate new cost
            new_cost = fvals(funcs, n_func);
            // record optim solution
            if (new_cost < optim_cost) {
                optim_cost = new_cost;
                copy_results(vars, n_var, res);
            }
            // acceptance probability
            p = exp(-(new_cost - current_cost) / temp);
            if (new_cost < current_cost || (double) rand() / RAND_MAX < p) {
                // accept update and update cost
                current_cost = new_cost;
            } else {
                // rollback to reject the update
                vars[i] = var;
            }
            // cooling down
            temp *= SA_COOLING_RATE;
            ++n_iter;
        }
        if (optim_cost == 0) break;
        // continue searching from the best solution so far
        temp0 *= SA_RESTART_TEMP;
        set_vars(vars, n_var, res);
#ifdef DEBUG_SIM_ANNEAL_OPTIM
        fprintf(stderr, "[DEBUG_SIM_ANNEAL_OPTIM::%s] optimization after %d iteration [%d attempts] with a minimum fval: %.6f\n",
                __func__, iter, n_iter, optim_cost);
#endif

    }
#ifdef DEBUG_SIM_ANNEAL_OPTIM
    fprintf(stderr, "[DEBUG_SIM_ANNEAL_OPTIM::%s] simulated annealing search finished after %d attempts with a minimum fval: %.6f\n",
            __func__, n_iter, optim_cost);
#endif
    return optim_cost;
}

int adjust_sequence_copy_number_by_graph_layout(asg_t *asg, double seq_coverage, double *_adjusted_cov, int *copy_number, int max_copy, int max_round)
{
    uint32_t i, j, n_seg, n_group, a_g, *arc_group;
    uint64_t link_id;
    int updated, round;
    asmg_t *g;
    asmg_arc_t *a;

    if (_adjusted_cov)
        *_adjusted_cov = seq_coverage;
    if (max_round == 0)
        max_round = 1;

    updated = 0;

    g = asg->asmg;
    n_seg = asg->n_seg;
    arc_group = asmg_uext_arc_group(g, &n_group);

#ifdef DEBUG_SEG_COV_ADJUST
    fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] identified %u arc groups\n", __func__, n_group);
    for (i = 0; i < g->n_arc; ++i) {
        a = &g->arc[i];
        if (a->del) continue;
        link_id = a->link_id;
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s%c -> %s%c arc group %u\n", __func__,
                asg->seg[a->v>>1].name, "+-"[a->v&1],
                asg->seg[a->w>>1].name, "+-"[a->w&1],
                arc_group[link_id]);
    }
#endif

    if (n_group == 0) {
        free(arc_group);
        return 0;
    }

#ifdef DEBUG_SEG_COV_ADJUST
    fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] sequence copy number BEFORE adjusted by graph layout\n", __func__);
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s %lu %u %d\n", __func__, 
                asg->seg[i].name, g->vtx[i].len, g->vtx[i].cov, copy_number[i]);
    }
#endif

    uint32_t *arc_copy_lb, *arc_copy_ub, vlb, wlb, lb, ub;
    MYCALLOC(arc_copy_lb, n_group); // arc copy number lower bound
    MYCALLOC(arc_copy_ub, n_group); // arc copy number upper bound

    // calculate lower and upper boundary of arc copy number for each arc group
    for (i = 0; i < g->n_arc; ++i) {
        a = &g->arc[i];
        if (a->del) continue;
        link_id = a->link_id;
        a_g = arc_group[link_id];
        // if a is the only outgoing/incoming edge from/to v/w
        // then lower bound is the v/w copy number
        // otherwise 0
        // upper bound is always the v/w copy number
        vlb = asmg_arc_n1(g, a->v) == 1? copy_number[a->v>>1] : 0;
        wlb = asmg_arc_n1(g, a->w^1) == 1? copy_number[a->w>>1] : 0;
        lb = MIN(vlb, wlb);
        ub = MAX(copy_number[a->v>>1], copy_number[a->w>>1]);
        // relax boundary by a factor of 1/3 - at least one copy
        lb = (uint32_t) ((double) lb * 2 / 3);
        ub = (uint32_t) ((double) ub * 4 / 3) + 1;
        ub = MIN(ub, (uint32_t) max_copy);
        arc_copy_lb[a_g] = MIN(lb, arc_copy_lb[a_g]);
        arc_copy_ub[a_g] = MAX(ub, arc_copy_ub[a_g]);
    }

    // make variable list
    var_t **VAR;
    MYCALLOC(VAR, n_group);
    for (i = 0; i < n_group; ++i) {
        var_t *var0;
        MYCALLOC(var0, 1);
        lb = arc_copy_lb[i];
        ub = arc_copy_ub[i];
        var0->B = 0;
        var0->D = lb;
        VAR[i] = var0;
        for (j = 1; j <= ub - lb; ++j) {
            var_t *var1;
            MYCALLOC(var1, 1);
            var1->B = j;
            var1->D = lb + j;
            var0->next = var1;
            var1->prev = var0;
            var0 = var1;
        }
        var0->next = VAR[i];
        VAR[i]->prev = var0;
    }

    // make objective function list
    kvec_t(func_t) funcs;
    kv_init(funcs);
    func_t *FUN;
    uint32_t k, v, na;
    int *funcmap; // the objective function index for each seg
    asmg_arc_t *av;
    MYMALLOC(funcmap, n_seg);
    for (i = 0; i < n_seg; ++i) {
        funcmap[i] = -1;
        if (g->vtx[i].del)
            continue;
#ifdef BALANCE_IN_OUT
        // with BALNCE_IN_OUT defined
        // the incoming and outgoing arcs of a vertex
        // are put in the same objective function unit
        // and marked by the last bit of V
        kvec_t(int) V = {0, 0 ,0};
        for (k = 0; k < 2; ++k) {
            v = i << 1 | k;
            na = asmg_arc_n(g, v);
            av = asmg_arc_a(g, v);
            for (j = 0; j < na; ++j) {
                if (av[j].del)
                    continue;
                link_id = av[j].link_id;
                a_g = arc_group[link_id];
                assert(a_g != (uint32_t) -1);
                kv_push(int, V, (int) (a_g << 1 | k));
            }
        }
        if (V.n) {
            MYREALLOC(V.a, V.n);
            funcmap[i] = funcs.n;
            kv_pushp(func_t, funcs, &FUN);
            FUN->weight = log10(g->vtx[i].len);
            FUN->v_exp = g->vtx[i].cov / seq_coverage; // copy_number[i];
            FUN->VAR = VAR;
            FUN->N = V.n;
            FUN->V = V.a;
            FUN->VAL = 0.;
        }
#else
        for (k = 0; k < 2; ++k) {
            v = i << 1 | k;
            na = asmg_arc_n(g, v);
            av = asmg_arc_a(g, v);
            kvec_t(int) V = {0, 0 ,0};
            for (j = 0; j < na; ++j) {
                if (av[j].del)
                    continue;
                link_id = av[j].link_id;
                a_g = arc_group[link_id];
                assert(a_g != (uint32_t) -1);
                kv_push(int, V, (int) a_g);
            }
            if (V.n) {
                MYREALLOC(V.a, V.n);
                funcmap[i] = funcs.n;
                kv_pushp(func_t, funcs, &FUN);
                FUN->weight = log10(g->vtx[i].len);
                FUN->v_exp = g->vtx[i].cov / seq_coverage; // copy_number[i];
                FUN->VAR = VAR;
                FUN->N = V.n;
                FUN->V = V.a;
                FUN->VAL = 0.;
            }
        }
#endif
    }

    // do optimization
    int *arc_copy;
    double adjusted_cov, min_avg_cov;
    double optim_cost, curr_cost;
    int64_t sol_space_size;
    int new_copy[2];
    
    min_avg_cov = graph_sequence_coverage_lower_bound(asg, 0.3);
    adjusted_cov = seq_coverage;

    MYCALLOC(arc_copy, n_group);
    sol_space_size = 1;
    for (i = 0; i < n_group && sol_space_size <= BRUTE_FORCE_N_LIM; ++i)
        sol_space_size *= (arc_copy_ub[i] - arc_copy_lb[i] + 1);
    
    optim_cost = fvals(funcs.a, funcs.n);
    round = 0;
    while (round++ < max_round) {
#ifdef DEBUG_SEG_COV_ADJUST
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] adjusting copy number round %d of %d\n", __func__, round, max_round);
#endif
        if (sol_space_size <= BRUTE_FORCE_N_LIM) {
            // do brute force optimization
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] run brute force searching: %ld\n",
                    __func__, sol_space_size);
#endif
            curr_cost = estimate_arc_copy_number_brute_force_impl(funcs.a, funcs.n, VAR, n_group, arc_copy, sol_space_size);
        } else {
            // do simulated annealing optimization
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] run simulated annealing optimization: %ld\n",
                    __func__, sol_space_size);
#endif
            curr_cost = estimate_arc_copy_number_siman_impl(funcs.a, funcs.n, VAR, n_group, arc_copy);
        }

#ifdef DEBUG_SEG_COV_ADJUST
        for (i = 0; i < n_group; ++i)
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] arc group %u optimum copy number %d [%u %u]\n",
                    __func__, i, arc_copy[i], arc_copy_lb[i], arc_copy_ub[i]);
#endif

        if (curr_cost >= optim_cost)
            // the results is no better
            break;
        optim_cost = curr_cost;
        
        // update average sequence coverage
        double total_covs, total_lens, new_adjusted_cov, copies;
        total_covs = total_lens = 0;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del)
                continue;
            copies = 0;
            for (k = 0; k < 2; ++k) {
                v = i << 1 | k;
                na = asmg_arc_n(g, v);
                av = asmg_arc_a(g, v);
                for (j = 0; j < na; ++j) {
                    if (av[j].del)
                        continue;
                    link_id = av[j].link_id;
                    a_g = arc_group[link_id];
                    assert(a_g != (uint32_t) -1);
                    copies += arc_copy[a_g];
                }
            }
            total_lens += (double) g->vtx[i].len * copies / 2;
            if (copies)
                total_covs += (double) g->vtx[i].len * g->vtx[i].cov;
        }
    
        if (total_lens < FLT_EPSILON) {
            // no sequences included - not good
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] zero copies for all sequences\n", __func__);
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] adjusting copy number terminated in round %d of %d\n",
                    __func__, round, max_round);
#endif
            break;
        }

        new_adjusted_cov = total_covs / total_lens;
        new_adjusted_cov = MAX(new_adjusted_cov, min_avg_cov);

        if (fabs(new_adjusted_cov - adjusted_cov) < FLT_EPSILON) {
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] adjusting copy number converaged in round %d of %d\n",
                    __func__, round, max_round);
#endif
            break; // converaged
        }

        // update sequence copy number using the arc copy number information
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del)
                continue;
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] #### sequence %s ####\n", __func__, asg->seg[i].name);
#endif
        for (k = 0; k < 2; ++k) {
                v = i << 1 | k;
                na = asmg_arc_n(g, v);
                av = asmg_arc_a(g, v);
                new_copy[k] = 0;
#ifdef DEBUG_SEG_COV_ADJUST
                fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] #### %sing arcs\n", __func__, k? "incom" : "outgo");
#endif
                for (j = 0; j < na; ++j) {
                    if (av[j].del)
                        continue;
                    link_id = av[j].link_id;
                    a_g = arc_group[link_id];
                    assert(a_g != (uint32_t) -1);
                    new_copy[k] += arc_copy[a_g];
#ifdef DEBUG_SEG_COV_ADJUST
                    fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s%c -> %s%c: %d\n", __func__,
                            asg->seg[av[j].v>>1].name, "+-"[av[j].v&1],
                            asg->seg[av[j].w>>1].name, "+-"[av[j].w&1],
                            arc_copy[a_g]);
#endif
                }
            }
#ifdef DEBUG_SEG_COV_ADJUST
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] indegree: %d;  outdegree: %d; diff: %d\n", __func__,
                    new_copy[1], new_copy[0], abs(new_copy[1] - new_copy[0]));
#endif
            // only update copy number if indegree matches outdegree
            // TODO better strategy to update sequence copy number
            if (new_copy[0] == new_copy[1] && copy_number[i] != new_copy[0]) {
#ifdef DEBUG_SEG_COV_ADJUST
                fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] sequence %s copy number updated %d -> %d\n",
                        __func__, asg->seg[i].name, copy_number[i], new_copy[0]);
#endif
                copy_number[i] = new_copy[0];
                updated = 1;
            }
        }
#ifdef DEBUG_SEG_COV_ADJUST
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] sequence copy number AFTER adjusted by graph layout\n", __func__);
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] %s %lu %u %d\n", __func__,
                    asg->seg[i].name, g->vtx[i].len, g->vtx[i].cov, copy_number[i]);
        }
#endif

        // update adjusted sequence coverage
        adjusted_cov = new_adjusted_cov;
        if (_adjusted_cov)
            *_adjusted_cov = adjusted_cov;

        // update objective functions
        // more precisely the expected copy number
        for (i = 0; i < n_seg; ++i) {
            if (funcmap[i] == -1)
                continue;
            funcs.a[funcmap[i]].v_exp = g->vtx[i].cov / adjusted_cov;
            funcs.a[funcmap[i]].VAL = 0.;
        }
        // reset VARs
        for (i = 0; i < n_group; ++i) {
            var_t *var = VAR[i];
            while (var->B) var = var->next;
            VAR[i] = var;
        }

#ifdef DEBUG_SEG_COV_ADJUST
        fprintf(stderr, "[DEBUG_SEG_COV_ADJUST::%s] adjusting copy number round %d of %d DONE with adjusted coverage: %.3f\n",
                __func__, round, max_round, adjusted_cov);
#endif
    }

    for (i = 0; i < n_group; ++i) {
        var_t *var0, *var1;
        var0 = VAR[i];
        while (var0->B != 0)
            var0 = var0->next;
        var0 = var0->next;
        v = 0;
        while (!v) {
            var1 = var0->next;
            v = var0->B == 0;
            free(var0);
            var0 = var1;
        }
    }
    free(VAR);
    for (i = 0; i < funcs.n; ++i)
        free(funcs.a[i].V);
    free(funcs.a);
    free(funcmap);
    free(arc_group);
    free(arc_copy_lb);
    free(arc_copy_ub);
    free(arc_copy);

    return updated;
}

kh_u32_t *sequence_duplication_by_copy_number(asg_t *asg, int *copy_number, int allow_del)
{
    uint32_t i, n_seg;
    int copy;
    kh_u32_t *seg_dups;
    asmg_t *g;

    n_seg = asg->n_seg;
    g = asg->asmg;
    seg_dups = kh_u32_init();
    for (i = 0; i < n_seg; ++i) {
        if (g->vtx[i].del) continue;
        copy = copy_number[i];
        if (copy > 1) {
            // do copy
            make_seg_dups(asg, seg_dups, i, copy);
#ifdef DEBUG_SEG_COPY
            fprintf(stderr, "[DEBUG_SEG_COPY::%s] make %d extra cop%s of %s [%lu %u]\n", __func__, 
                    copy, copy > 1? "ies" : "y", asg->seg[i].name, g->vtx[i].len, g->vtx[i].cov);
#endif
        } else if (copy == 0 && allow_del) {
            /*** TODO is this safe? **/
            asmg_vtx_del(g, i, 1);
#ifdef DEBUG_SEG_COPY
            fprintf(stderr, "[DEBUG_SEG_COPY::%s] delete seg %s [%lu %u]\n", __func__,
                    asg->seg[i].name, g->vtx[i].len, g->vtx[i].cov);
#endif
        }
    }

#ifdef DEBUG_SEG_COPY
    fprintf(stderr, "[DEBUG_SEG_COPY::%s] graph after expansion\n", __func__);
    asg_print(asg, stderr, 1);
    fprintf(stderr, "[DEBUG_SEG_COPY::%s] expanded graph stats\n", __func__);
    asg_stat(asg, stderr);
#endif
    return seg_dups;
}

/**********************
 * BFS graph explorer *
 * ********************/

typedef struct list_node {
    uint32_t v; // sid << 31 | oriv
    struct list_node *prev;
    size_t n_m, n_n;
    struct list_node **next;
} llnode;

typedef llnode* llnodep;

static llnode *new_node(int v) {
    llnode *node;
    MYMALLOC(node, 1);
    node->v = v;
    node->n_m = 4;
    node->n_n = 0;
    node->prev = NULL;
    MYMALLOC(node->next, node->n_m);
    return node;
}

static void add_next(llnode *node, struct list_node *next)
{
    if (node->n_m == node->n_n)
        MYEXPAND(node->next, node->n_m);
    node->next[node->n_n] = next;
    node->n_n++;
}

static void llnode_destroy(llnode *node)
{
    if (!node) return;
    size_t i;
    for (i = 0; i < node->n_n; ++i)
        llnode_destroy(node->next[i]);
    free(node->next);
    free(node);
}

KDQ_INIT(llnodep)

static int path_contain_vertex(llnode *node, uint32_t v)
{
    int contained = 0;
    while (node != NULL) {
        if ((node->v >> 1) == (v >> 1)) {
            contained = 1;
            break;
        }
        node = node->prev;
    }

    return contained;
}

static llnode **graph_path_extension_bfs(asmg_t *g, llnode *root, kh_u32_t *seg_dups, int max_path, uint64_t *n_leaf, int *_exceed_limit)
{
    uint64_t i, j, v, w, nv;
    int skip, exceed_limit;
    llnode *next_node, *node;
    kdq_llnodep_t *node_q;
    void *km;
    asmg_arc_t *av;
    khint32_t k32;

    kvec_t(llnodep) leaf_node;
    kvec_t(uint32_t) dups;

    kv_init(leaf_node);
    kv_init(dups);
    km = km_init();
    node_q = kdq_init_llnodep(km);

    exceed_limit = 0;
    *kdq_pushp_llnodep(node_q) = root;

    while (kdq_size(node_q) > 0) {
        node = *kdq_shift_llnodep(node_q);
        v = node->v;
        nv = asmg_arc_n(g, v);
        av = asmg_arc_a(g, v);
        dups.n = 0;
        for (i = 0; i < nv; ++i) {
            if (av[i].del) continue;

            w = av[i].w;
            skip = 0;
            k32 = kh_u32_get(seg_dups, w >> 1);
            if (k32 < kh_end(seg_dups)) {
                // dup seg
                // check if already processed
                for (j = 0; j < dups.n; ++j) {
                    if (dups.a[j] == kh_val(seg_dups, k32)) {
                        skip = 1;
                        break;
                    }
                }
            }

            if (!skip && !path_contain_vertex(node, w)) {
                next_node = new_node(w);
                next_node->prev = node;
                add_next(node, next_node);
                *kdq_pushp_llnodep(node_q) = next_node;
                
                if (k32 < kh_end(seg_dups))
                    kv_push(uint32_t, dups, kh_val(seg_dups, k32));
            }
        }

        if (node->n_n == 0) {
            // leaf node
            kv_push(llnodep, leaf_node, node);
        }
        
        if (kdq_size(node_q) + leaf_node.n > max_path) {
            exceed_limit = 1;
            break;
        }
    }

    kdq_destroy_llnodep(node_q);
    km_destroy(km);
    kv_destroy(dups);

    if (exceed_limit) {
        *n_leaf = 0;
        *_exceed_limit = 1;
        kv_destroy(leaf_node);
        return 0;
    } else {
        *n_leaf = leaf_node.n;
        *_exceed_limit = 0;
        return leaf_node.a;
    }
}

static inline void cmp_array(uint32_t *arr, uint32_t n)
{
    if (!arr) return;
    size_t i;
    for (i = 0; i < n; i++)
        arr[i] ^= 1;
}


static inline void rev_array(uint32_t *arr, uint32_t n)
{
    if (!arr) return;
    uint32_t tmp, i, j;
    if (n == 0) return;
    for (i = 0, j = n - 1; i < j; ++i, --j) {
        tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }
}

static inline void rev_path(path_t *path)
{
    if (!path) return;
    size_t i;
    rev_array(path->v, path->nv);
    for (i = 0; i < path->nv; ++i)
        path->v[i] ^= 1;
}

static uint64_t find_source_vtx(asmg_t *g, int use_max_scc)
{
    uint64_t i, s, n_seg;
    s = UINT64_MAX;
    
    if (!use_max_scc) {
        // naively select the longest sequence
        uint64_t m_len, len;
        n_seg = g->n_vtx;
        m_len = 0;
        for (i = 0; i < n_seg; ++i) {
            if (g->vtx[i].del) continue;
            len = g->vtx[i].len * g->vtx[i].cov;
            if (m_len < len) {
                m_len = len;
                s = i;
            }
        }
    } else {
        // select the largest sequence from the largest SCC
        int c, m_c, n_scc, *scc;
        uint64_t m_len, len, *lens;
        n_seg = g->n_vtx * 2;

        MYMALLOC(scc, n_seg);
        // find strongly connected components
        n_scc = asmg_tarjans_scc(g, scc);
        
        MYCALLOC(lens, n_scc);
        // find larges SCC
        for (i = 0; i < n_seg; ++i) {
            if (scc[i] < 0) continue; // deleted vtx
            if (scc[i] != scc[i^1] || (i & 1))
                lens[scc[i]] += g->vtx[i>>1].len * g->vtx[i>>1].cov;
        }
        m_len = 0;
        m_c = -1;
        for (c = 0; c < n_scc; ++c) {
            if (m_len < lens[c]) {
                m_len = lens[c];
                m_c = c;
            }
        }
        if (m_c >= 0) {
            m_len = 0;
            for (i = 0; i < n_seg; ++i) {
                if (scc[i] != m_c) continue;
                len = g->vtx[i>>1].len * g->vtx[i>>1].cov;
                if (m_len < len) {
                    m_len = len;
                    s = i;
                }
            }
            s >>= 1;
        }
        free(scc);
        free(lens);
    }

    return s;
}

path_t *graph_path_finder_bfs(asg_t *asg, kh_u32_t *seg_dups, double min_cfrac, int64_t source, int64_t target, int max_path)
{
    uint64_t i, j, s, n_leaf;
    int circ, exceed_limit;
    asmg_t *g;
    path_t *best_cpath, *best_lpath, *best_path;
    llnode *root, *next_node, *node, **leaves;
    kvec_t(llnodep) leaf_node, root_node;
    khint32_t k32;

    g = asg->asmg;
    best_path = NULL;
    
    // find source vertex
    if (source >= 0) {
        // find a copy of source node
        s = UINT64_MAX;
        for (k32 = (khint32_t) 0; k32 < kh_end(seg_dups); k32++) {
            if (kh_exist(seg_dups, k32) && 
                kh_val(seg_dups, k32) == (source >> 1) &&
                !g->vtx[s=kh_key(seg_dups, k32)].del)
                break;
        }
        if (s != UINT64_MAX)
            s = s << 1 | (source & 1);
        else
            s = source;
    } else {
        s = find_source_vtx(g, 1);
        if (s == UINT64_MAX) return NULL;
        s <<= 1;
        if (target >= 0)
            fprintf(stderr, "[W::%s] cannot specify target node without source node\n", __func__);
    }

    kv_init(root_node);
    kv_init(leaf_node);

    root = new_node(s);
    kv_push(llnodep, root_node, root);

    n_leaf = 0;
    exceed_limit = 0;
    leaves = graph_path_extension_bfs(g, root, seg_dups, max_path, &n_leaf, &exceed_limit);
    
    // for linear paths do extension from the other direction of root node
    // no - should do this even if the path is circular
    // n_leaf will be zero if path number exceeds limit
    for (i = 0; i < n_leaf; ++i) {
        node = leaves[i];
        // circ = asmg_arc1(g, node->v, s << 1) != 0;
        // if (circ) {
        //     kv_push(llnodep, leaf_node, node);
        //     continue;
        // }
        
        // make a new path tracing back from leaf to root
        root = new_node(node->v^1);
        kv_push(llnodep, root_node, root);
        while (node->prev != NULL) {
            next_node = new_node(node->prev->v^1);
            next_node->prev = root;
            add_next(root, next_node);
            root = next_node;
            node = node->prev;
        }

        assert(root->v == (s^1)); // should be always ended with seq s
        
        if (source < 0) {
            // only need to do this if source has not been specified
            uint64_t n = 0;
            exceed_limit = 0;
            llnode **tmp_nodes = graph_path_extension_bfs(g, root, seg_dups, max_path, &n, &exceed_limit);
            
            for (j = 0; j < n; ++j)
                kv_push(llnodep, leaf_node, tmp_nodes[j]);
            free(tmp_nodes);

            if (exceed_limit || leaf_node.n > max_path) {
                exceed_limit = 1;
                break;
            }
        } else
            kv_push(llnodep, leaf_node, root);
    }
    free(leaves);

    if (exceed_limit) {
        fprintf(stderr, "[W::%s] path exploration exceeds limit %d\n", __func__, max_path);
        fprintf(stderr, "[W::%s] consider an larger value of '-N'\n", __func__);
        goto final_clean;
    }

#ifdef DEBUG_PATH_FINDER
    fprintf(stderr, "[DEBUG_PATH_FINDER::%s] number leaf nodes: %lu\n", __func__, leaf_node.n);
#endif

    kvec_t(uint32_t) path;
    kv_init(path);
    best_cpath = best_lpath = NULL;
    for (i = 0; i < leaf_node.n; ++i) {
        node = leaf_node.a[i];
        path.n = 0;
        while (node != NULL) {
            kv_push(uint32_t, path, node->v);
            node = node->prev;
        }
        //rev_array(path.a, path.n);
        cmp_array(path.a, path.n);

        if (target >= 0) {
            // clip at the last target node
            int64_t k = path.n - 1;
            while (k >= 0) {
                if (target == path.a[k])
                    break;
                k32 = kh_u32_get(seg_dups, path.a[k]>>1);
                if (k32 < kh_end(seg_dups) &&
                    target == (kh_val(seg_dups, k32)<<1 | (path.a[k]&1)))
                    break;
                k -= 1;
            }
            if (k >= 0)
                path.n = k + 1;
            else 
                continue;
        }

        circ = asmg_arc1(g, path.a[path.n-1], path.a[0]) != 0;

#ifdef DEBUG_PATH_FINDER
        uint64_t n = path.n;
        fprintf(stderr, "[DEBUG_PATH_FINDER::%s] Path %lu [%s] (%lu): %s%c", __func__, i, 
                circ? "circle" : "linear", n, asg->seg[path.a[0]>>1].name, "+-"[path.a[0]&1]);
        for (j = 1; j < n; ++j)
            fprintf(stderr, ",%s%c", asg->seg[path.a[j]>>1].name, "+-"[path.a[j]&1]);
        fprintf(stderr, "\n");
#endif

        // calculate sequence length
        asmg_arc_t *a;
        uint64_t l, l1;
        uint32_t ec;
        double cov, wl;
        
        l = g->vtx[path.a[0]>>1].len;
        cov = g->vtx[path.a[0]>>1].cov;
        wl = cov * l;
        ec = 0;
        for (j = 1; j < path.n; ++j) {
            a = asmg_arc1(g, path.a[j-1], path.a[j]);
            assert(!!a);
            l1 = g->vtx[path.a[j]>>1].len - a->ls;
            cov = g->vtx[path.a[j]>>1].cov;
            l += l1;
            wl += cov * l1;
            ec += a->cov;
        }
        
        if (circ) {
            a = asmg_arc1(g, path.a[path.n-1], path.a[0]);
            assert(!!a);
            l1 = a->ls;
            cov = g->vtx[path.a[0]>>1].cov;
            l -= l1;
            wl -= cov * l1;
            ec += a->cov;
        }

        int replace = 0;
        if (circ) {
            if (best_cpath == NULL)
                replace = 1;
            else
                if (best_cpath->wlen != wl)
                    replace = best_cpath->wlen < wl;
                else if (best_cpath->len != l)
                    replace = best_cpath->len < l;
                else if (best_cpath->ec != ec)
                    replace = best_cpath->ec < ec;
                else replace = best_cpath->nv < path.n;
        } else {
            if (best_lpath == NULL)
                replace = 1;
            else
                if (best_lpath->wlen != wl)
                    replace = best_lpath->wlen < wl;
                else if (best_lpath->len != l)
                    replace = best_lpath->len < l;
                else if (best_lpath->ec != ec)
                    replace = best_lpath->ec < ec;
                else replace = best_lpath->nv < path.n;
        }

        if (replace) {
            // replace copied vertices with originals
            for (j = 0; j < path.n; ++j) {
                k32 = kh_u32_get(seg_dups, path.a[j]>>1);
                if (k32 < kh_end(seg_dups))
                    path.a[j] = kh_val(seg_dups, k32)<<1 | (path.a[j]&1);
            }

            best_path = circ? best_cpath : best_lpath;
            if (best_path == NULL) {
                MYCALLOC(best_path, 1);
                if (circ) best_cpath = best_path;
                else best_lpath = best_path;
            }
            if (best_path->nv < path.n)
                MYREALLOC(best_path->v, path.n);
            memcpy(best_path->v, path.a, sizeof(uint32_t) * path.n);
            best_path->circ = circ;
            best_path->nv = path.n;
            best_path->len = l;
            best_path->wlen = wl;
            best_path->ec = ec;
        }
    }

    if (best_cpath == NULL)
        best_path = best_lpath;
    else if (best_lpath == NULL)
        best_path = best_cpath;
    else {
        if (best_cpath->wlen + DBL_EPSILON >= best_lpath->wlen * min_cfrac) {
            best_path = best_cpath;
            path_destroy(best_lpath);
            free(best_lpath);
        } else {
            best_path = best_lpath;
            path_destroy(best_cpath);
            free(best_cpath);
        }
    }

    kv_destroy(path);

final_clean:
    for (i = 0; i < root_node.n; ++i)
        llnode_destroy(root_node.a[i]);

    kv_destroy(root_node);
    kv_destroy(leaf_node);

    return best_path;
}
/******END BFS ********/


/**********************
 * DFS graph explorer *
 * ********************/

typedef struct { size_t n; double l; uint32_t *a; } u32_v;
typedef struct { size_t n, m; u32_v *a; } u32_vv;

static int PSORT(const void *a, const void *b)
{
    double x, y;
    x = ((u32_v *) a)->l;
    y = ((u32_v *) b)->l;
    return ((x < y) - (x > y));
}

static int dfs_core(asmg_t *g, uint32_t *vtex, uint32_t vtxn, int *visited, int *copy_number, u32_vv *paths, int64_t target, int max_path)
{
    if (paths->n >= max_path) return 1;

    uint32_t i, s, nv;
    uint64_t v, w;
    asmg_arc_t *av;
    u32_v *path;

    v = vtex[vtxn-1];
    s = 0;
    // if v is target and has been visited completely
    // then stop extension
    if (v != target || visited[v>>1] < copy_number[target>>1]) {
        nv = asmg_arc_n(g, v);
        av = asmg_arc_a(g, v);
        for (i = 0; i < nv; ++i) {
            if (av[i].del) continue;
            
            w = av[i].w;
            if (visited[w>>1] >= copy_number[w>>1])
                continue; // already visited
            
            vtex[vtxn] = w;
            vtxn++;
            visited[w>>1]++;
            dfs_core(g, vtex, vtxn, visited, copy_number, paths, target, max_path);
            visited[w>>1]--;
            vtxn--;
            s++;
        }
    }

    if (!s) {
        // no extension
        // add to path
        kv_pushp(u32_v, *paths, &path);
        path->n = vtxn;
        MYMALLOC(path->a, vtxn);
        memcpy(path->a, vtex, sizeof(uint32_t) * vtxn);
    }

    return (paths->n >= max_path);
}

static int graph_path_extension_dfs(asmg_t *g, uint32_t *vs, uint32_t vn, uint32_t *vtex, int *visited, int *copy_number, u32_vv *paths, int64_t target, int max_path)
{
    uint32_t i;

    MYBZERO(visited, g->n_vtx);
    for (i = 0; i < vn; i++)
        visited[vs[i]>>1] += 1;
    memcpy(vtex, vs, sizeof(uint32_t) * vn);

    return dfs_core(g, vtex, vn, visited, copy_number, paths, target, max_path);
}

path_t *graph_path_finder_dfs(asg_t *asg, int *copy_number, double min_cfrac, int64_t source, int64_t target, int max_path)
{
    uint32_t i, j, vn, s, *vs, *vtex;
    int *visited, max_vtxn, circ, exceed_limit;
    asmg_t *g;
    path_t *best_cpath, *best_lpath, *best_path;
    u32_vv paths;

    g = asg->asmg;
    best_path = NULL;
    
    // find source vertex
    if (source >= 0)
        s = source;
    else {
        s = find_source_vtx(g, 1);
        if (s >= g->n_vtx) return NULL;
        s <<= 1;
        if (target >= 0)
            fprintf(stderr, "[W::%s] cannot specify target node without source node\n", __func__);
    }

    max_vtxn = 0;
    for (i = 0; i < asg->n_seg; i++)
        max_vtxn += copy_number[i];

    kv_init(paths);
    kv_resize(u32_v, paths, 1024);
    MYMALLOC(visited, asg->n_seg);
    MYMALLOC(vtex, max_vtxn);
    exceed_limit = graph_path_extension_dfs(g, &s, 1, vtex, visited, copy_number, &paths, target, max_path);

    if (source < 0 && !exceed_limit) {
        // need to extend the path from the other direction
        double wl;
        asmg_arc_t *a;
        uint32_t n = paths.n;
        int64_t t = target;
        if (t >= 0) t ^= 1;
        for (i = 0; i < n; i++) {
            vs = paths.a[i].a;
            vn = paths.a[i].n;
            wl = (double) g->vtx[vs[0]>>1].cov * g->vtx[vs[0]>>1].len;
            for (j = 1; j < vn; j++) {
                a = asmg_arc1(g, vs[j-1], vs[j]);
                //assert(!!a);
                wl += (double) g->vtx[vs[j]>>1].cov * (g->vtx[vs[j]>>1].len - a->ls);
            }
            paths.a[i].l = wl;
        }
        // sort paths by weighted length
        qsort(paths.a, paths.n, sizeof(u32_v), PSORT);

        for (i = 0; i < n; i++) {
            vs = paths.a[i].a;
            vn = paths.a[i].n;

            rev_array(vs, vn);
            cmp_array(vs, vn);

            exceed_limit = graph_path_extension_dfs(g, vs, vn, vtex, visited, copy_number, &paths, t, max_path<<1);
            if (exceed_limit) break;

            // this is not needed any more
            free(vs);
            paths.a[i].a = NULL;
        }

        for (i = 0; i < paths.n; i++) {
            vs = paths.a[i].a;
            vn = paths.a[i].n;

            rev_array(vs, vn);
            cmp_array(vs, vn);
        }
    }

    free(vtex);
    free(visited);

    if (exceed_limit) {
        fprintf(stderr, "[W::%s] path exploration exceeds limit %d\n", __func__, max_path);
        fprintf(stderr, "[W::%s] the results might be suboptimal\n", __func__);
        fprintf(stderr, "[W::%s] consider a larger value of '-N'\n", __func__);
    }

    double cov, wl;
    asmg_arc_t *a;
    uint64_t l, l1;
    uint32_t ec;
    best_cpath = best_lpath = NULL;
    for (i = 0; i < paths.n; ++i) {
        vs = paths.a[i].a;
        vn = paths.a[i].n;

        if (!vs) continue;

        if (target >= 0) {
            // clip at the last target node
            int k = vn - 1;
            while (k >= 0) {
                if (target == vs[k])
                    break;
                k -= 1;
            }
            if (k >= 0)
                vn = k + 1;
            else 
                continue;
        }

        circ = asmg_arc1(g, vs[vn-1], vs[0]) != 0;

#ifdef DEBUG_PATH_FINDER
        fprintf(stderr, "[DEBUG_PATH_FINDER::%s] Path %lu [%s] (%lu): %s%c", __func__, i, 
                circ? "circle" : "linear", n, asg->seg[vs[0]>>1].name, "+-"[vs[0]&1]);
        for (j = 1; j < vn; ++j)
            fprintf(stderr, ",%s%c", asg->seg[vs[j]>>1].name, "+-"[vs[j]&1]);
        fprintf(stderr, "\n");
#endif

        // calculate sequence length
        l = g->vtx[vs[0]>>1].len;
        cov = g->vtx[vs[0]>>1].cov;
        wl = cov * l;
        ec = 0;
        for (j = 1; j < vn; ++j) {
            a = asmg_arc1(g, vs[j-1], vs[j]);
            //assert(!!a);
            l1 = g->vtx[vs[j]>>1].len - a->ls;
            cov = g->vtx[vs[j]>>1].cov;
            l += l1;
            wl += cov * l1;
            ec += a->cov;
        }
        
        if (circ) {
            a = asmg_arc1(g, vs[vn-1], vs[0]);
            //assert(!!a);
            l1 = a->ls;
            cov = g->vtx[vs[0]>>1].cov;
            l -= l1;
            wl -= cov * l1;
            ec += a->cov;
        }

        int replace = 0;
        if (circ) {
            if (best_cpath == NULL)
                replace = 1;
            else
                if (best_cpath->wlen != wl)
                    replace = best_cpath->wlen < wl;
                else if (best_cpath->len != l)
                    replace = best_cpath->len < l;
                else if (best_cpath->ec != ec)
                    replace = best_cpath->ec < ec;
                else replace = best_cpath->nv < vn;
        } else {
            if (best_lpath == NULL)
                replace = 1;
            else
                if (best_lpath->wlen != wl)
                    replace = best_lpath->wlen < wl;
                else if (best_lpath->len != l)
                    replace = best_lpath->len < l;
                else if (best_lpath->ec != ec)
                    replace = best_lpath->ec < ec;
                else replace = best_lpath->nv < vn;
        }

        if (replace) {
            best_path = circ? best_cpath : best_lpath;
            if (best_path == NULL) {
                MYCALLOC(best_path, 1);
                if (circ) best_cpath = best_path;
                else best_lpath = best_path;
            }
            free(best_path->v);
            best_path->v = vs;
            best_path->circ = circ;
            best_path->nv = vn;
            best_path->len = l;
            best_path->wlen = wl;
            best_path->ec = ec;
            paths.a[i].a = NULL; // vs
        }
    }

    if (best_cpath == NULL)
        best_path = best_lpath;
    else if (best_lpath == NULL)
        best_path = best_cpath;
    else {
        if (best_cpath->wlen + DBL_EPSILON >= best_lpath->wlen * min_cfrac) {
            best_path = best_cpath;
            path_destroy(best_lpath);
            free(best_lpath);
        } else {
            best_path = best_lpath;
            path_destroy(best_cpath);
            free(best_cpath);
        }
    }

    for (i = 0; i < paths.n; i++)
        free(paths.a[i].a);
    kv_destroy(paths);
    
    return best_path;
}

/******END DFS ********/

void path_report(asg_t *asg, path_t *path)
{
    if (path == NULL || path->nv == 0)
        return;

    uint32_t i, vn, ec, circ, *vs;
    uint64_t l, l1;
    double cov, wl;
    asmg_t *g;
    asmg_arc_t *a;

    g = asg->asmg;
    vs = path->v;
    vn = path->nv;

    circ = asmg_arc1(g, vs[vn-1], vs[0]) != 0;

    // calculate sequence length
    l = g->vtx[vs[0]>>1].len;
    cov = g->vtx[vs[0]>>1].cov;
    wl = cov * l;
    ec = 0;
    for (i = 1; i < vn; ++i) {
        a = asmg_arc1(g, vs[i-1], vs[i]);
        //assert(!!a);
        l1 = g->vtx[vs[i]>>1].len - a->ls;
        cov = g->vtx[vs[i]>>1].cov;
        l += l1;
        wl += cov * l1;
        ec += a->cov;
    }
    
    if (circ) {
        a = asmg_arc1(g, vs[vn-1], vs[0]);
        //assert(!!a);
        l1 = a->ls;
        cov = g->vtx[vs[0]>>1].cov;
        l -= l1;
        wl -= cov * l1;
        ec += a->cov;
    }

    path->circ = circ;
    path->len = l;
    path->wlen = wl;
    path->ec = ec;
}

static int pcmpfunc(const void *a, const void *b)
{
    path_t x, y;
    x = *(path_t *) a;
    y = *(path_t *) b;
    if (x.wlen != y.wlen)
        return (x.wlen < y.wlen) - (x.wlen > y.wlen);
    if (x.len != y.len)
        return (x.len < y.len) - (x.len > y.len);
    if (x.ec != y.ec)
        return (x.ec < y.ec) - (x.ec > y.ec);
    if (x.nv != y.nv)
        return (x.nv < y.nv) - (x.nv > y.nv);

    return 0;
}

static char *strdup1(char *src, uint32_t n)
{
    char *dst;
    MYMALLOC(dst, n+1);
    memcpy(dst, src, n);
    dst[n] = '\0';
    return dst;
}

path_t make_path_from_str(asg_t *asg, char *path_str, char *sid)
{
    int circ;
    uint32_t v, cov;
    uint64_t len, len1;
    double wlen;
    kvec_t(uint32_t) vt;
    char *ptr, *s;

    kv_init(vt);

    while (isspace(*path_str) && *path_str != '\0')
        ++path_str;

    while (*path_str != '\0') {
        ptr = path_str;
        while (!isspace(*ptr) && *ptr != ',' && *ptr !='\0')
            ++ptr;

        if (*(ptr-1) == '+' || *(ptr-1) == '-') {
            s = strdup1(path_str, ptr - path_str - 1);
            v = asg_name2id(asg, s);
            if (v == UINT32_MAX) {
                fprintf(stderr, "[E::%s] sequence does not exist: %s\n", __func__, s);
                exit(EXIT_FAILURE);
            }
            v = (v<<1) | (*(ptr-1)=='-');
            kv_push(uint32_t, vt , v);
            free(s);
        } else {
            fprintf(stderr, "[E::%s] invalid path string: %s\n", __func__, path_str);
            exit(EXIT_FAILURE);
        }

        if (isspace(*ptr) || *ptr == '\0') break;
        path_str = ptr + 1;
    }
    
    if (vt.n == 0) {
        fprintf(stderr, "[E::%s] invalid path string: %s\n", __func__, path_str);
        exit(EXIT_FAILURE);
    }

    size_t i;
    asmg_t *g;
    asmg_arc_t *a;
    uint32_t ec;
    g = asg->asmg;
    a = asmg_arc1(g, vt.a[vt.n-1], vt.a[0]);
    circ = !!a;
    len = g->vtx[vt.a[0]>>1].len;
    cov = g->vtx[vt.a[0]>>1].cov;
    wlen = (double) cov * len;
    ec = 0;
    if (circ) len -= a->ls, wlen -= cov * a->ls, ec += a->cov;
    for (i = 1; i < vt.n; ++i) {
        len1 = g->vtx[vt.a[i]>>1].len;
        cov = g->vtx[vt.a[i]>>1].cov;
        len += len1;
        wlen += (double) cov * len1;
        a = asmg_arc1(g, vt.a[i-1], vt.a[i]);
        if (!a) {
            fprintf(stderr, "[W::%s] gap introduced as link does not exist: %s%c -> %s%c\n", __func__,
                    asg->seg[vt.a[i-1]>>1].name, "+-"[vt.a[i-1]&1],
                    asg->seg[vt.a[i]>>1].name, "+-"[vt.a[i]&1]);
        } else {
            len -= a->ls;
            wlen -= (double) cov * a->ls;
            ec += a->cov;
        }
    }

    path_t path = {sid? strdup1(sid, strlen(sid)) : 0, vt.n, circ, 0, vt.a, len, wlen, ec};
    
    return path;
}

void path_sort(path_v *paths)
{
    // sort by wlen -> len -> circ -> nv
    qsort(paths->a, paths->n, sizeof(path_t), pcmpfunc);

    size_t i;
    double b_ll, b_cl;
    b_ll = b_cl = .0;
    for (i = 0 ; i < paths->n; ++i) {
        if (!paths->a[i].circ && paths->a[i].wlen > b_ll)
            b_ll = paths->a[i].wlen;
        if (paths->a[i].circ && paths->a[i].wlen > b_cl)
            b_cl = paths->a[i].wlen;
    }
    if (b_cl >= b_ll)
        b_ll = DBL_MAX;
    // find best
    for (i = 0 ; i < paths->n; ++i) {
        if (!paths->a[i].circ && paths->a[i].wlen >= b_ll)
            paths->a[i].best = 1;
        if (paths->a[i].circ && paths->a[i].wlen >= b_cl)
            paths->a[i].best = 1;
    }
}

void path_stats(asg_t *asg, path_v *paths, FILE *fo)
{
    uint32_t i, j, n;
    int nd_i, nd_n, nd_l, nd_w, nd_e;
    uint64_t max_l;
    uint32_t max_e;
    double max_wl;
    path_t path;

    nd_n = 0;
    max_l = 0;
    max_e = 0;
    max_wl = .0;
    for (i = 0; i < paths->n; ++i) {
        if (paths->a[i].nv > nd_n)
            nd_n = paths->a[i].nv;
        if (paths->a[i].len > max_l)
            max_l = paths->a[i].len;
        if (paths->a[i].ec > max_e)
            max_e = paths->a[i].ec;
        if (paths->a[i].wlen > max_wl)
            max_wl = paths->a[i].wlen;
    }
    nd_i = floor(log10(fabs(paths->n))) + 1;
    nd_n = floor(log10(fabs(nd_n))) + 1;
    nd_l = floor(log10(fabs(max_l))) + 1;
    nd_e = floor(log10(fabs(max_e))) + 1;
    nd_w = floor(log10(fabs(max_wl))) + 3;

    for (i = 0; i < paths->n; ++i) {
        path = paths->a[i];
        n = path.nv;
        fprintf(fo, "%c %-*u %s %-*u %-*u %-*.1f %-*u %s%c", path.best? '*' : '#', nd_i, i, path.circ? "circle" : "linear",
                nd_n, n, nd_l, path.len, nd_w, path.wlen, nd_e, path.ec, asg->seg[path.v[0]>>1].name, "+-"[path.v[0]&1]);
        for (j = 1; j < n; ++j)
            fprintf(fo, ",%s%c", asg->seg[path.v[j]>>1].name, "+-"[path.v[j]&1]);
        fprintf(fo, "\n");
    }
}

static char comp_table[] = {
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
    48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', 123, 124, 125, 126, 127
};

static void put_chars(char *seq, int len, int rv, int ow, FILE *fo, uint64_t *l, int line_wd)
{
    int i;
    if (!rv) {
        // forward
        for (i = ow; i < len; ++i) {
            fputc(seq[i], fo);
            if (++(*l) % line_wd == 0)
                fputc('\n', fo);
        }
    } else {
         // reverse
        for (i = len - ow - 1; i >= 0; --i) {
            fputc(comp_table[(int) seq[i]], fo);
            if (++(*l) % line_wd == 0)
                fputc('\n', fo);
        }
    }
}

static void make_gap(FILE *fo, uint64_t *l, int line_wd, int gap_size)
{
    int i;
    for (i = 0; i < gap_size; ++i) {
        fputc('N', fo);
        if (++(*l) % line_wd == 0)
            fputc('\n', fo);
    }
}

void print_seq(asg_t *asg, path_t *path, FILE *fo, int id, int force_linear, int line_wd, int gap_size)
{
    uint32_t i, n, v, lo, cov, n_gap;
    uint64_t l;
    asmg_t *g;
    asmg_arc_t *a;
    
    n = path->nv;
    if (n == 0) return;

    for (i = 0; i < n; ++i) {
        v = path->v[i];
        if (!asg->seg[v>>1].seq) {
            fprintf(stderr, "[E::%s] cannot make FASTA output: sequence not included in the GFA file\n", __func__);
            return;
        }
    }

    g = asg->asmg;
    lo = 0;
    cov = 0;
    if (path->circ && force_linear) {
        a = asmg_arc1(g, path->v[n-1], path->v[0]);
        assert(!!a);
        lo = a->ls;
        cov = g->vtx[path->v[0]>>1].cov;
    }
    
    if (path->sid)
        fprintf(fo, ">%s\tlength=%u wlength=%.1f nv=%u ec=%u circular=%s path=%s%c", path->sid, path->len + lo, path->wlen + (double) cov * lo,
                path->nv, path->ec, (force_linear || !path->circ)? "false" : "true", asg->seg[path->v[0]>>1].name, "+-"[path->v[0]&1]);
    else
        fprintf(fo, ">ctg%06d%c\tlength=%u wlength=%.1f nv=%u ec=%u circular=%s path=%s%c", id, (force_linear || !path->circ)? 'l' : 'c', 
                path->len + lo, path->wlen + (double) cov * lo, path->nv, path->ec, (force_linear || !path->circ)? "false" : "true", 
                asg->seg[path->v[0]>>1].name, "+-"[path->v[0]&1]);

    for (i = 1; i < n; ++i)
        fprintf(fo, ",%s%c", asg->seg[path->v[i]>>1].name, "+-"[path->v[i]&1]);
    fprintf(fo, "\n");

    l = 0;
    v = path->v[0];
    lo = (force_linear || !path->circ)? 0 : asmg_arc1(g, path->v[n-1], v)->ls;
    put_chars(asg->seg[v>>1].seq, asg->seg[v>>1].len, v&1, lo, fo, &l, line_wd);

    n_gap = 0;
    for (i = 1; i < n; ++i) {
        v = path->v[i];
        a = asmg_arc1(g, path->v[i-1], v);
        if (!!a) {
            put_chars(asg->seg[v>>1].seq, asg->seg[v>>1].len, v&1, a->ls, fo, &l, line_wd);
        } else {
            make_gap(fo, &l, line_wd, gap_size);
            put_chars(asg->seg[v>>1].seq, asg->seg[v>>1].len, v&1, 0, fo, &l, line_wd);
            ++n_gap;
        }
    }

    if (!path->circ || !force_linear)
        assert(l - (uint64_t) n_gap * gap_size == path->len);
    if (l % line_wd != 0)
        fputc('\n', fo);
}

void print_all_best_seqs(asg_t *g, path_v *paths, FILE *fo)
{
    if (paths->n == 0)
        return;

    size_t i;
    for (i = 0; i < paths->n; ++i)
        if (!!paths->a[i].best)
            print_seq(g, &paths->a[i], fo, i, 0, 60, 100);
}



/******************
 * Assembly Graph *
 * ****************/

KSTREAM_INIT(gzFile, gzread, 65536)
KHASHL_MAP_INIT(KH_LOCAL, kh_sdict_t, kh_sdict, kh_cstr_t, uint32_t, kh_hash_str, kh_eq_str)
typedef kh_sdict_t sdhash_t;

asg_t *asg_init()
{
    asg_t *g;
    MYCALLOC(g, 1);
    g->h_seg = kh_sdict_init();
    MYCALLOC(g->asmg, 1);
    return g;
}

static void asg_seg_destroy(asg_seg_t *seg)
{
    if (seg->name) free(seg->name);
    if (seg->seq) free(seg->seq);
}

void asg_destroy(asg_t *g)
{
    if (!g) return;
    uint64_t i;
    for (i = 0; i < g->n_seg; ++i)
        asg_seg_destroy(&g->seg[i]);
    free(g->seg);
    // key has been freed by asg_seg_destroy
    if (g->h_seg) kh_sdict_destroy(g->h_seg);
    if (g->asmg) asmg_destroy(g->asmg);
    free(g);
}

uint32_t asg_add_seg(asg_t *g, char *name, int allow_dups)
{
    if (!name)
        return UINT32_MAX;
    sdhash_t *h = g->h_seg;
    khint_t k;
    int absent;
    k = kh_sdict_put(h, name, &absent);
    if (absent) {
        asg_seg_t *s;
        if (g->n_seg == g->m_seg)
            MYEXPAND(g->seg, g->m_seg);
        s = &g->seg[g->n_seg];
        s->len = 0;
        s->seq = 0;
        s->cov = 0;
        kh_key(h, k) = s->name = strdup(name);
        kh_val(h, k) = g->n_seg++;
    } else if (!allow_dups) {
        fprintf(stderr, "[E::%s] duplicate segment '%s'\n", __func__, name);
        exit(EXIT_FAILURE);
    }
    return kh_val(h, k);
}

uint32_t asg_add_seg1(asg_t *g, char *name, char *seq, uint32_t len, uint64_t cov, int allow_dups)
{
    uint32_t k = asg_add_seg(g, name, allow_dups);
    asg_seg_t *s = &g->seg[k];
    s->seq = strdup(seq);
    s->len = len;
    s->cov = cov;
    return k;
}

uint32_t asg_name2id(asg_t *g, char *name)
{
    sdhash_t *h = g->h_seg;
    khint_t k;
    k = kh_sdict_get(h, name);
    return k == kh_end(h)? UINT32_MAX : kh_val(h, k);
}

static void asg_update_seg_seq(asg_seg_t *seg, uint32_t l, char *s)
{
    if (seg->seq) free(seg->seq);
    MYMALLOC(seg->seq, l+1);
    memcpy(seg->seq, s, l);
    seg->seq[l] = 0;
    seg->len = l;
    seg->cov = 0;
}

asmg_t *asg_make_asmg_copy(asmg_t *g, asmg_t *_g)
{
    asmg_t *g1;
    
    if (_g)
        g1 = _g;
    else
        MYCALLOC(g1, 1);
    
    g1->n_vtx = g1->m_vtx = g->n_vtx;
    MYCALLOC(g1->vtx, g1->n_vtx);
    // this works as a and seq are always 0
    memcpy(g1->vtx, g->vtx, sizeof(asmg_vtx_t) * g1->n_vtx);
    g1->n_arc = g1->m_arc = g->n_arc;
    MYCALLOC(g1->arc, g1->n_arc);
    memcpy(g1->arc, g->arc, sizeof(asmg_arc_t) * g1->n_arc);

    uint64_t n_vtx = asmg_vtx_n(g1);
    MYCALLOC(g1->idx_p, n_vtx);
    memcpy(g1->idx_p, g->idx_p, sizeof(uint64_t) * n_vtx);
    MYCALLOC(g1->idx_n, n_vtx);
    memcpy(g1->idx_n, g->idx_n, sizeof(uint64_t) * n_vtx);
    
    return g1;
}

asg_t *asg_make_copy(asg_t *asg)
{
    uint64_t i;
    khint_t k;
    int absent;

    // seg sequence not copied
    asg_t *asg1 = asg_init();
    asg_seg_t *seg, *seg1;
    asg1->m_seg = asg1->n_seg = asg->n_seg;
    MYCALLOC(asg1->seg, asg1->n_seg);
    kh_sdict_t *h_seg = (kh_sdict_t *) asg1->h_seg;
    for (i = 0; i < asg1->n_seg; ++i) {
        seg = &asg->seg[i];
        seg1 = &asg1->seg[i];
        seg1->name = strdup(seg->name);
        seg1->len = seg->len;
        seg1->cov = seg->cov;
        k = kh_sdict_put(h_seg, seg1->name, &absent);
        kh_val(h_seg, k) = i;
    }
    asg_make_asmg_copy(asg->asmg, asg1->asmg);

    return asg1;
}

int clean_graph_by_sequence_coverage(asg_t *asg, double min_cf, int max_copy, int verbose)
{
    uint32_t i, j, n_seg, nv;
    uint8_t *visited;
    double avg_cov;
    kvec_t(uint32_t) rm_v;
    asmg_t *g;

    g = asg->asmg;
    n_seg = asg->n_seg;
    kv_init(rm_v);

    MYCALLOC(visited, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (visited[i] || g->vtx[i].del) continue;
        asg->asmg = asg_make_asmg_copy(g, 0);
        asmg_subgraph(asg->asmg, &i, 1, 0, 0, 0, 1);
        avg_cov = graph_sequence_coverage_precise(asg, min_cf, 0, max_copy, 0,0,0,0, 0);
        if (verbose > 1)
            fprintf(stderr, "[M::%s] subgraph seeding from %s per copy average coverage: %.3f\n", 
                    __func__, asg->seg[i].name, avg_cov);
        for (j = 0; j < n_seg; ++j) {
            if (asg->asmg->vtx[j].del) continue;
            if (avg_cov && asg->asmg->vtx[j].cov / avg_cov < min_cf)
                kv_push(uint32_t, rm_v, j);
            visited[j] = 1;
        }
        asmg_destroy(asg->asmg);
    }
    asg->asmg = g; // roll back

    nv = rm_v.n;
    for (i = 0; i < nv; ++i)
        asmg_vtx_del(g, rm_v.a[i], 1);
    asmg_finalize(g, 0);

    if (verbose > 1) {
        if (verbose > 2) {
            fprintf(stderr, "[M::%s] graph after cleaning\n", __func__);
            asg_print(asg, stderr, 1);
        }
        fprintf(stderr, "[M::%s] number sequence cleaned: %u\n", __func__, nv);
        for (i = 0; i < nv; ++i)
            fprintf(stderr, "[M::%s] %s removed\n", __func__, asg->seg[rm_v.a[i]].name);
        fprintf(stderr, "[M::%s] graph stats after cleaning\n", __func__);
        asg_stat(asg, stderr);
    }
    
    kv_destroy(rm_v);
    free(visited);

    return nv;
}

double sequence_covered_by_path(asg_t *asg, path_t *path, uint32_t len)
{
    uint32_t i, n, l, *v;
    uint8_t *flag;
    MYCALLOC(flag, asg->n_seg);
    v = path->v;
    l = 0;
    for (i = 0, n = path->nv; i < n; ++i) {
        if (!flag[v[i]>>1]) {
            l += asg->seg[v[i]>>1].len;
            flag[v[i]>>1] = 1;
        }
    }
    free(flag);
    return (double) l / len;
}

/****************
 * GFA File I/O *
 * **************/
// adapted from https://github.com/lh3/gfatools

#define aux_decimal_val(type, p) *((type *) (p))
// always as double
static double gfa_aux_decimal_value(uint8_t *aux)
{
    double val = .0;
    if (*aux == 'c') val = aux_decimal_val(int8_t, aux + 1);
    else if (*aux == 'C') val = aux_decimal_val(uint8_t, aux + 1);
    else if (*aux == 's') val = aux_decimal_val(int16_t, aux + 1);
    else if (*aux == 'S') val = aux_decimal_val(uint16_t, aux + 1);
    else if (*aux == 'i') val = aux_decimal_val(int64_t, aux + 1);
    else if (*aux == 'I') val = aux_decimal_val(uint64_t, aux + 1);
    else if (*aux == 'f') val = aux_decimal_val(double, aux + 1);
    return val;
}

static inline int gfa_aux_type2size(int x)
{
    if (x == 'C' || x == 'c' || x == 'A') return 1;
    else if (x == 'S' || x == 's') return 2;
    else if (x == 'I' || x == 'i' || x == 'f') return 8;
    else return 0;
}

#define __skip_tag(s) do { \
    int type = *(s); \
    ++(s); \
    if (type == 'Z') { while (*(s)) ++(s); ++(s); } \
    else if (type == 'B') (s) += 5 + gfa_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
    else (s) += gfa_aux_type2size(type); \
} while(0)

static uint8_t *gfa_aux_get(int l_data, const uint8_t *data, const char tag[2])
{
    const uint8_t *s = data;
    int y = tag[0]<<8 | tag[1];
    while (s < data + l_data) {
        int x = (int)s[0]<<8 | s[1];
        s += 2;
        if (x == y) return (uint8_t*)s;
        __skip_tag(s);
    }
    return 0;
}

// GFA tag
char TAG_ARC_COV[4]; // arc coverage EC:i
char TAG_SEQ_COV[4]; // seq coverage SC:f
char TAG_SBP_COV[4]; // seq total base coverage KC:i FC:i

int is_valid_gfa_tag(const char *tag)
{
    int is_valid = 0;
    if (strlen(tag) != 4)
        return is_valid;
    is_valid = isalpha(tag[0]) && 
        (isalpha(tag[1]) || isdigit(tag[1])) && 
        (tag[2] == ':') &&
        (tag[3] == 'A' || tag[3] == 'i' || tag[3] == 'f' || tag[3] == 'Z' || tag[3] == 'B');
    return is_valid;
}

static int gfa_aux_parse(char *s, uint8_t **data, int *max)
{
    char *q, *p;
    kstring_t str;
    if (s == 0) return 0;
    str.l = 0, str.m = *max, str.s = (char*)*data;
    if (*s == '\t') ++s;
    for (p = q = s;; ++p) {
        if (*p == 0 || *p == '\t') {
            int c = *p;
            *p = 0;
            if (p - q >= 5 && q[2] == ':' && q[4] == ':' && (q[3] == 'A' || q[3] == 'i' || q[3] == 'f' || q[3] == 'Z' || q[3] == 'B')) {
                int type = q[3];
                kputsn_(q, 2, &str);
                q += 5;
                if (type == 'A') {
                    kputc_('A', &str);
                    kputc_(*q, &str);
                } else if (type == 'i') {
                    int64_t x;
                    x = strtoll(q, &q, 10);
                    kputc_(type, &str); kputsn_((char*)&x, 8, &str);
                } else if (type == 'f') {
                    double x;
                    x = strtod(q, &q);
                    kputc_('f', &str); kputsn_(&x, 8, &str);
                } else if (type == 'Z') {
                    kputc_('Z', &str); kputsn_(q, p - q + 1, &str); // note that this include the trailing NULL
                } else if (type == 'B') {
                    type = *q++; // q points to the first ',' following the typing byte
                    if (p - q >= 2 && (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I' || type != 'f')) {
                        int32_t n;
                        char *r;
                        for (r = q, n = 0; *r; ++r)
                            if (*r == ',') ++n;
                        kputc_('B', &str); kputc_(type, &str); kputsn_(&n, 4, &str);
                        // TODO: to evaluate which is faster: a) aligned array and then memmove(); b) unaligned array; c) kputsn_()
                        if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0);  kputc_(x, &str); }
                        else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtol(q + 1, &q, 0);  kputc_(x, &str); }
                        else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0);  kputsn_(&x, 2, &str); }
                        else if (type == 'S') while (q + 1 < p) { uint16_t x = strtol(q + 1, &q, 0);  kputsn_(&x, 2, &str); }
                        else if (type == 'i') while (q + 1 < p) { int64_t  x = strtoll(q + 1, &q, 0); kputsn_(&x, 8, &str); }
                        else if (type == 'I') while (q + 1 < p) { uint64_t x = strtoll(q + 1, &q, 0); kputsn_(&x, 8, &str); }
                        else if (type == 'f') while (q + 1 < p) { double   x = strtod(q + 1, &q);     kputsn_(&x, 8, &str); }
                    }
                } // should not be here, as we have tested all types
            }
            q = p + 1;
            if (c == 0) break;
        }
    }
    if (str.l > 0 && str.l == str.m) ks_resize(&str, str.l + 1);
    if (str.s) str.s[str.l] = 0;
    *max = str.m, *data = (uint8_t*)str.s;
    return str.l;
}

#define PARSE_A_ERR -1
#define PARSE_Q_ERR -2
#define PARSE_S_ERR -3
#define PARSE_L_ERR -4

static int asg_parse_fa_hdr(asg_t *g, char *s, asg_seg_t **seg)
{
    uint64_t i;
    char *c = s;
    while (*c != 0 && !isspace(*c))
        ++c;
    *c = 0;
    if (c - s == 1) // empty
        return PARSE_A_ERR;
    i = asg_add_seg(g, s + 1, 0);
    *seg = &g->seg[i];
    return 0;
}

static inline int gfa_parse_S(asg_t *g, char *s, int min_cov)
{
    if (*s != 'S') return PARSE_S_ERR;

    int i, c, is_ok;
    char *p, *q, *seg, *seq, *rest;
    uint32_t sid, len;
    
    seg = seq = rest = 0;
    len = 0;
    for (i = 0, p = q = s + 2;; ++p) {
        if (*p == 0 || *p == '\t') {
            c = *p;
            *p = 0;
            if (i == 0) seg = q;
            else if (i == 1) {
                seq = q[0] == '*'? 0 : strdup(q);
                is_ok = 1, rest = c? p + 1 : 0;
                break;
            }
            ++i, q = p + 1;
            if (c == 0) break;
        }
    }

    if (is_ok) { // all mandatory fields read
        // parse sequence coverage if presented
        int l_aux, m_aux = 0;
        uint32_t LN = 0;
        uint8_t *aux = 0, *s_LN = 0;
        asg_seg_t *s;
        l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
        s_LN = l_aux? gfa_aux_get(l_aux, aux, "LN") : 0;
        if (s_LN && s_LN[0] == 'i')
            LN = *(int64_t*)(s_LN + 1);
        if (seq == 0) {
            if (LN > 0) len = LN;
        } else {
            len = strlen(seq);
        }
        if (LN > 0 && len != LN)
            fprintf(stderr, "[W::%s] for segment '%s', LN:i:%u tag is different from sequence length %d\n", __func__, seg, LN, len);
        sid = asg_add_seg(g, seg, 0);
        s = &g->seg[sid];
        s->len = len, s->seq = seq;
        if (l_aux > 0) {
            uint8_t *s_SBP_COV = 0, *s_SEQ_COV = 0;
            double dv = 0;
            char tag[2];
            if (TAG_SBP_COV[0] != 0) {
                memcpy(tag, TAG_SBP_COV, 2);
                s_SBP_COV = gfa_aux_get(l_aux, aux, tag);
                if (s_SBP_COV && *s_SBP_COV == TAG_SBP_COV[3]) {
                    dv = gfa_aux_decimal_value(s_SBP_COV);
                    s->cov = len > 0? dv/len : dv;
                } else {
                    fprintf(stderr, "[W::%s] for segment '%s', %s tag is absent\n", __func__, seg, tag);
                }
            } else if (TAG_SEQ_COV[0] != 0) {
                memcpy(tag, TAG_SEQ_COV, 2);
                s_SEQ_COV = gfa_aux_get(l_aux, aux, tag);
                if (s_SEQ_COV && *s_SEQ_COV == TAG_SEQ_COV[3]) {
                    s->cov = gfa_aux_decimal_value(s_SEQ_COV);
                } else {
                    fprintf(stderr, "[W::%s] for segment '%s', %s tag is absent\n", __func__, seg, tag);
                }
            } else {
                // check KC and FC
                s_SBP_COV = gfa_aux_get(l_aux, aux, "KC");
                if (s_SBP_COV && *s_SBP_COV == 'i') {
                    dv = *(int64_t*)(s_SBP_COV + 1);
                } else {
                    s_SBP_COV = gfa_aux_get(l_aux, aux, "FC");
                    if (s_SBP_COV && *s_SBP_COV == 'i')
                        dv = *(int64_t*)(s_SBP_COV + 1);
                }
                s->cov = len > 0? dv/len : dv;
            }
        }
        if (s->cov == 0 && min_cov) {
            fprintf(stderr, "[W::%s] the coverage of segment '%s' is zero\n", __func__, seg);
            s->cov = 1;
        }
        free(aux);
    } else return PARSE_S_ERR;

    return 0;
}

static int gfa_parse_L(asg_t *g, char *s, int min_cov)
{
    if (*s != 'L') return PARSE_L_ERR;

    int i, is_ok = 0;
    char *p, *q, *segv, *segw, *rest;
    uint64_t ov, ow, oriv, oriw;
    
    segv = segw = rest = 0;
    ov = ow = oriv = oriw = UINT64_MAX;
    for (i = 0, p = q = s + 2;; ++p) {
        if (*p == 0 || *p == '\t') {
            int c = *p;
            *p = 0;
            if (i == 0) {
                segv = q;
            } else if (i == 1) {
                if (*q != '+' && *q != '-') return PARSE_L_ERR;
                oriv = (*q != '+');
            } else if (i == 2) {
                segw = q;
            } else if (i == 3) {
                if (*q != '+' && *q != '-') return PARSE_L_ERR;
                oriw = (*q != '+');
            } else if (i == 4) {
                if (*q == '*') {
                    ov = ow = 0;
                } else if (*q == ':') {
                    ov = UINT64_MAX;
                    ow = isdigit(*(q+1))? strtoul(q+1, &q, 10) : UINT64_MAX;
                } else if (isdigit(*q)) {
                    char *r;
                    ov = strtol(q, &r, 10);
                    if (isupper(*r)) { // CIGAR
                        ov = ow = 0;
                        do {
                            long l;
                            l = strtol(q, &q, 10);
                            if (*q == 'M' || *q == 'D' || *q == 'N') ov += l;
                            if (*q == 'M' || *q == 'I' || *q == 'S') ow += l;
                            ++q;
                        } while (isdigit(*q));
                    } else if (*r == ':') { // overlap lengths
                        ow = isdigit(*(r+1))? strtoul(r+1, &r, 10) : UINT64_MAX;
                    } else break;
                } else break;
                is_ok = 1, rest = c? p + 1 : 0;
                break;
            }
            ++i, q = p + 1;
            if (c == 0) break;
        }
    }
    if (i == 4 && is_ok == 0) ov = ow = 0, is_ok = 1; // no overlap field
    if (is_ok) {
        uint64_t v, w;
        int l_aux, m_aux = 0;
        uint8_t *aux = 0;
        asmg_arc_t *arc;
        v = asg_add_seg(g, segv, 1) << 1 | oriv;
        w = asg_add_seg(g, segw, 1) << 1 | oriw;
        arc = asmg_arc_add(g->asmg, v, w, 0, ov, UINT64_MAX, 0, 0);
        l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
        if (l_aux) {
            uint8_t *s_ARC_COV = 0;
            char tag[2];
            if (TAG_ARC_COV[0] != 0) {
                memcpy(tag, TAG_ARC_COV, 2);
                s_ARC_COV = gfa_aux_get(l_aux, aux, tag);
                if (s_ARC_COV && *s_ARC_COV == TAG_ARC_COV[3])
                    arc->cov = gfa_aux_decimal_value(s_ARC_COV);
                else
                    fprintf(stderr, "[W::%s] for arc '%s%c' -> '%s%c', %s tag is absent\n", __func__, segv, "+-"[oriv], segw, "+-"[oriw], tag);
            } else {
                // check EC
                s_ARC_COV = gfa_aux_get(l_aux, aux, "EC");
                if (s_ARC_COV && *s_ARC_COV == 'i')
                    arc->cov = *(int64_t*)(s_ARC_COV + 1);
            }
        }
        if (arc->cov == 0 && min_cov) {
            fprintf(stderr, "[W::%s] the coverage of arc '%s%c' -> '%s%c' is zero\n", __func__, segv, "+-"[oriv], segw, "+-"[oriw]);
            arc->cov = 1;
        }
        free(aux);
    } else return PARSE_L_ERR;
    return 0;
}

static void asg_finalize_asmg(asg_t *g)
{
    if (g->asmg == 0)
        MYCALLOC(g->asmg, 1);
    asmg_t *asmg = g->asmg;
    uint64_t i;
    asmg_vtx_t *vtx;
    asmg->n_vtx = asmg->m_vtx = g->n_seg;
    MYCALLOC(asmg->vtx, g->n_seg);
    for (i = 0; i < g->n_seg; ++i) {
        vtx = &asmg->vtx[i];
        // will never need vtx->a
        // seg and vtx indices are consistent
        //vtx->n = 1;
        //MYMALLOC(vtx->a, 1);
        //vtx->a[0] = i;
        vtx->len = g->seg[i].len;
        vtx->cov = g->seg[i].cov;
    }
    asmg_finalize(asmg, 0);
}

asg_t *asg_read(const char *fn, int min_s_cov, int min_l_cov)
{
    int dret, ret, is_fa, is_fq, is_gfa;
    gzFile fp;
    kstream_t *ks;
    asg_seg_t *fa_seg;
    uint64_t lineno;
    asg_t *g;
    kstring_t s = {0, 0, 0}, fa_seq = {0, 0, 0};

    fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
    if (fp == 0) return 0;
    ks = ks_init(fp);
    g = asg_init();
    lineno = 0;
    fa_seg = 0;
    is_fa = is_fq = is_gfa = 0;
    
    while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
        ++lineno;
        if (s.l == 0) // empty line
            continue;

        ret = 0;
        if (!is_gfa && s.s[0] == '>') { // FASTA header
            is_fa = 1;
            if (fa_seg) asg_update_seg_seq(fa_seg, fa_seq.l, fa_seq.s);
            // parse header
            asg_parse_fa_hdr(g, s.s, &fa_seg);
            fa_seq.l = 0;
        } else if (!is_gfa && s.s[0] == '@') { // FASTQ header
            is_fq = 1;
            // parse and write header
            asg_parse_fa_hdr(g, s.s, &fa_seg);
            // parse sequence
            if(ks_getuntil(ks, KS_SEP_LINE, &s, &dret) < 0)
                ret = PARSE_Q_ERR;
            else
                asg_update_seg_seq(fa_seg, s.l, s.s);
            ++lineno;
            // skip quality score lines
            if (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) < 0 ||
                    ks_getuntil(ks, KS_SEP_LINE, &s, &dret) < 0)
                ret = PARSE_Q_ERR;
            lineno += 2;
        } else if (is_fa) { // FASTA sequence line
            // append sequence to seg
            kputsn(s.s, s.l, &fa_seq);
        } else {
            is_gfa = 1;
            if (s.s[0] == 'S')
                ret = gfa_parse_S(g, s.s, min_s_cov);
            else if (s.s[0] == 'L')
                ret = gfa_parse_L(g, s.s, min_l_cov);
        }

        if (ret < 0) {
            fprintf(stderr, "[E::%s] failed to parse %s file: %c-line at line %lu (error code %d)\n",
                    __func__, is_fa? "FASTA" : (is_fq? "FASTQ" : "GFA"), s.s[0], lineno, ret);
            exit(EXIT_FAILURE);
        }
    }
    if (is_fa && fa_seg) asg_update_seg_seq(fa_seg, fa_seq.l, fa_seq.s);

    free(fa_seq.s);
    free(s.s);
    ks_destroy(ks);
    gzclose(fp);

    // add vtx, fix symmetric arcs, sort and index etc
    asg_finalize_asmg(g);

    return g;
}

void asg_stat(asg_t *asg, FILE *fo)
{
    uint64_t i, nv, n_vtx, n_seg, max_deg, tot_seg_len, n_link, n_arc, tot_deg;
    asmg_t *g;

    n_vtx = n_seg = max_deg = tot_seg_len = n_link = n_arc = tot_deg = 0;
    g = asg->asmg;
    for (i = 0; i < asg->n_seg; ++i) {
        if (g->vtx[i].del) continue;
        tot_seg_len += asg->seg[i].len;
        ++n_seg;
    }
    fprintf(fo, "Number of segments: %lu\n", n_seg);
    fprintf(fo, "Total segment length: %lu\n", tot_seg_len);
    if (n_seg) fprintf(fo, "Average segment length: %.3f\n", (double) tot_seg_len / n_seg);

    for (i = 0; i < g->n_arc; ++i) {
        if (!g->arc[i].del) {
            ++n_arc;
            if (!g->arc[i].comp)
                ++n_link;
        }
    }
    fprintf(fo, "Number of links: %lu\n", n_link);
    fprintf(fo, "Number of arcs: %lu\n", n_arc);

    n_vtx = asmg_vtx_n(g);
    for (i = 0; i < n_vtx; ++i) {
        nv = asmg_arc_n1(g, i);
        if (nv > max_deg)
            max_deg = nv;
        tot_deg += nv;
    }
    fprintf(fo, "Max degree: %lu\n", max_deg);
    if (n_seg > 0) fprintf(fo, "Average degree: %.3f\n", (double) tot_deg / n_seg / 2);
}

void asg_print(asg_t *g, FILE *fo, int no_seq)
{
    uint32_t i, cov;
    uint64_t k;
    asmg_t *asmg = g->asmg;

    fprintf(fo, "H\tVN:Z:1.0\n");
    for (i = 0; i < g->n_seg; ++i) {
        const asg_seg_t *s = &g->seg[i];
        // seg and vtx indices are interchangeable
        if (asmg && asmg->vtx[i].del) continue;
        cov = asmg? asmg->vtx[i].cov : s->cov;
        fprintf(fo, "S\t%s\t", s->name);
        if (s->seq && !no_seq) fputs(s->seq, fo);
        else fputc('*', fo);
        fprintf(fo, "\tLN:i:%u\tKC:i:%lu\tSC:f:%.3f\n", s->len, (uint64_t)s->len*cov, (double)cov);
    }

    if (asmg == 0) return;
    for (k = 0; k < asmg->n_arc; ++k) {
        const asmg_arc_t *a = &asmg->arc[k];
        if (a->del || a->comp) continue;
        fprintf(fo, "L\t%s\t%c\t%s\t%c\t%luM\tEC:i:%u\n", g->seg[a->v>>1].name, "+-"[a->v&1], 
                g->seg[a->w>>1].name, "+-"[a->w&1], a->ls, a->cov);
    }
}

void asg_print_fa(asg_t *g, FILE *fo, int line_wd)
{
    uint64_t i, l;
    for (i = 0; i < g->n_seg; ++i) {
        if (g->asmg && g->asmg->vtx[i].del)
            continue;
        if (g->seg[i].seq == 0)
            fprintf(stderr, "[W::%s] skip empty sequence: %s\n", __func__, g->seg[i].name);
        fprintf(fo, ">%s\n", g->seg[i].name);
        l = 0;
        put_chars(g->seg[i].seq, g->seg[i].len, 0, 0, fo, &l, line_wd);
        if (l % line_wd != 0) fputc('\n', fo);
    }
}

char **asg_vtx_name_list(asg_t *g, uint64_t *_n)
{
    if (_n) *_n = 0;
    uint64_t i, n = 0;
    uint64_t *vlist = asmg_vtx_list(g->asmg, &n);
    char **names = 0;
    if (n) {
        MYMALLOC(names, n);
        for (i = 0; i < n; ++i)
            names[i] = g->seg[vlist[i]].name;
    }
    free(vlist);
    if (_n) *_n = n;
    return names;
}


