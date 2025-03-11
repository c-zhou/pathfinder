/*********************************************************************************
 * MIT License                                                                   *
 *                                                                               *
 * Copyright (c) 2022 Chenxi Zhou <chnx.zhou@gmail.com>                          *
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
 * 03/08/22 - Chenxi Zhou: Created                                               *
 *                                                                               *
 *********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "ketopt.h"
#include "kvec.h"

#include "path.h"
#include "graph.h"

#define PATHFINDER_VERSION "0.1"

#define DEFAULT_MAX_PATH 10000000
#define DEFAULT_MAX_COPY 10

int VERBOSE = 0;

typedef struct { size_t n, m; uint32_t *a; } v_u32_t;
typedef struct { size_t n, m; v_u32_t *a; } lv_u32_t;

static void lv_u32_destroy(lv_u32_t *lv)
{
    size_t i;
    for (i = 0; i < lv->n; ++i)
        kv_destroy(lv->a[i]);
    kv_destroy(*lv);
}

static lv_u32_t parse_subgraph(asg_t *asg)
{
    uint32_t i, j, n_seg, *vlist, nv;
    lv_u32_t lv_list;
    asmg_t *g;
    uint8_t *visited;

    g = asg->asmg;
    n_seg = asg->n_seg;
    kv_init(lv_list);
    MYCALLOC(visited, n_seg);
    for (i = 0; i < n_seg; ++i) {
        if (visited[i] || g->vtx[i].del) continue;
        vlist = asmg_subgraph(g, &i, 1, 0, 0, &nv, 0);
        assert(nv > 0);
        v_u32_t lv = {nv, nv, vlist};
        kv_push(v_u32_t, lv_list, lv);
        for (j = 0; j < nv; ++j)
            visited[vlist[j]] = 1;
    }
    free(visited);

    return lv_list;
}

static int64_t parse_node(asg_t *asg, char *name)
{
    if (!name) return -1;
    int l;
    char s;
    int64_t v, d;
    
    v = asg_name2id(asg, name);
    if (v != UINT32_MAX)
        return (v<<1);

    l = strlen(name);
    s = '\0';
    d = 0;
    if (l > 1) {
        if (name[l-1] == '+') {
            s = '+';
            name[l-1] = '\0';
        } else if (name[l-1] == '-') {
            s = '-';
            name[l-1] = '\0';
            d = 1;
        }
    }

    if (s != '\0') {
        v = asg_name2id(asg, name);
        if (v == UINT32_MAX) {
            fprintf(stderr, "[E::%s] cannot find the node '%s' in the graph\n", __func__, name);
            exit (1);
        }
        v = v << 1 | d;
        name[l-1] = s;    
    }

    return v;
}

static v_u32_t *find_subgraph(lv_u32_t *lv, uint32_t v)
{
    uint32_t i, j, n, *a;

    for (i = 0; i < lv->n; i++) {
        a = lv->a[i].a;
        n = lv->a[i].n;
        for (j = 0; j < n && v != a[j]; j++);
        if (j < n) return (&lv->a[i]);
    }

    return NULL;
}

static void print_copy_number(asg_t *asg, double avg_coverage, int *copy_number, int mstr)
{
    uint32_t i, j, nv, v, *copy;
    asg_seg_t *seg;
    fprintf(stderr, "[M::%s] estimated haplotype coverage: %.2f\n", __func__, avg_coverage);
    fprintf(stderr, "[M::%s] copy number estimation results\n", __func__);
    fprintf(stderr, "[M::%s] node name length coverage copyEstd\n", __func__);
    for (i = j = 0, nv = asg->n_seg; i < nv; i++) {
        if (asg->asmg->vtx[i].del)
            continue;
        seg = &asg->seg[i];
        fprintf(stderr, "[M::%s] [%4u] %*s %10u %10u %4d\n", __func__, ++j, mstr, seg->name, seg->len, seg->cov, copy_number[i]);
    }
}

static void print_path(path_t *path, asg_t *asg, int *copy_number, char *s_source, char *s_target, int mstr, FILE *out)
{
    uint32_t i, j, nv, v, *copy;
    asg_seg_t *seg;
    
    MYCALLOC(copy, asg->n_seg);
    for (i = 0, nv = path->nv; i < nv; i++)
        copy[path->v[i]>>1] += 1;
    
    fprintf(out, "SUBGRAPH node name length coverage copyEstd copyIncl copyDiff\n");
    for (i = j = 0, nv = asg->n_seg; i < nv; i++) {
        if (asg->asmg->vtx[i].del)
            continue;
        seg = &asg->seg[i];
        fprintf(out, "[%4u] %*s %10u %10u %4d %4u %4d\n", ++j, mstr, seg->name, seg->len, seg->cov, copy_number[i], copy[i], (int32_t)copy[i]-copy_number[i]);
    }
    fprintf(out, "PATH node name [source=%s target=%s nv=%u len=%u wlen=%.0f ec=%u circ=%s]\n", 
        s_source? s_source : "NULL", s_target? s_target : "NULL", path->nv, path->len, path->wlen, path->ec, path->circ? "true" : "false");
    for (i = 0, nv = path->nv; i < nv; i++) {
        v = path->v[i];
        fprintf(out, "[%4u] %*s%c\n", i+1, mstr, asg->seg[v>>1].name, "+-"[v&1]);
    }
    free(copy);
}

static uint64_t *parse_domi_path(uint64_t *domi_tree, uint64_t source, uint64_t target, asmg_t *g, int *copy_number, int *n)
{
    kvec_t(uint64_t) path;
    kv_init(path);

    kv_push(uint64_t, path, target);
    while ((target = domi_tree[target]) != UINT64_MAX) {
        if (target == source)
            break;
        if (copy_number[target>>1] == 1 &&
            asmg_arc_n1(g, target) == 1 &&
            asmg_arc_n1(g, target^1) == 1)
            kv_push(uint64_t, path, target);
    }
    if (target == source)
        kv_push(uint64_t, path, target);
    else {
        path.a[1] = source;
        path.n = 2;
    }
    *n = path.n;
    return path.a;
}

static int pathfinder(char *asg_file, int min_copy, int max_copy, double min_cfrac, int max_path, int do_part, int do_adjust, FILE *out_file, char *s_source, char *s_target, int VERBOSE)
{   
    asg_t *asg;
    int64_t source, target;
    int i, j, mstr, ret = 0;
    
    asg = asg_read(asg_file);
    if (asg == 0) {
        fprintf(stderr, "[E::%s] failed to read the graph: %s\n", __func__, asg_file);
        return 1;
    }

    mstr = 0;
    for (i = 0; i < asg->n_seg; i++)
        if ((j = strlen(asg->seg[i].name)) > mstr)
            mstr = j;

    source = parse_node(asg, s_source);
    target = parse_node(asg, s_target);

    lv_u32_t subgraph_nodes = parse_subgraph(asg);
    v_u32_t *sub_gs = NULL;
    int sub_gn = 0;
    // locate source and target node
    if (source >= 0) {
        sub_gs = find_subgraph(&subgraph_nodes, source>>1);
        sub_gn = 1;
        if (target >= 0 && find_subgraph(&subgraph_nodes, target>>1) != sub_gs) {
            fprintf(stderr, "[E::%s] nodes %s and %s are not in the same subgraph\n", __func__, s_source, s_target);
            ret = 1;
            goto do_clean;
        }
    } else {
        sub_gs = subgraph_nodes.a;
        sub_gn = subgraph_nodes.n;
    }

    // now process each subgraph
    for (i = 0; i < sub_gn; i++) {
        int *copy_number = 0;
        uint32_t *nodes = (sub_gs + i)->a;
        uint32_t nnode = (sub_gs + i)->n;
        
        fprintf(stderr, "[M::%s] process subgraph [%d]:", __func__, i+1);
        for (j = 0; j < nnode; j++)
            fprintf(stderr, " %s", asg->seg[nodes[j]].name);
        fputc('\n', stderr);

        asg_t *asg_copy = asg_make_copy(asg);
        // modify the grpah to get the subgraph
        asmg_subgraph(asg_copy->asmg, nodes, nnode, 0, 0, 0, 1);
        
        double avg_coverage, adjusted_avg_coverage;
        // initial guess of sequence copy numbers
        avg_coverage = graph_sequence_coverage_precise(asg_copy, 0, min_copy, max_copy, &copy_number);
        if (VERBOSE > 1) {
            fprintf(stderr, "[M::%s] initial copy number estimation\n", __func__);
            print_copy_number(asg_copy, avg_coverage, copy_number, mstr);
        }
        // adjust estimation considering graph structure
        if (do_adjust) {
            adjust_sequence_copy_number_by_graph_layout(asg_copy, avg_coverage, &adjusted_avg_coverage, copy_number, max_copy, 10);
            if (VERBOSE > 1) {
                fprintf(stderr, "[M::%s] adjusted copy number estimation\n", __func__);
                print_copy_number(asg_copy, adjusted_avg_coverage, copy_number, mstr);
            }
            avg_coverage = adjusted_avg_coverage;
        }
        
        // now do path finding
        path_t *path = NULL;
        if (!do_part)
            path = graph_path_finder_dfs(asg_copy, copy_number, min_cfrac, source, target, max_path);
        else {
            uint64_t s, t, *domi_tree, *domi_path;
            int domi_n = 0;
            domi_tree = asmg_lengauer_tarjans_domtree(asg_copy->asmg, source);
            domi_path = parse_domi_path(domi_tree, source, target, asg->asmg, copy_number, &domi_n);
            if (domi_n > 2) {
                fprintf(stderr, "[M::%s] graph partitioned into %d parts\n", __func__, domi_n-1);
                path_t **paths, *p;
                int success = 1;
                MYCALLOC(paths, domi_n-1);
                for (j = domi_n; --j;) {
                    s = domi_path[j];
                    t = domi_path[j-1];
                    fprintf(stderr, "[M::%s] processing partition %s%c -> %s%c\n", __func__,
                            asg_copy->seg[s>>1].name, "+-"[s&1], asg_copy->seg[t>>1].name, "+-"[t&1]);
                    p = graph_path_finder_dfs(asg_copy, copy_number, min_cfrac, s, t, max_path);
                    if (!p) {
                        success = 0;
                        break;
                    } else paths[j-1] = p;
                }
                if (success) {
                    fprintf(stderr, "[M::%s] concatenate %d graph partitions into one\n", __func__, domi_n-1);
                    uint32_t *vs, nv = 0;
                    for (j = 0; j < domi_n-1; j++)
                        nv += paths[j]->nv;
                    nv = nv + 2 - domi_n; // remove overlaps
                    MYMALLOC(vs, nv);
                    nv = 0;
                    for (j = domi_n; --j; ) {
                        p = paths[j-1];
                        memcpy(vs+nv, p->v, sizeof(uint32_t)*p->nv);
                        nv += p->nv - 1;
                    }
                    nv += 1;
                    MYCALLOC(path, 1);
                    path->v = vs;
                    path->nv = nv;
                    path_report(asg_copy, path);
                }
                for (j = 0; j < domi_n - 1; j++) {
                    path_destroy(paths[j]);
                    free(paths[j]);
                }
                free(paths);
            } else {
                fprintf(stderr, "[M::%s] graph partition is not possible\n", __func__);
                path = graph_path_finder_dfs(asg_copy, copy_number, min_cfrac, source, target, max_path);
            }
            free(domi_tree);
            free(domi_path);
        }
            
        if (path) {
            if (out_file)
                print_path(path, asg, copy_number, s_source, s_target, mstr, out_file);
            else
                print_path(path, asg, copy_number, s_source, s_target, mstr, stderr);
        } else {
            print_copy_number(asg_copy, avg_coverage, copy_number, mstr);
            fprintf(stderr, "[W::%s] no valid path found in the graph\n", __func__);
        }
        asg_destroy(asg_copy);
        path_destroy(path);
        free(path);
        free(copy_number);
    }

do_clean:
    lv_u32_destroy(&subgraph_nodes);
    asg_destroy(asg);

    return ret;
}

static ko_longopt_t long_options[] = {
    { "edge-c-tag",     ko_required_argument, 301 },
    { "kmer-c-tag",     ko_required_argument, 302 },
    { "seq-c-tag",      ko_required_argument, 303 },
    { "max-copy",       ko_required_argument, 'c' },
    { "max-path",       ko_required_argument, 'N' },
    { "verbose",        ko_required_argument, 'v' },
    { "version",        ko_no_argument,       'V' },
    { "help",           ko_no_argument,       'h' },
    { 0, 0, 0 }
};

int main(int argc, char *argv[])
{
    const char *opt_str = "ac:d:hN:o:pv:V";
    ketopt_t opt = KETOPT_INIT;
    int c, min_copy, max_copy, max_path, do_part, do_adjust, ret = 0;
    FILE *fp_help;
    char *out_file, *ec_tag, *kc_tag, *sc_tag, *source, *target;
    double min_cfrac;

    sys_init();

    out_file = NULL;
    source = NULL;
    target = NULL;
    fp_help = stderr;
    ec_tag = 0;
    kc_tag = 0;
    sc_tag = 0;
    do_part = 0;
    do_adjust = 0;
    min_cfrac = 1.;
    min_copy = 1;
    max_copy = DEFAULT_MAX_COPY;
    max_path = DEFAULT_MAX_PATH;

    while ((c = ketopt(&opt, argc, argv, 1, opt_str, long_options)) >=0) {
        if (c == 'c') min_copy = atoi(opt.arg);
        else if (c == 'C') max_copy = atoi(opt.arg);
        else if (c == 'd') min_cfrac = atof(opt.arg);
        else if (c == 'N') max_path = atoi(opt.arg);
        else if (c == 'p') do_part = 1;
        else if (c == 'a') do_adjust = 1;
        else if (c == 'o') out_file = opt.arg;
        else if (c == 301) ec_tag = opt.arg;
        else if (c == 302) kc_tag = opt.arg;
        else if (c == 303) sc_tag = opt.arg;
        else if (c == 'v') VERBOSE = atoi(opt.arg);
        else if (c == 'h') fp_help = stdout;
        else if (c == 'V') {
            puts(PATHFINDER_VERSION);
            return 0;
        }
        else if (c == '?') {
            fprintf(stderr, "[E::%s] unknown option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
        else if (c == ':') {
            fprintf(stderr, "[E::%s] missing option: \"%s\"\n", __func__, argv[opt.i - 1]);
            return 1;
        }
    }

    if (argc == opt.ind || fp_help == stdout) {
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Usage: pathfinder [options] <file>[.gfa[.gz]] [<source>[+|-] [<target>[+|-]]]\n");
        fprintf(fp_help, "Options:\n");
        fprintf(fp_help, "  -c INT               minimum copy number of sequences to consider [%d]\n", min_copy);
        fprintf(fp_help, "  -C INT               maximum copy number of sequences to consider [%d]\n", max_copy);
        fprintf(fp_help, "  -d FLOAT             prefer a circular path if length >= FLOAT * linear length [%.2f]\n", min_cfrac);
        fprintf(fp_help, "  -p                   do graph partitioning if possible\n");
        fprintf(fp_help, "  -a                   adjust seuqnece copy number estimation by graph structure\n");
        fprintf(fp_help, "  -N INT               maximum number of graph paths to explore [%d]\n", max_path);
        fprintf(fp_help, " \n");
        fprintf(fp_help, "  --edge-c-tag  STR    edge coverage tag in the GFA file [EC:i] \n");
        fprintf(fp_help, "  --kmer-c-tag  STR    kmer coverage tag in the GFA file [KC:i] \n");
        fprintf(fp_help, "  --seq-c-tag   STR    sequence coverage tag in the GFA file [SC:f]\n");
        fprintf(fp_help, " \n");
        fprintf(fp_help, "  -o FILE              write output to a file [stdout]\n");
        fprintf(fp_help, "  -v INT               verbose level [%d]\n", VERBOSE);
        fprintf(fp_help, "  --version            show version number\n");
        fprintf(fp_help, "\n");
        fprintf(fp_help, "Example: ./pathfinder input.gfa\n\n");
        return fp_help == stdout? 0 : 1;
    }

    if (argc - opt.ind < 1) {
        fprintf(stderr, "[E::%s] missing input: please specify the GFA file\n", __func__);
        return 1;
    }

    if (argc - opt.ind > 1)
        source = argv[opt.ind + 1];
    if (argc - opt.ind > 2)
        target = argv[opt.ind + 2];

    if (ec_tag != 0) {
        if(is_valid_gfa_tag(ec_tag)) {
            memcpy(TAG_ARC_COV, ec_tag, 4);
        } else {
            fprintf(stderr, "[E::%s] invalid GFA tag (Regexp: [A-Za-z][A-Za-z0-9]:[A|i|f|Z|B]): %s\n", __func__, ec_tag);
            return 1;
        }
    }

    if (kc_tag != 0) {
        if(is_valid_gfa_tag(kc_tag)) {
            memcpy(TAG_SBP_COV, kc_tag, 4);
        } else {
            fprintf(stderr, "[E::%s] invalid GFA tag (Regexp: [A-Za-z][A-Za-z0-9]:[A|i|f|Z|B]): %s\n", __func__, kc_tag);
            return 1;
        }
    }

    if (sc_tag != 0) {
        if(is_valid_gfa_tag(sc_tag)) {
            memcpy(TAG_SEQ_COV, sc_tag, 4);
        } else {
            fprintf(stderr, "[E::%s] invalid GFA tag (Regexp: [A-Za-z][A-Za-z0-9]:[A|i|f|Z|B]): %s\n", __func__, sc_tag);
            return 1;
        }
    }

    if (do_part && (source == NULL || target == NULL)) {
        do_part = 0;
        fprintf(stderr, "[W::%s] graph partitioning only applied with source and target nodes set\n", __func__);
    }

    if (out_file != NULL && strcmp(out_file, "-") != 0) {
        if (freopen(out_file, "wb", stdout) == NULL) {
            fprintf(stderr, "[ERROR]\033[1;31m failed to write the output to file '%s'\033[0m: %s\n", out_file, strerror(errno));
            exit(1);
        }
    }

    ret = pathfinder(argv[opt.ind], min_copy, max_copy, min_cfrac, max_path, do_part, do_adjust, stdout, source, target, VERBOSE);
    
    if (ret) {
        fprintf(stderr, "[E::%s] failed to analysis the GFA file\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (fflush(stdout) == EOF) {
        fprintf(stderr, "[E::%s] failed to write the results\n", __func__);
        exit(EXIT_FAILURE);
    }

    if (VERBOSE >= 0) {
        fprintf(stderr, "[M::%s] Version: %s\n", __func__, PATHFINDER_VERSION);
        fprintf(stderr, "[M::%s] CMD:", __func__);
        int i;
        for (i = 0; i < argc; ++i)
            fprintf(stderr, " %s", argv[i]);
        fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    }

    return 0;
}
