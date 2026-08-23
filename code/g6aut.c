/* Read graph6 on stdin; write "graph6 aut_order" per line. */
#define MAXN 64
#include "gtools.h"

int main(void)
{
    char *line = NULL;
    size_t cap = 0;
    ssize_t len;
    graph g[MAXN * SETWORDSNEEDED(MAXN)];
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    int n, m;

    options.writeautoms = FALSE;
    options.getcanon = FALSE;

    while ((len = getline(&line, &cap, stdin)) >= 0) {
        if (len > 0 && line[len - 1] == '\n')
            line[--len] = '\0';
        if (len > 0 && line[len - 1] == '\r')
            line[--len] = '\0';
        if (len == 0 || line[0] == '>' || line[0] == '#')
            continue;
        n = graphsize(line);
        if (n < 1 || n > MAXN) {
            fprintf(stderr, "bad n=%d\n", n);
            exit(1);
        }
        m = SETWORDSNEEDED(n);
        nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
        stringtograph(line, g, m);
        densenauty(g, lab, ptn, orbits, &options, &stats, m, n, NULL);
        fputs(line, stdout);
        fputc(' ', stdout);
        writegroupsize(stdout, stats.grpsize1, stats.grpsize2);
        fputc('\n', stdout);
    }
    free(line);
    return 0;
}
