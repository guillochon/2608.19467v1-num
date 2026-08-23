#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#ifndef MAX_SUBSET_N
#define MAX_SUBSET_N 20
#endif
#ifndef MAX_C5_N
#define MAX_C5_N 80
#endif

/* Stock GMP has mpq_mul/mpq_div but not these mixed-type helpers. */
static void mpq_mul_ui(mpq_t rop, const mpq_t op, unsigned long n) {
    mpz_mul_ui(mpq_numref(rop), mpq_numref(op), n);
    if (rop != op) mpz_set(mpq_denref(rop), mpq_denref(op));
    mpq_canonicalize(rop);
}
static void mpq_div_z(mpq_t rop, const mpq_t op, const mpz_t z) {
    if (rop != op) mpz_set(mpq_numref(rop), mpq_numref(op));
    mpz_mul(mpq_denref(rop), mpq_denref(op), z);
    mpq_canonicalize(rop);
}

static int number(const char *s, int *v) {
    char *e; long x; errno = 0; x = strtol(s, &e, 10);
    if (errno || !*s || *e || x < 0 || x > INT_MAX) return 0;
    *v = (int)x; return 1;
}
static void factorial(mpz_t z, int n) {
    int i; mpz_set_ui(z, 1); for (i = 2; i <= n; ++i) mpz_mul_ui(z, z, (unsigned long)i);
}
static unsigned long choose_ul(int n, int k) {
    unsigned long r = 1; int i;
    if (k < 0 || k > n) return 0;
    if (k > n-k) k = n-k;
    for (i = 1; i <= k; ++i) r = r * (unsigned long)(n-k+i) / (unsigned long)i;
    return r;
}
static int pop64(uint64_t x) { return __builtin_popcountll((unsigned long long)x); }
static int ctz64(uint64_t x) { return __builtin_ctzll((unsigned long long)x); }

/* Integer-scaled ordinary subset recurrence, used for complete bipartite graphs. */
static int subset_eval(int n, const uint64_t *adj, const mpz_t aut, mpq_t out) {
    size_t count, mask, t; int r, d;
    mpz_t *tab, p[MAX_SUBSET_N+1], coef[MAX_SUBSET_N+1], den, tmp, lcm, z;
    if (n < 1 || n > MAX_SUBSET_N || n >= 64 || n >= (int)(8*sizeof(size_t))) return 0;
    count = (size_t)1 << n;
    if (count > SIZE_MAX / sizeof(*tab)) return 0;
    fprintf(stderr, "subset recurrence: n=%d, %zu states\n", n, count);
    tab = malloc(count * sizeof(*tab)); if (!tab) return 0;
    for (t=0; t<count; ++t) mpz_init(tab[t]);
    for (r=0; r<=n; ++r) { mpz_init(p[r]); mpz_init(coef[r]); }
    mpz_inits(den,tmp,lcm,z,NULL); mpz_set_ui(p[1],1);
    for (r=2; r<=n; ++r) {
        mpz_set_ui(lcm,1);
        for (d=0; d<r; ++d) { mpz_set_ui(z,choose_ul(r-1,d)); mpz_lcm(lcm,lcm,z); }
        mpz_mul(p[r],p[r-1],lcm);
    }
    for (mask=1; mask<count; ++mask) {
        int k=pop64((uint64_t)mask);
        if (k==1) { mpz_set_ui(tab[mask],1); continue; }
        for (d=0; d<k; ++d) {
            mpz_divexact(coef[d],p[k],p[k-1]);
            mpz_divexact_ui(coef[d],coef[d],choose_ul(k-1,d));
        }
        for (t=mask; t; t &= t-1) {
            int v=ctz64((uint64_t)t); size_t sub=mask ^ ((size_t)1<<v);
            d=pop64(adj[v] & (uint64_t)sub);
            mpz_mul(tmp,tab[sub],coef[d]); mpz_add(tab[mask],tab[mask],tmp);
        }
    }
    factorial(den,n); mpz_mul(den,den,aut); mpz_mul(den,den,p[n]);
    mpq_set_num(out,tab[count-1]); mpq_set_den(out,den); mpq_canonicalize(out);
    mpz_clears(den,tmp,lcm,z,NULL);
    for (r=0; r<=n; ++r) { mpz_clear(p[r]); mpz_clear(coef[r]); }
    for (t=0; t<count; ++t) mpz_clear(tab[t]); free(tab); return 1;
}

static int make_bipartite(int a, int b, int *n, uint64_t *adj, mpz_t aut) {
    int i,j; mpz_t x,y;
    if (a<0 || b<0 || a+b<1 || a+b>MAX_SUBSET_N) return 0;
    *n=a+b; for(i=0;i<*n;++i) adj[i]=0;
    for(i=0;i<a;++i) for(j=a;j<*n;++j) { adj[i]|=UINT64_C(1)<<j; adj[j]|=UINT64_C(1)<<i; }
    mpz_inits(x,y,NULL); factorial(x,a); factorial(y,b); mpz_mul(aut,x,y);
    if(a==b) mpz_mul_ui(aut,aut,2); mpz_clears(x,y,NULL); return 1;
}

/* Exact quotient-state recurrence for positive C5 independent-set blow-ups. */
static int c5_eval(const int m[5], mpq_t out) {
    size_t stride[5], states=1, idx; int i, total=0, stab=0;
    mpq_t *f; mpz_t aut, bin, den; mpq_t term;
    for(i=0;i<5;++i) {
        if(m[i]<1 || total > MAX_C5_N-m[i]) return 0;
        stride[i]=states;
        if(states > SIZE_MAX/(size_t)(m[i]+1)) return 0;
        states *= (size_t)(m[i]+1); total += m[i];
    }
    if(states > SIZE_MAX/sizeof(*f)) return 0;
    fprintf(stderr,"C5 quotient recurrence: n=%d, %zu states\n",total,states);
    f=malloc(states*sizeof(*f)); if(!f) return 0;
    for(idx=0;idx<states;++idx) mpq_init(f[idx]);
    mpq_set_ui(f[0],1,1); mpz_inits(aut,bin,den,NULL); mpq_init(term);
    for(idx=1;idx<states;++idx) {
        int c[5], k=0; size_t q=idx;
        for(i=4;i>=0;--i) { c[i]=(int)(q/stride[i]); q%=stride[i]; k+=c[i]; }
        for(i=0;i<5;++i) if(c[i]) {
            int d=c[(i+4)%5]+c[(i+1)%5];
            mpz_bin_uiui(bin,(unsigned long)(k-1),(unsigned long)d);
            mpq_set(term,f[idx-stride[i]]);
            mpq_mul_ui(term,term,(unsigned long)c[i]);
            mpq_div_z(term,term,bin);
            mpq_add(f[idx],f[idx],term);
        }
    }
    mpz_set_ui(aut,1);
    for(i=0;i<5;++i) { factorial(den,m[i]); mpz_mul(aut,aut,den); }
    for(i=0;i<5;++i) {
        int s;
        for(s=-1;s<=1;s+=2) {
            int ok=1,j;
            for(j=0;j<5;++j) if(m[j]!=m[(i+s*j+10)%5]) ok=0;
            if(ok) ++stab;
        }
    }
    mpz_mul_ui(aut,aut,(unsigned long)stab); factorial(den,total); mpz_mul(den,den,aut);
    mpq_set(out,f[states-1]); mpq_div_z(out,out,den); mpq_canonicalize(out);
    mpq_clear(term); mpz_clears(aut,bin,den,NULL);
    for(idx=0;idx<states;++idx) mpq_clear(f[idx]); free(f); return 1;
}

/* Recurrence (6): O(n^2) exact evaluation of L(K_{a,n-a}). */
static int kbal_eval(int n, mpq_t out) {
    int t, a, a0;
    mpq_t *cur, *nxt, tmp;
    mpz_t fac, bin;
    if (n < 1) return 0;
    cur = malloc((size_t)(n + 1) * sizeof(*cur));
    nxt = malloc((size_t)(n + 1) * sizeof(*nxt));
    if (!cur || !nxt) { free(cur); free(nxt); return 0; }
    for (a = 0; a <= n; ++a) { mpq_init(cur[a]); mpq_init(nxt[a]); }
    mpq_init(tmp); mpz_inits(fac, bin, NULL);
    mpq_set_ui(cur[0], 1, 1);
    for (t = 1; t <= n; ++t) {
        mpq_set_ui(nxt[0], 1, 1);
        mpq_set_ui(nxt[t], 1, 1);
        for (a = 1; a < t; ++a) {
            mpz_bin_uiui(bin, (unsigned long)(t - 1), (unsigned long)(a - 1));
            mpq_set_z(tmp, bin);
            mpq_div(nxt[a], cur[a - 1], tmp);
            mpz_bin_uiui(bin, (unsigned long)(t - 1), (unsigned long)a);
            mpq_set_z(tmp, bin);
            mpq_div(tmp, cur[a], tmp);
            mpq_add(nxt[a], nxt[a], tmp);
        }
        for (a = 0; a <= t; ++a) mpq_set(cur[a], nxt[a]);
    }
    a0 = n / 2;
    factorial(fac, n);
    if ((n & 1) == 0) mpz_mul_ui(fac, fac, 2);
    mpq_set_z(tmp, fac);
    mpq_div(out, cur[a0], tmp);
    mpq_canonicalize(out);
    mpq_clear(tmp); mpz_clears(fac, bin, NULL);
    for (a = 0; a <= n; ++a) { mpq_clear(cur[a]); mpq_clear(nxt[a]); }
    free(cur); free(nxt);
    fprintf(stderr, "bipartite recurrence (6): n=%d, K_%d,%d\n", n, a0, n - a0);
    return 1;
}

static int graph6_decode(const char *s, int *n, uint64_t *adj) {
    int i, j, k, N, maxe, bit;
    unsigned char c;
    if (!s || !*s) return 0;
    c = (unsigned char)s[0];
    if (c == 126) return 0; /* n >= 63 not needed */
    N = (int)c - 63;
    if (N < 1 || N > MAX_SUBSET_N) return 0;
    *n = N;
    for (i = 0; i < N; ++i) adj[i] = 0;
    maxe = N * (N - 1) / 2;
    bit = 0;
    s++;
    k = 0;
    for (j = 1; j < N; ++j) {
        for (i = 0; i < j; ++i) {
            if (bit == 0) {
                if (!s[k]) return 0;
                c = (unsigned char)s[k++] - 63;
                bit = 6;
            }
            bit--;
            if (c & (1u << bit)) {
                adj[i] |= UINT64_C(1) << j;
                adj[j] |= UINT64_C(1) << i;
            }
        }
    }
    (void)maxe;
    return 1;
}

int main(int argc, char **argv) {
    int a,b,n=0,i,m[5],ok=0; uint64_t adj[MAX_SUBSET_N]; mpz_t aut; mpq_t ans;
    if(argc<2) {
        fprintf(stderr,
            "usage: %s bipartite a b | kbal n | c5blowup m | c5parts a b c d e | graph6 <g6> <aut>\n",
            argv[0]);
        return 2;
    }
    mpz_init(aut); mpq_init(ans);
    if(!strcmp(argv[1],"bipartite") && argc==4 && number(argv[2],&a) && number(argv[3],&b))
        ok=make_bipartite(a,b,&n,adj,aut) && subset_eval(n,adj,aut,ans);
    else if(!strcmp(argv[1],"kbal") && argc==3 && number(argv[2],&n))
        ok=kbal_eval(n,ans);
    else if(!strcmp(argv[1],"c5blowup") && argc==3 && number(argv[2],&a)) {
        for(i=0;i<5;++i) m[i]=a; ok=c5_eval(m,ans); n=5*a;
    } else if(!strcmp(argv[1],"c5parts") && argc==7) {
        ok=1; for(i=0;i<5;++i) if(!number(argv[i+2],&m[i])) ok=0;
        if(ok) { for(i=0,n=0;i<5;++i) n+=m[i]; ok=c5_eval(m,ans); }
    } else if(!strcmp(argv[1],"graph6") && argc==4 && number(argv[3],&a)) {
        ok=graph6_decode(argv[2],&n,adj);
        if(ok) { mpz_set_ui(aut,(unsigned long)a); ok=subset_eval(n,adj,aut,ans); }
    }
    if(!ok) { fprintf(stderr,"invalid parameters, allocation failure, or configured limit exceeded\n"); mpq_clear(ans); mpz_clear(aut); return 1; }
    gmp_printf("n=%d likelihood=%Qd\n",n,ans);
    mpq_clear(ans); mpz_clear(aut); return 0;
}
