// uses NTL
//   www.shoup.net/ntl

#include<NTL/ZZ.h>
#include<NTL/pair.h>
using namespace NTL;

#define	TRYDIV_BOUND (1<<16)
#define MR_NUM_TRIAL 20
#define RHO_TIME_OUT 5

long IsPrimePower(ZZ&, const ZZ&, long);
long brent_rho(ZZ&, const ZZ&, double);
long mpqs(ZZ&, const ZZ&);

void factor_(Vec<Pair<ZZ, long> >& f, const ZZ& n)
// input:
//   n = odd integer, n>=3
// output:
//   f = prime factorization of n (appended to f)
{
    long i,j,k(f.length());
    ZZ p,q;
    if(j = IsPrimePower(p, n, MR_NUM_TRIAL)) {
        f.SetLength(k+1);
        f[k].a = p;
        f[k].b = j;
        return;
    }
    Vec<Pair<ZZ, long> > g,h;
    if(brent_rho(p, n, RHO_TIME_OUT) == 0);
    else if(mpqs(p,n) == 0);
    else Error("factor not found");
    div(q,n,p);
    factor_(g,p);
    factor_(h,q);
    for(i=j=0; i<g.length() || j<h.length(); k++) {
        f.SetLength(k+1);
        if(j==h.length() || i<g.length() && g[i].a < h[j].a)
            f[k] = g[i++];
        else if(i==g.length() || g[i].a > h[j].a)
            f[k] = h[j++];
        else {
            f[k] = g[i++];
            f[k].b += h[j++].b;
        }
    }
}

void factor(Vec<Pair<ZZ, long> >& f, const ZZ& n)
// input:
//   n = integer
// output:
//   f = prime factorization of |n|
//       vector of (prime, exponent) pair
//       in increasing order of primes
{
    long i(0),j,p;
    ZZ m;
    abs(m,n);
    if(IsZero(m) || IsOne(m)) {
        f.SetLength(0);
        return;
    }
    if(j = MakeOdd(m)) {
        f.SetLength(1);
        f[0].a = 2;
        f[0].b = j;
        if(IsOne(m)) return;
        i++;
    }
    PrimeSeq ps;
    ps.reset(3);
    while((p = ps.next()) <= TRYDIV_BOUND) {
        for(j=0; divide(m,m,p); j++);
        if(j==0) continue;
        f.SetLength(i+1);
        f[i].a = p;
        f[i].b = j;
        if(IsOne(m)) return;
        i++;
    }
    factor_(f,m);
}
