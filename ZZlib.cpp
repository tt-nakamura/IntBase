// uses NTL
//   www.shoup.net/ntl

#include<NTL/ZZ.h>
#include<NTL/pair.h>
using namespace NTL;

void factor(Vec<Pair<long, long> >& f, long n)
// input:
//   n = integer
// output:
//   f = prime factorization of |n| by trial division
//       vector of (prime, exponent) pair
//       in increasing order of primes
{
    long i(0),j,m,p;

    if(n<0) n=-n;
    if(n==0 || n==1) {
        f.SetLength(0);
        return;
    }
    m = SqrRoot(n);
    PrimeSeq ps;
    while((p = ps.next()) <= m) {
        if(p==0) Error("too large n");
        for(j=0; n%p==0; j++) n/=p;
        if(j==0) continue;
        f.SetLength(i+1);
        f[i].a = p;
        f[i].b = j;
        i++;
        m = SqrRoot(n);
    }
    if(n==1) return;
    f.SetLength(i+1);
    f[i].a = n;
    f[i].b = 1;
}

long conductor(long D)
// input:
//   D = discriminant, D=0 or 1 (mod 4)
// return:
//   largest f such that f^2 divides D or D/4
//   D/4 if s odd or t=3 (mod 4), where D = 2^s t (t odd)
//   D otherwise
{
    long f(1),i,j;
    Vec<Pair<long, long> > p;
    factor(p,D);
    for(i=0, j=D; (j&1)==0; i++) j>>=1;
    if(i&1 || j&2) p[0].b -= 2;
    for(i=0; i<p.length(); i++)
        if(p[i].b >= 2)
            f *= power_long(p[i].a, p[i].b>>1);
    return f;
}

long Jacobi(long a, long b)
// input:
//   a,b = integers, 0 <= a < b, b odd
// return:
//   Jacobi symbol (a/b)
{
    long c,t(0);
    while(a) {
        for(c=0; (a&1)==0; c++) a>>=1;
        if(c&1) {
            c = b&7;
            if(c==3 || c==5) t ^= 1;
        }
        if(a&b&2) t ^= 1;
        c = b%a;
        b = a;
        a = c;
    }
    if(b!=1) return 0;
    else if(t) return -1;
    else return 1;
}

long SqrRootMod(long a, long p)
// input:
//   a = integer, 0 <= a < p
//   p = odd prime
// return:
//   x such that x^2 = a (mod p)
//   by Cipolla method
{
    if(a==0) return 0;
    long d, s((p+1)>>1);
    long b,c0,c1,x0,x1;
    do {
        b = rand()%p;
        d = (b*b-a)%p;
        if(d<0) d+=p;
    } while(Jacobi(d,p) >= 0);
    c0=1; c1=0; x0=b; x1=1;
    for(;;) {
        if(s&1) {
            b = (c0*x0 + c1*x1%p*d)%p;
            c1 = (c0*x1 + c1*x0)%p;
            c0 = b;
        }
        s>>=1;
        if(s==0) return b;
        b = (x0*x0 + x1*x1%p*d)%p;
        x1 = (x0*x1<<1)%p;
        x0 = b;
    }
}

long RootInt(long n, long k)
// input:
//   n,k = positive integers
// return:
//   floor(n^{1/k})
{
    long x, y(long(round(pow(n,1./k))));
    do {
        x = y;
        y = ((k-1)*x + n/power_long(x,k-1))/k;
    } while(x>y);
    return x;
}

long IsPrimePower(ZZ& p, const ZZ& n, long NumTrials)
// input:
//   n = odd integer, n>=3
//   NumTrials = number of times to perform Miller-Rabin test
// output:
//   p = prime factor of n if n = p^k (k>=1)
//       with error probability < 4^{-NumTrials}
// return:
//   k if n = p^k (k>=1)
//   0 otherwise
{
    long k(0);
    ZZ a,d;
    for(p=n; ProbPrime(p, NumTrials)==0; p=d) {
        do RandomBnd(a,p);
        while(MillerWitness(p,a)==0);
        PowerMod(d,a,p,p);
        d -= a;
        GCD(d,d,p);
        if(IsOne(d) || d==p) return 0;
    }
    for(d=n; divide(d,d,p); k++);
    return (IsOne(d) ? k:0);
}
