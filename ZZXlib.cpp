// uses NTL
//   www.shoup.net/ntl

#include<NTL/ZZXFactoring.h>
using namespace NTL;

void PowerMod(ZZX& x, const ZZX& a, long n, const ZZX& f)
// input:
//   a = polynomial
//   n = integer, n>=0
//   f = monic polynomial
// output:
//   x = a^n mod f
{
    if(&x==&a) { ZZX b(a); PowerMod(x,b,n,f); return; }
    if(n==0 || IsOne(a)) { set(x); return; }
    long m(1L<<(NumBits(n)-1));
    x=a;
    for(m>>=1; m; m>>=1) {
        SqrMod(x,x,f);
        if(n&m) MulMod(x,x,a,f);
    }
}

long IsIrred(const ZZX& f)
// input:
//   f = polynomial
// return:
//   1 if f irreducible, 0 otherwise
{
    ZZ c;
    Vec<Pair<ZZX, long> > h;
    factor(c,h,f);
    return deg(h[0].a) == deg(f);
}

void BuildIrred(ZZX& f, long n, long b)
// input:
//   n,b = positive integers
// output:
//   f = random irreducible monic polynomial of degree n
//       | coefficients of f | < b
{
    long i;
    clear(f);
    SetCoeff(f,n);
    do for(i=0; i<n; i++) {
        f[i] = RandomBnd(b);
        if(RandomBits_long(1)) negate(f[i], f[i]);
    } while(!IsIrred(f));
}
