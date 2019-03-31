// uses NTL
//   www.shoup.net/ntl

#include<NTL/ZZX.h>
#include<NTL/pair.h>
using namespace NTL;

long conductor(long);
long RootInt(long, long);

long IntBase(Vec<Pair<ZZX,long> >& w, const ZZX& T)
// input:
//   T = monic irreducible polynomial defining number field K
// output:
//   w = integral bases of K
//       vector of (monic polynomial, denominator) pair
//       in Hermite normal form
// return:
//   discriminant of K
{
    long i,j,D,F,n(deg(T));
    ZZ a;
    ZZX f;

    discriminant(a,T);
    if(!a.SinglePrecision())
        Error("too large discriminant");
    conv(D,a);
    F = conductor(D);
    D /= F*F;
    w.SetLength(n);
    set(w[0].a);
    w[0].b = 1;
    for(i=1; i<n; i++) {
        clear(w[i].a);
        SetCoeff(w[i].a, i);
        for(w[i].b = RootInt(F, n-i); w[i].b >= w[i-1].b; w[i].b--) {
            if(w[i].b == 1) break;
            if(F%power_long(w[i].b, n-i)) continue;
            do {
                CharPolyMod(f, w[i].a, T);
                for(j=n-1, a=w[i].b; j>=0; j--, a*=w[i].b)
                    if(!divide(f[j], a)) break;
                if(j<0) break;
                for(j=0; j<i; j++)
                    if(++(w[i].a[j]) < w[i].b/w[j].b) break;
                    else clear(w[i].a[j]);
            } while(j<i);
            if(j<0) break;
        }
        F /= w[i].b;
    }
    return D*F*F;
}
