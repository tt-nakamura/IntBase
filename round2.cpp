// uses NTL
//   www.shoup.net/ntl

#include<NTL/ZZX.h>
#include<NTL/HNF.h>
#include<NTL/lzz_pXFactoring.h>
using namespace NTL;

void factor(Vec<Pair<ZZ, long> >&, const ZZ&);
void PowerMod(ZZX&, const ZZX&, long, const ZZX&);
void SuplBase(mat_zz_p&, const mat_zz_p&);

void round2(ZZ& D, mat_ZZ& W, ZZ& d, const ZZX& T)
// input:
//   T = monic irreducible polynomial defining number field K
// output:
//   D = discriminant of K
//   W/d = integral bases of K in Hermite normal form
//         by Zassenhaus round2 method
// reference:
//   H. Cohen "A Course in Computational Algebraic Number Theory"
//     (Springer) section 6.1
{
    long i,j,k,l,p,q,m,n(deg(T));
    ZZ a,b;
    zz_p c;
    ZZX f,f1,w[n],u[n];
    zz_pX g,g1,h,Tp;
    mat_ZZ A,B;
    mat_zz_p V;
    Vec<Pair<ZZ, long> > p1;
    Vec<Pair<zz_pX, long> > t;

    set(d);
    ident(W,n);
    for(i=0; i<n; i++) w[i].SetLength(i+1);

    discriminant(D,T);
    factor(p1,D);
    if(MakeOdd(a=D)&1 || a%4==3) p1[0].b -= 2;

    for(l=0; l<p1.length(); l++) {
        if(p1[l].b < 2) continue;
        if(!p1[l].a.SinglePrecision())
            Error("too large discriminant");
        p1[l].b >>= 1;
        conv(p, p1[l].a);
        zz_p::init(p);

        conv(Tp,T);
        CanZass(t,Tp);
        set(g);
        for(i=0; i<t.length(); i++) g *= t[i].a;
        div(h,Tp,g);
        conv(f,g);
        conv(f1,h);
        f *= f1;
        f -= T;
        f /= p;
        conv(g1,f);
        GCD(g1,g1,g);
        GCD(g1,g1,h);// Z
        m = deg(g1);
        if(m==0) continue;

        div(g,Tp,g1);
        conv(f,g);// U
        A.SetDims(n+m,n);
        clear(A);
        for(i=0; i<m; i++) {
            for(j=0; j<=i; j++) w[i][j] = W[i][j];
            MulMod(f1,f,w[i],T);
            for(j=0; j<=deg(f1); j++) A[i][j] = f1[j];
        }
        for(i=0; i<n; i++)
            for(j=0; j<=i; j++) mul(A[m+i][j], W[i][j], p);
        mul(b,p,d);
        power(a,b,n);
        HNF(W,A,a);
        if(IsOne(W[n-1][n-1])) d = b;
        else for(i=0; i<n; i++)
            for(j=0; j<=i; j++) W[i][j] /= p;
        
        if(m >= p1[l].b) continue;

        for(q=p; q<n; q*=p);

        for(;;) {
            A.SetDims(n,n);
            clear(A);
            for(i=0; i<n; i++) {
                for(j=0; j<=i; j++) w[i][j] = W[i][j];
                PowerMod(f,w[i],q,T);
                for(j=0; j<=deg(f); j++) A[i][j] = f[j];
            }
            inv(a,B,W);
            power(b,d,q-1);
            A *= B;
            a *= b;
            for(i=0; i<n; i++)
                for(j=0; j<n; j++) A[i][j] /= a;
            conv(V,A);
            kernel(V,V);
            m = V.NumRows();
            SuplBase(V,V);
            conv(A,V);
            A *= W;
            for(i=m; i<n; i++) A[i] *= p;

            for(i=0; i<n; i++)
                for(j=0; j<n; j++) SetCoeff(u[i], j, A[i][j]);
            inv(b,A,A);
            b *= d;
            V.SetDims(n,n*n);
            for(k=0; k<n; k++) {
                clear(B);
                for(i=0; i<n; i++) {
                    MulMod(f,w[i],u[k],T);
                    for(j=0; j<=deg(f); j++) B[i][j] = f[j];
                }
                B *= A;
                for(i=0; i<n; i++) {
                    for(j=0; j<n; j++) {
                        zz_p& v(V[i][k*n+j]);
                        if(divide(a, B[i][j], b))
                            conv(v,a);
                        else {
                            GCD(a, B[i][j], b);
                            B[i][j] /= a;
                            div(a,b,a);
                            conv(v, B[i][j]);
                            conv(c,a);
                            v /= c;
                        }
                    }
                }
            }
            kernel(V,V);
            m = V.NumRows();
            if(m==0) break;

            conv(B,V);
            B *= W;
            B.SetDims(n+m,n);
            for(i=0; i<n; i++)
                for(j=0; j<n; j++) mul(B[m+i][j], W[i][j], p);
            mul(b,p,d);
            power(a,b,n);
            HNF(A,B,a);
            a = d;
            if(IsOne(A[n-1][n-1])) d = b;
            else for(i=0; i<n; i++)
                for(j=0; j<n; j++) A[i][j] /= p;
            if(d==a && A==W) break;
            W = A;
        }
    }
    for(i=0; i<n; i++) {
        div(a,d,W[i][i]);
        sqr(a,a);
        D /= a;
    }
}

void round2(ZZ& D, Vec<Pair<ZZX,long> >& w, const ZZX& T)
// input:
//   T = monic irreducible polynomial defining number field K
// output:
//   D = discriminant of K
//   w = integral bases of K
//       vector of (monic polynomial, denominator) pair
//       in Hermite normal form
{
    long i,j;
    ZZ a,d;
    mat_ZZ W;
    round2(D,W,d,T);
    w.SetLength(W.NumRows());
    for(i=0; i<w.length(); i++) {
        div(a, d, W[i][i]);
        conv(w[i].b, a);
        w[i].a.SetLength(i+1);
        for(j=0; j<=i; j++)
            div(w[i].a[j], W[i][j], W[i][i]);
    }
}
