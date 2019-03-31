// uses NTL
//   www.shoup.net/ntl

#include<NTL/mat_lzz_p.h>
using namespace NTL;

void SuplBase(mat_zz_p& B, const mat_zz_p& M)
// input:
//   M = m by n matrix of rank m (m<n)
// output:
//   B = n by n invertible matrix
//       first m rows of B are those of M
{
    long i,j,k,l,m(M.NumRows()), n(M.NumCols());
    zz_p a,b;
    mat_zz_p A(M),I;
    ident(I,n);
    for(k=0; k<m; k++) {
        for(l=k; l<n; l++)
            if(!IsZero(A[k][l])) break;
        if(l==n) Error("rows of M are dependent");
        if(l>k) {
            swap(I[l], I[k]);
            swap(A[k][l], A[k][k]);
        }
        inv(a, A[k][k]);
        for(i=k+1; i<m; i++) {
            if(l>k) swap(A[i][k], A[i][l]);
            if(IsZero(A[i][k])) continue;
            A[i][k] *= a;
            for(j=0; j<n; j++) {
                if(j==k || IsZero(A[k][j])) continue;
                mul(b, A[i][k], A[k][j]);
                A[i][j] -= b;
            }
        }
    }
    if(&B!=&M) B=M;
    B.SetDims(n,n);
    for(i=m; i<n; i++) B[i] = I[i];
}
