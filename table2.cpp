// uses NTL
//   www.shoup.net/ntl

#include<NTL/ZZX.h>
#include<NTL/pair.h>
using namespace NTL;

void round2(ZZ&, Vec<Pair<ZZX,long> >&, const ZZX&);
long IsIrred(const ZZX&);

main() {
    long n,b(32);
    ZZX T;
    ZZ D;
    Vec<Pair<ZZX,long> > w;
    for(n=4; n<=15; n++) {
        do {
            clear(T);
            SetCoeff(T, n);
            SetCoeff(T, 0, RandomBnd(2*b-1)-b+1);
            SetCoeff(T, RandomBnd(n-1)+1, RandomBnd(2*b-1)-b+1);
        } while(!IsIrred(T));
        round2(D,w,T);
        std::cout << T << ' ';
        std::cout << D << ' ';
        std::cout << w << '\n';
    }
}
