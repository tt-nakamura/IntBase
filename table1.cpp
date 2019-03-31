// uses NTL
//   www.shoup.net/ntl

#include<NTL/ZZX.h>
#include<NTL/pair.h>
using namespace NTL;

long IntBase(Vec<Pair<ZZX,long> >&, const ZZX&);
void round2(ZZ&, Vec<Pair<ZZX,long> >&, const ZZX&);
void BuildIrred(ZZX&, long, long);

main() {
    long n(3),b(32),N(10),D,i;
    ZZX T;
    ZZ D1;
    Vec<Pair<ZZX,long> > w,w1;
    for(i=0; i<N; i++) {
        BuildIrred(T,n,b);
        D = IntBase(w,T);
        round2(D1,w1,T);
        if(D!=D1 || w!=w1) Error("wrong");
        std::cout << T << ' ';
        std::cout << D << ' ';
        std::cout << w << '\n';
    }
}
