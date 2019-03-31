NTL = -lntl -lgmp
OBJ = round2.o SuplBase.o ZZFactoring.o rho.o mpqs.o ZZlib.o ZZXlib.o

table1: table1.o IntBase.o $(OBJ)
	g++ table1.o IntBase.o $(OBJ) $(NTL)

table2: table2.o $(OBJ)
	g++ table2.o $(OBJ) $(NTL)
