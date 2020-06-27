CC=g++
EXE=./testgrid.exe


compile:$(EXE)

$(EXE):./testgrid.cpp ../src/lib/grid.h ../src/lib/grid.cpp
	$(CC) -o $(EXE) $(CFLAG) ./testgrid.cpp ../src/lib/grid.cpp


test:$(EXE)
	$(EXE)
