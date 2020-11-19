OBJS = testomp.o

CXX 	= g++
COMMON 	= -std=c++11 -fopenmp
CFLAGS	= $(COMMON) -c
LDFLAGS	= $(COMMON)

all: testomp
	@echo TestOMP has been compiled

testomp: $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o testomp
	
testomp.o: testomp.cpp
	$(CXX) $(CFLAGS) testomp.cpp

.PHONY: clean
clean: 
	rm -f testomp *.o
