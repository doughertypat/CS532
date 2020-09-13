CXX = mpicxx
CXXFLAGS = -std=c++11 

# clear any suffixes
.SUFFIXES:
# set suffixes
.SUFFIXES: .cpp .o

all: hw1_DoughertyP 
	@echo Patrick Dougherty\'s HW1 program has been compiled

hw1_DoughertyP: hw1_DoughertyP.o
	$(CXX) $(CXXFLAGS) $? -o $@

# create object file (*.o) from source file (*.cpp)
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<


.PHONY: clean
clean:
	$(RM) -r hw1_DoughertyP *.o *.dSYM
