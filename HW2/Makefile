CXX = mpicxx
CXXFLAGS = -std=c++11 

# clear any suffixes
.SUFFIXES:
# set suffixes
.SUFFIXES: .cpp .o

all: hw2_DoughertyP 
	@echo Patrick Dougherty\'s HW2 program has been compiled

hw1_DoughertyP: hw2_DoughertyP.o
	$(CXX) $(CXXFLAGS) $? -o $@

# create object file (*.o) from source file (*.cpp)
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<


.PHONY: clean
clean:
	$(RM) -r hw2_DoughertyP *.o *.dSYM
