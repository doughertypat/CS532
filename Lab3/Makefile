CXX = mpicxx
CXXFLAGS = -std=c++11 

# clear any suffixes
.SUFFIXES:
# set suffixes
.SUFFIXES: .cpp .o

all: lab3 
	@echo The Lab3 program has been compiled

lab3: Lab3_DoughertyP.o
	$(CXX) $(CXXFLAGS) $? -o $@

# create object file (*.o) from source file (*.cpp)
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<


.PHONY: clean
clean:
	$(RM) -r lab3 *.o *.dSYM
