STD=-std=c++17
CXX=g++ $(STD)
DEBUG=-g -O3
CXXFLAGS=-Wall $(DEBUG)
LFLAGS=-Wall $(DEBUG)
DEPS=main.o

main: $(DEPS)
	$(CXX) $(LFLAGS) $^ -o $@

%.o: %.c
	$(CXX) $(CXXFLAGS)

clean:
	rm *.o 
