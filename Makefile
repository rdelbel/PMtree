FLAGS = -O3 -std=c++11
OBJS = main.o pmtree.o
CXX = g++

all: $(OBJS)
	$(CXX) $(FLAGS) -o dtree $(OBJS) -lrt

main.o: main.cpp
	$(CXX) -c $(FLAGS) main.cpp

pmtree.o: pmtree.cpp
	$(CXX) -c $(FLAGS) pmtree.cpp
clean: 
	rm *.o
