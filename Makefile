FLAGS = -O3
OBJS = main.o
CXX = g++

all: $(OBJS)
	$(CXX) $(FLAGS) -o dtree $(OBJS) -lrt

main.o: main.cpp
	$(CXX) -c $(FLAGS) main.cpp

clean: 
	rm *.o
