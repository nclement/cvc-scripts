CXX = g++
#FLAGS = -g
FLAGS = -O3 -Wall

all : conformationGenerator converter clashDetector conformationGenerator_single

conformationGenerator : conformationGenerator.cpp
	${CXX} ${FLAGS} -o $@ $<

converter : converter.cpp
	${CXX} ${FLAGS} -o $@ $<

clashDetector : clashDetector.cpp
	${CXX} ${FLAGS} -o $@ $<

% : %.cpp
	${CXX} ${FLAGS} -o $@ $<

