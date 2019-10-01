CC=/usr/bin/g++-9

proj1:
	$(CC) P1.cpp -o P1 -lpthread -std=c++17 -gdwarf -lboost_system -Og
