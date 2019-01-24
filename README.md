# roadefChallenge2018S21

This is a submission of team S21 to the ROADEF/EURO Challenge 2018 (http://www.roadef.org/challenge/2018/en/index.php)

Hybrid heuristic algorithm is implemented in C++.

Compilation is simple, just use the c++11 standard and include pthread library:

g++ src/roadef2018.cpp solution_generator/make_solutionfile.cpp algo/*.cpp -o challengeSG -std=c++11 -pthread -O3
