biginteger: biginteger.out
	@./biginteger.out

biginteger.out: biginteger.cpp biginteger.h
	@g++ -std=c++20 -o biginteger.out biginteger.cpp -Wall -Wextra -Werror -fsanitize=address
