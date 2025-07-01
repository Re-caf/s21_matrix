NAME = s21_matrix
OPEN_CMD = xdg-open
GCC = gcc
CC_FLAGS = -Wall -Werror -Wextra -std=c11
RFLAGS = -lcheck -lpthread -lm -lsubunit
SOURCE = $(wildcard s21_*.c)
TESTS_SOURCE = tests.c

all: clean s21_matrix.a

s21_matrix.a: s21_matrix.o
	ar rcs $(NAME).a *.o
	ranlib $(NAME).a
	@rm -rf *.o

s21_matrix.o: 
	$(GCC) $(CC_FLAGS) -c s21_matrix.c

test: ${TESTS_SOURCE} s21_matrix.a
	make clang
	${CC} $(CC_FLAGS) -o test $^ $(RFLAGS)
	./test

gcov_report: clean s21_matrix.a
	gcc ${CC_FLAGS} --coverage ${SOURCE} tests.c -o s21_test ${RFLAGS}
	./s21_test
	lcov -t "s21_test" -o s21_test.info -c -d .
	genhtml -o report s21_test.info
	$(OPEN_CMD) ./report/index.html


clang:
	cp ../materials/linters/.clang-format .
	clang-format -i *.c *.h
	rm -rf .clang-format

clean:
	-rm -rf *.o && rm -rf *.gcno
	-rm -rf *.a && rm -rf *.gcda
	-rm -rf *.info && rm -rf *.gcov
	-rm -rf ./test && rm -rf ./gcov_report
	-rm -rf ./report/
	-rm -rf s21_test
	-rm -rf *.out

