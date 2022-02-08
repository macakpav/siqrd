# Makefile
#-----------------------------------------------------------------------------------------
# Variable for the selected compiler
CC=g++# g++ clang++

# set to false to debug, true to optimize
optimize=true

# folder system
bin_folder=
obj_folder=obj/
src_folder=src/

# Possible debug flags per compiler
# Additional debug options in src/debug_levels.hpp
CFLAGS_g++ = -std=c++17 -g -Wall -Wextra -Werror -DDLVL2 #-DDMETHODS
CFLAGS_clang++ = -std=c++17 -g -Wall -Wextra -Werror -DDLVL2

# Optimization flags
CFLAGS_g++_opt = -std=c++17 -O3 -DNDEBUG
CFLAGS_clang++_opt = -std=c++17 -O3 -DNDEBUG

# Select the right flags for the current compiler
ifeq ($(optimize) , true)
	CFLAGS=$(CFLAGS_$(CC)_opt)
else
	CFLAGS=$(CFLAGS_$(CC))
endif

# Items to remove with make clean
r=$(obj_folder)*.o $(bin_folder)*.exe

#-----------------------------------------------------------------------------------------
default:
	@ echo "Only make all, make allrun, make run1,2,3,4 and make clean are possible."
	@ echo "	 (For benchmarking: make time, mem and prof.)"

all: clean simulation solvertest estimation1 estimation2
allrun: all run1 run2 run3 run4
clean:
	@ rm -f $(r)
	@ clear

./$(obj_folder)simulation.o: ./$(src_folder)simulation.cpp ./$(src_folder)debug_levels.hpp ./$(src_folder)siqrd/*.hpp ./$(src_folder)ode/*.hpp ./$(src_folder)saving/*.hpp
	$(CC) -c  $(CFLAGS) ./$(src_folder)simulation.cpp -o ./$(obj_folder)simulation.o

simulation: ./$(obj_folder)simulation.o
	$(CC) -o ./$(bin_folder)simulation.exe ./$(obj_folder)simulation.o

run1: simulation
	./$(bin_folder)simulation.exe 100 100

./$(obj_folder)solvertest.o: ./$(src_folder)solvertest.cpp ./$(src_folder)debug_levels.hpp ./$(src_folder)siqrd/*.hpp ./$(src_folder)ode/*.hpp ./$(src_folder)saving/*.hpp
	$(CC) -c  $(CFLAGS) ./$(src_folder)solvertest.cpp -o ./$(obj_folder)solvertest.o

solvertest: ./$(obj_folder)solvertest.o
	$(CC) -o ./$(bin_folder)solvertest.exe ./$(obj_folder)solvertest.o

run2: solvertest
	./$(bin_folder)solvertest.exe 50000 500

./$(obj_folder)estimation1.o: ./$(src_folder)estimation1.cpp ./$(src_folder)debug_levels.hpp ./$(src_folder)siqrd/*.hpp ./$(src_folder)ode/*.hpp ./$(src_folder)saving/*.hpp ./$(src_folder)optimization/*.hpp
	$(CC) -c  $(CFLAGS) $(src_folder)estimation1.cpp -o ./$(obj_folder)estimation1.o

estimation1: ./$(obj_folder)estimation1.o
	$(CC) -o ./$(bin_folder)estimation1.exe ./$(obj_folder)estimation1.o

run3: estimation1
	./$(bin_folder)estimation1.exe

./$(obj_folder)bench_time.o: ./$(src_folder)bench_time.cpp ./$(src_folder)debug_levels.hpp ./$(src_folder)siqrd/*.hpp ./$(src_folder)ode/*.hpp ./$(src_folder)saving/*.hpp ./$(src_folder)optimization/*.hpp
	$(CC) -c  $(CFLAGS_$(CC)_opt) -DNINFO ./$(src_folder)bench_time.cpp -o ./$(obj_folder)bench_time.o

bench_time: ./$(obj_folder)bench_time.o
	$(CC) -o ./$(bin_folder)bench_time.exe ./$(obj_folder)bench_time.o

time: clean bench_time
	./$(bin_folder)bench_time.exe


./$(obj_folder)bench_mem.o: ./$(src_folder)bench_mem.cpp ./$(src_folder)debug_levels.hpp ./$(src_folder)siqrd/*.hpp ./$(src_folder)ode/*.hpp ./$(src_folder)saving/*.hpp ./$(src_folder)optimization/*.hpp
	g++ -c -std=c++17 -Wall -ggdb3 -DNDEBUG ./$(src_folder)bench_mem.cpp -o ./$(obj_folder)bench_mem.o

bench_mem: ./$(obj_folder)bench_mem.o
	g++ -pg -o ./$(bin_folder)bench_mem.exe ./$(obj_folder)bench_mem.o

mem: clean bench_mem
	@valgrind ./$(bin_folder)bench_mem.exe

prof: clean bench_mem
	./$(bin_folder)bench_mem.exe
	gprof bench_mem.exe > analysis.txt

./$(obj_folder)estimation2.o: ./$(src_folder)estimation2.cpp ./$(src_folder)debug_levels.hpp ./$(src_folder)siqrd/*.hpp ./$(src_folder)ode/*.hpp ./$(src_folder)saving/*.hpp ./$(src_folder)optimization/*.hpp
	$(CC) -c  $(CFLAGS) ./$(src_folder)estimation2.cpp -o ./$(obj_folder)estimation2.o

estimation2: ./$(obj_folder)estimation2.o
	$(CC) -o ./$(bin_folder)estimation2.exe ./$(obj_folder)estimation2.o

run4: estimation2
	./$(bin_folder)estimation2.exe

# pdf: plot

# plot:
# 	@ pdflatex plot.tex > /dev/null


