
BUILD=../build/
LIBS= -lm -pthread
# your_command_string="ls $(SRC)"
# output=$(eval "$your_command_string")

all: *.h build
	./$(BUILD)main
	rm $(BUILD)main

build:
	g++ *.cpp  -o $(BUILD)main


#========================================   Style   ====================================================================
style:
	clang-format -i $(SRC)*.cpp $(SRC)*.h
	clang-format -n $(SRC)*.cpp $(SRC)*.h
