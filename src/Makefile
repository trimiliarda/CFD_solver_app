SRC=modul_cxx/
BUILD=../build/
# your_command_string="ls $(SRC)"
# output=$(eval "$your_command_string")

all:
	g++ $(SRC)*.cpp  -o $(BUILD)main
	./$(BUILD)main
	rm $(BUILD)main

#======================================== Style    ====================================================================
style:
	clang-format -i $(SRC)*.cpp $(SRC)*.h
	clang-format -n $(SRC)*.cpp $(SRC)*.h
