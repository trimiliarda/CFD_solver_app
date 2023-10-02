
BUILD=../build/
# your_command_string="ls $(SRC)"
# output=$(eval "$your_command_string")

all: *.h
	g++ Polygon.cpp Cell.cpp Mesh.cpp Functions.cpp main.cpp -o $(BUILD)main
	./$(BUILD)main
	rm $(BUILD)main
# all: *.h
# 	gcc *.cpp  -o $(BUILD)main
# 	./$(BUILD)main
# 	rm $(BUILD)main

build:
	g++ *.cpp  -o $(BUILD)main


#======================================== Style    ====================================================================
style:
	clang-format -i $(SRC)*.cpp $(SRC)*.h
	clang-format -n $(SRC)*.cpp $(SRC)*.h
