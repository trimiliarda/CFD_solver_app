
BUILD=../build/
LIBS= -lm -pthread

FILES=
#======================================== 
UNAME:=$(shell uname -s)
ifeq ($(UNAME),Linux)
	FILES=*.cpp
endif
ifeq ($(UNAME),Darwin)
	FILES=Polygon.cpp Cell.cpp Mesh.cpp Functions.cpp main.cpp
endif

# your_command_string="ls $(SRC)"
# output=$(eval "$your_command_string")

all: build
	./$(BUILD)main
	rm $(BUILD)main

build:
	g++ $(FILES) -o $(BUILD)main


#========================================   Style   ====================================================================
style:
	clang-format -i $(SRC)*.cpp $(SRC)*.h
	clang-format -n $(SRC)*.cpp $(SRC)*.h
