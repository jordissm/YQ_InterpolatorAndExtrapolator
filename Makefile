# Makefile for compiling main.cpp with different compilers and flags for macOS and Linux

# Define the source directory, source file, and the output executable
SRC_DIR = src
SRC = $(SRC_DIR)/main.cpp
EXE = YQ_InterpolatorAndExtrapolator

# Compiler and flags for macOS
MACOS_COMPILER = clang++
MACOS_FLAGS = -Xclang -fopenmp -std=c++17 -O3 -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp

# Compiler and flags for Linux
LINUX_COMPILER = g++
LINUX_FLAGS = -fopenmp -std=c++17

# Determine the operating system
UNAME_S := $(shell uname -s)

# Default target
all: $(EXE)

# Compilation rule
$(EXE): $(SRC)
ifeq ($(UNAME_S), Darwin)
	$(MACOS_COMPILER) $(MACOS_FLAGS) $(SRC) -o $(EXE)
else ifeq ($(UNAME_S), Linux)
	$(LINUX_COMPILER) $(LINUX_FLAGS) $(SRC) -o $(EXE)
else
	@echo "Unsupported OS: $(UNAME_S)"
	exit 1
endif

# Clean rule
clean:
	rm -f $(EXE)
