# Makefile for Fortran Project in hshm directory

# Compiler
FC = gfortran

# Compiler flags
FFLAGS = #-Wall -Wextra -pedantic -std=f2008
DEBUG_FLAGS = -g -fcheck=all -fbacktrace
OPTIMIZE_FLAGS = -O3

# Directories
SRC_DIR = .
ROOT_DIR = ..
BUILD_DIR = build
BIN_DIR = bin
MOD_DIR = $(BUILD_DIR)

# Source files - 添加存在性检查
ROOT_SRC = $(wildcard $(ROOT_DIR)/gdrtesv0.f90)
HSHM_SRC = $(SRC_DIR)/hshm.f90
EXAMPLE_SRC = $(SRC_DIR)/example.f90

# 检查根目录源文件是否存在
ifeq ($(ROOT_SRC),)
$(error gdrtesv0.f90 not found in $(ROOT_DIR))
endif

# Object files
ROOT_OBJ = $(BUILD_DIR)/gdrtesv0.o
HSHM_OBJ = $(BUILD_DIR)/hshm.o
EXAMPLE_OBJ = $(BUILD_DIR)/example.o

# Module files
ROOT_MOD = $(MOD_DIR)/gdrtesv0.mod
HSHM_MOD = $(MOD_DIR)/hshmv4.mod

# Executable
TARGET = $(BIN_DIR)/example

# Default target
all: $(TARGET)

# Create directories
$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

# Compile root module (gdrtesv0.f90)
$(ROOT_OBJ) $(ROOT_MOD): $(ROOT_SRC) | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(OPTIMIZE_FLAGS) -cpp -J$(MOD_DIR) -c $< -o $(ROOT_OBJ)

# Compile hshm module (depends on gdrtesv0.mod)
$(HSHM_OBJ) $(HSHM_MOD): $(HSHM_SRC) $(ROOT_MOD) | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(OPTIMIZE_FLAGS) -cpp -J$(MOD_DIR) -c $< -o $(HSHM_OBJ)

# Compile example program (depends on both modules)
$(EXAMPLE_OBJ): $(EXAMPLE_SRC) $(ROOT_MOD) $(HSHM_MOD) | $(BUILD_DIR)
	$(FC) $(FFLAGS) $(OPTIMIZE_FLAGS) -cpp -J$(MOD_DIR) -c $< -o $(EXAMPLE_OBJ)

# Link all objects to create executable
$(TARGET): $(ROOT_OBJ) $(HSHM_OBJ) $(EXAMPLE_OBJ) | $(BIN_DIR)
	$(FC) $(FFLAGS) $(OPTIMIZE_FLAGS) -I$(MOD_DIR) $^ -o $@

# Debug build
debug: FFLAGS += $(DEBUG_FLAGS)
debug: $(TARGET)

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Run the program
run: $(TARGET)
	./$(TARGET)

# Show help
help:
	@echo "Available targets:"
	@echo "  all      - Build the program (default)"
	@echo "  debug    - Build with debug flags"
	@echo "  clean    - Remove build artifacts"
	@echo "  run      - Build and run the program"
	@echo "  help     - Show this help message"

# Phony targets
.PHONY: all debug clean run help

# 添加对模块依赖关系的处理
$(HSHM_OBJ): $(ROOT_MOD)
$(EXAMPLE_OBJ): $(ROOT_MOD) $(HSHM_MOD)