# 1. 编译器基本设置
FC = gfortran
BUILD_DIR = build
BIN_DIR = bin

# 基础编译参数（指定 module 存放目录）
FFLAGS = -cpp -J$(BUILD_DIR)

# 2. 区分 Release 和 Debug 的编译参数
# 严格的调试参数：包含代码行回溯、全边界检查、禁用优化、开启警告
DEBUG_FLAGS = -g -O0 -fcheck=all -fbacktrace -Wall -Wextra -pedantic
# 优化参数
OPTIMIZE_FLAGS = -O3

# 3. 文件路径定义
ROOT_SRC = ../gdrtes.f90
SJHM_SRC = sjhm.f90
SJHMTEST_SRC  = sjhmtest.f90

ROOT_OBJ = $(BUILD_DIR)/gdrtes.o
SJHM_OBJ = $(BUILD_DIR)/sjhm.o
SJHMTEST_OBJ  = $(BUILD_DIR)/sjhmtest.o

TARGET   = $(BIN_DIR)/sjhmtest

# 4. 核心编译规则
# 默认编译为 Release 版本（使用优化）
all: FFLAGS += $(OPTIMIZE_FLAGS)
all: $(TARGET)

# 调试编译目标：追加严格检查参数
debug: FFLAGS += $(DEBUG_FLAGS)
debug: $(TARGET)

# 创建必要的目录
$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

# 编译第一个模块 (gdrtes)
$(ROOT_OBJ): $(ROOT_SRC) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# 编译第二个模块 (sjhm) - 依赖 ROOT_OBJ
$(SJHM_OBJ): $(SJHM_SRC) $(ROOT_OBJ) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# 编译主程序 (sjhmtest) - 依赖前两个 OBJ
$(SJHMTEST_OBJ): $(SJHMTEST_SRC) $(ROOT_OBJ) $(SJHM_OBJ) | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# 链接生成可执行文件
$(TARGET): $(ROOT_OBJ) $(SJHM_OBJ) $(SJHMTEST_OBJ) | $(BIN_DIR)
	$(FC) $(FFLAGS) -I$(BUILD_DIR) $^ -o $@

# 5. 清理
clean:
	rm -f $(ROOT_OBJ) $(SJHM_OBJ) $(SJHMTEST_OBJ)
	rm -f $(ROOT_MOD) $(SJHM_MOD)
	rm -f $(TARGET)

.PHONY: all debug clean
