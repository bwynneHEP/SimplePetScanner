BUILD_DIR=./build

# Compile, set up build dir if not done already
all: $(BUILD_DIR)/Makefile
	make -C $(BUILD_DIR)

# Set up the build dir
$(BUILD_DIR)/Makefile:
	make -C $(BUILD_DIR) -f initial

# Tidy for rebuild
clean:
	make -C $(BUILD_DIR) clean
	rm -rf build/CMakeFiles/
	rm -f build/cmake_install.cmake
	rm -f build/Makefile
	rm -f build/CMakeCache.txt
