BUILD_DIR ?= build
CONFIG ?= Release

.PHONY: all configure build clean

all: build

configure:
	cmake -S . -B $(BUILD_DIR) -DCMAKE_BUILD_TYPE=$(CONFIG)

build: configure
	cmake --build $(BUILD_DIR) --config $(CONFIG)

clean:
	cmake -E rm -rf $(BUILD_DIR)
