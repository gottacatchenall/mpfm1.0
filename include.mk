# ========================================
# Source
# ========================================

ROOT_DIR:=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

MAIN =              $(ROOT_DIR)/src/core/main.cpp

INCLUDE_DIRS = 		$(ROOT_DIR)/src/core/						\
					$(ROOT_DIR)/src/core/AlleleTracker 			\
					$(ROOT_DIR)/src/core/EnvFactor 				\
					$(ROOT_DIR)/src/core/Individual 			\
					$(ROOT_DIR)/src/core/GenomeDict	 			\
					$(ROOT_DIR)/src/core/Patch

SRCS =				$(ROOT_DIR)/src/core/init.cpp 										\
					$(ROOT_DIR)/src/core/AlleleTracker/AlleleTracker.cpp				\
					$(ROOT_DIR)/src/core/EnvFactor/EnvFactor.cpp 						\
					$(ROOT_DIR)/src/core/GenomeDict/GenomeDict.cpp						\
					$(ROOT_DIR)/src/core/Individual/Individual.cpp						\
					$(ROOT_DIR)/src/core/Patch/Patch.cpp								\
					$(ROOT_DIR)/src/core/random_gen.cpp



# ========================================
# Tests
# ========================================
TEST_MAIN =         $(ROOT_DIR)/test/main.cc
TEST_SRCS = 		$(ROOT_DIR)/test/tests.cc

GTEST_LIB = /usr/local/lib/libgtest.a
GMOCK_LIB = /usr/local/lib/libgmock.a
