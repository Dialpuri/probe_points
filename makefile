# nautilus makefile - Jordan Dialpuri Edit 12/10/22

LIBS = \
-lccp4c \
-lclipper-ccp4 \
-lclipper-contrib \
-lclipper-core \
-lclipper-minimol \
-lclipper-mmdb \
-lfftw \
-lmmdb2 \
-lrfftw \
-lm \
-lstdc++

FLAGS = \
-std=c++11 \
-O2 \
-Wall \
-Wno-sign-compare \
-fPIC \
-ftemplate-depth-50 \
-D_GLIBCXX_USE_CXX11_ABI=0 \
-g

CCP4 = "../ccp4-8.0"

# INCDIR = -I${CCP4}/include
INCDIR = -Iinclude
LIBDIR = -L${CCP4}/lib
#LIBDIR = -Llib
CFLAGS = ${FLAGS} ${INCDIR}
LFLAGS = ${FLAGS} ${LIBDIR} ${LIBS}

SRC_DIR = src
BIN_DIR = bin

BUILD_PRINT = @echo "\e[1;34mBuilding $<\e[0m"
COMPLETE_PRINT = @echo "\e[1;32mBuilding complete!\e[0m"

MKDIR_P = mkdir -p

#Source CCP4
#IGNORE := $(shell bash -c "source ../ccp4-8.0/bin/ccp4.setup-sh; env | sed 's/=/:=/' | sed 's/^/export /' > makeenv")
include makeenv

$(shell mkdir -p $(BIN_DIR))


all:
	@g++ src/hash.cpp -o hash_exec ${CFLAGS} ${LFLAGS}

clean:
	rm bin/*.o hash*
