FLAGS = \
-std=c++11 \
-O2 \
-Wall \
-Wno-sign-compare \
-fPIC \
-ftemplate-depth-50 \
-D_GLIBCXX_USE_CXX11_ABI=0

LIBS = \
-l:libccp4c.so.8.0 \
-l:libclipper-ccp4.so.2 \
-l:libclipper-contrib.so.2 \
-l:libclipper-core.so.2 \
-l:libclipper-minimol.so.2 \
-l:libclipper-mmdb.so.2 \
-l:libfftw.so.2 \
-l:libmmdb2.so.0 \
-l:librfftw.so.2 \
-lm \
-lstdc++ \

SHARED = \
	probe


OBJS = $(SHARED:=.o)

INCDIR = -Iinclude
LIBDIR = -L${CCP4}/lib

CFLAGS = ${FLAGS} ${INCDIR}
LFLAGS = ${FLAGS} ${LIBDIR} ${LIBS}

SRC_DIR = src
BIN_DIR = bin

TARGET_OBJS = $(addprefix ${BIN_DIR}/,${OBJS})
TARGET_SRC = $(addprefix ${SRC_DIR}/,${OBJS})

BUILD_PRINT = @echo "\e[1;34mBuilding $<\e[0m"
COMPLETE_PRINT = @echo "\e[1;32mBuilding complete!\e[0m"

MKDIR_P = mkdir -p

IGNORE := $(shell bash -c "source /opt/xtal/ccp4-8.0/bin/ccp4.setup-sh; env | sed 's/=/:=/' | sed 's/^/export /' > makeenv")
include makeenv

OBJS = $(SHARED:=.o)

TARGET_OBJS = $(addprefix ${BIN_DIR}/,${OBJS})

all:
	g++ src/probe.cpp src/probe.h -o probe ${CFLAGS} ${LFLAGS}

clean:
	rm bin/*.o probe*
