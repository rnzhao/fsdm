override CFLAGS += -O3 -Wall -Wno-unknown-pragmas -g -I./src
LDLIBS = -lm -lz

ifeq ($(OS),Windows_NT)
	OS_TYPE = Windows
else
	OS_TYPE = $(shell uname -s)
endif

ifeq ($(OS_TYPE),Darwin)
	CC = clang
	CFLAGS += -std=gnu99
else
	CC = gcc
endif

ifneq ($(DEBUG),1)
	DEFINES = -DNDEBUG
endif

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)

TARGET = fsdm

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(DEFINES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f *.o $(TARGET)
