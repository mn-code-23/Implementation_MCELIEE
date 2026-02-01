# Project: shrivastava
# Makefile created by ouzdeville
CC       = gcc -Wno-implicit-function-declaration
OBJ      = gf.o rng.o matrix.o key_gen.o encrypt.o decrypt.o util.o param.o main.o
LINKOBJ  = gf.o rng.o matrix.o key_gen.o encrypt.o decrypt.o util.o param.o main.o
INCS     = 
CXXINCS  = 
BIN      = MCELIECEMS2E2025
LFLAGS=
CFLAGS= -c -Wall -I. 
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) -O3 $(LINKOBJ) -o $(BIN) $(LIBS)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<
	
gen: clean-custom
	$(CC) main_genparams.c -o main_genparams
