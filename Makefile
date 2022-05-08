PROGRAM_NAME:=fft_program

SRC_DIR:=src
INCLUDE_DIR:=include
SCRIPTS_DIR:=scripts
DATA_DIR:=data

C_COMPILER:=mpicc
PYTHON:=python3

SRCS = \
  $(SRC_DIR)/main.c \
  $(SRC_DIR)/complex_own.c \
  $(SRC_DIR)/fft.c 

C_FLAGS := \
  -lm \
  -g3 \
  -Wall

all:
	mkdir -p $(DATA_DIR)
	$(C_COMPILER) $(SRCS) $(C_FLAGS) -I$(INCLUDE_DIR) -o $(PROGRAM_NAME)
	mpiexec -f nodes -n 4 ./$(PROGRAM_NAME) data/generated_data.csv

clean:
	rm -r fft_program data