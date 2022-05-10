PROGRAM_NAME := fft_program

MPI_PROC_COUNT ?= 8
FFT_FLTER_THRESHOLD ?= 0.2

SRC_DIR:=src
INCLUDE_DIR:=include
SCRIPTS_DIR:=scripts

DATA_DIR := data
DATA_SRC := $(DATA_DIR)/generated_data.csv
DATA_OUT_FFT := $(DATA_DIR)/fft_data.csv
DATA_OUT_IFFT := $(DATA_DIR)/rev_fft_data.csv

C_COMPILER:=mpicc
PYTHON:=python3

SRCS := \
  $(SRC_DIR)/main.c \
  $(SRC_DIR)/fft.c \
  $(SRC_DIR)/file_utils.c

C_FLAGS := \
  -lm \
  -g3 \
  -Wall

DEFINES := \
  -DFFT_FLTER_THRESHOLD=$(FFT_FLTER_THRESHOLD)

all: build data run_parallel

build:
	$(C_COMPILER) $(SRCS) $(C_FLAGS) -I$(INCLUDE_DIR) -o $(PROGRAM_NAME) $(DEFINES)

data:
	mkdir -p $(DATA_DIR)
	$(PYTHON) $(SCRIPTS_DIR)/datagen.py &

run_serial: build data
	./$(PROGRAM_NAME) $(DATA_SRC)
	$(PYTHON) $(SCRIPTS_DIR)/plotter.py $(DATA_OUT_FFT) $(DATA_OUT_IFFT)&

run_parallel: build data
	mpiexec -f nodes -n $(MPI_PROC_COUNT) ./$(PROGRAM_NAME) $(DATA_SRC)
	$(PYTHON) $(SCRIPTS_DIR)/plotter.py $(DATA_OUT_FFT) $(DATA_OUT_IFFT)&

demo:
	$(PYTHON) $(SCRIPTS_DIR)/demo.py &
	
clean:
	rm -r fft_program data