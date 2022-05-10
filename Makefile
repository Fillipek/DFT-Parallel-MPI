PROGRAM_NAME:=fft_program

SRC_DIR:=src
INCLUDE_DIR:=include
SCRIPTS_DIR:=scripts
DATA_DIR:=data

C_COMPILER:=mpicc
PYTHON:=python3

SRCS = \
  $(SRC_DIR)/main.c \
  $(SRC_DIR)/fft.c \
  $(SRC_DIR)/file_utils.c \

C_FLAGS := \
  -lm \
  -g3 \
  -Wall

all:
	mkdir -p $(DATA_DIR)
	$(PYTHON) $(SCRIPTS_DIR)/datagen.py
	$(C_COMPILER) $(SRCS) $(C_FLAGS) -I$(INCLUDE_DIR) -o $(PROGRAM_NAME)
	./$(PROGRAM_NAME) data/generated_data.csv
	$(PYTHON) $(SCRIPTS_DIR)/plotter.py &

demo:
	$(PYTHON) $(SCRIPTS_DIR)/demo.py &
	
clean:
	rm -r fft_program data