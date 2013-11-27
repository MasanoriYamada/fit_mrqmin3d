INCFLAGS = 
LDFLAGS =
LIBLIST = -lgfortran
INCLST = 
OPTION = -O3 
DOPTION = -g -Wall

CC = gcc
CXX = g++
GF = gfortran

############SRC LIST#################

SRC_1=DCOVSRT.f 
SRC_2=DGAULEG.f 
SRC_3=DGAUSSJ.f 
SRC_4=DMRQCOF.f 
SRC_5=DMRQMIN.f
SRC_LIST=$(SRC_1) $(SRC_2) $(SRC_3) $(SRC_4) $(SRC_5) 

#############OBJECT LIST################# 

OBJ_1=DCOVSRT.o 
OBJ_2=DGAULEG.o
OBJ_3=DGAUSSJ.o
OBJ_4=DMRQCOF.o 
OBJ_5=DMRQMIN.o
OBJ_LIST=$(OBJ_1) $(OBJ_2) $(OBJ_3) $(OBJ_4) $(OBJ_5) 

all:$(OBJ_LIST)
debug:Debug



$(OBJ_LIST):$(SRC_LIST)
	$(GF) $(OPTION) -c $(SRC_LIST)

Debug:$(SRC_LIST)
	$(GF) $(DOPTION) -c $(SRC_LIST)


clean: 
	@rm -rf *.out *.o *~ *.bak *.deps

