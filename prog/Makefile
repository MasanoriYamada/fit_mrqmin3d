INCFLAGS = 
LDFLAGS =
LIBLIST = -lgfortran

OPTION = -O3 
DOPTION = -g -Wall -D NDEBUG 

CC = gcc
CXX = g++
GF = gfortran

############TARGET LIST#################
TARGET= fit.out
INC_DIR= ../include
INCLIST = $(wildcard $(addsuffix /*.h,$(INC_DIR)))
############SRC LIST#################
MAIN_SRC = main_fit.cpp
SRC_DIR = fit_func
SRC_FUNC=$(wildcard $(addsuffix /*.cpp,$(SRC_DIR)))
#SRC_FUNC = $(wildcard $(SRC_DIR)))


FIT_SRC_LIST=$(MAIB_SRC) $(SRC_FUNC)


SRC_1=mrqmin/DCOVSRT.f 
SRC_2=mrqmin/DGAULEG.f 
SRC_3=mrqmin/DGAUSSJ.f 
SRC_4=mrqmin/DMRQCOF.f 
SRC_5=mrqmin/DMRQMIN.f
SRC_LIST=$(SRC_1) $(SRC_2) $(SRC_3) $(SRC_4) $(SRC_5) 

#############OBJECT LIST################# 
#MAIN_OBJ = main_fit.o
#FIT_OBJ = fit.o
#FIT_OBJ_1=fit_func/func.1g1y.o
#FIT_OBJ_2=
#FIT_OBJ_3=
#FIT_OBJ_4=
#FIT_OBJ_5=
#FIT_OBJ_6=
#FIT_OBJ_LIST=$(MAIN_OBJ) $(FIT_OBJ) $(FIT_OBJ_1) $(FIT_OBJ_2) $(FIT_OBJ_3) $(FIT_OBJ_4)  $(FIT_OBJ_5) $(FIT_OBJ_6) 
#FUNC_OBJ_LIST=$(FIT_OBJ_1) $(FIT_OBJ_2) $(FIT_OBJ_3) $(FIT_OBJ_4)  $(FIT_OBJ_5) $(FIT_OBJ_6) 

OBJ_1=mrqmin/DCOVSRT.o 
OBJ_2=mrqmin/DGAULEG.o
OBJ_3=mrqmin/DGAUSSJ.o
OBJ_4=mrqmin/DMRQCOF.o 
OBJ_5=mrqmin/DMRQMIN.o
OBJ_LIST=$(OBJ_1) $(OBJ_2) $(OBJ_3) $(OBJ_4) $(OBJ_5) 

all:$(TARGET)
$(TARGET):$(INCLIST) $(MAIN_SRC) $(FIT_SRC_LIST) $(SRC_LIST) $(OBJ_LIST) $(INCLIST)
	$(CXX) $(OPTION) -o $@ $(MAIN_SRC) $(FIT_SRC_LIST) $(OBJ_LIST) $(LIBLIST)
##add new target rule 


debug:Debug.out
Debug.out:$(INCLIST) $(MAIN_SRC) $(FIT_SRC_LIST) $(SRC_LIST) $(OBJ_LIST) $(INCLIST)
	$(CXX) $(DOPTION) -o $@ $(MAIN_SRC) $(FIT_SRC_LIST) $(OBJ_LIST) $(LIBLIST)



##
$(FIT_OBJ_LIST):$(FIT_SRC_LIST) $(INCLIST)
	$(CXX) $(OPTION) -c $(FIT_SRC_LIST)
$(OBJ_LIST):$(SRC_LIST) $(INCLIST)
	make -C mrqmin


clean: 
	@rm -rf *.out *.o *~ *.bak *.deps *.dSYM
	make clean -C mrqmin
