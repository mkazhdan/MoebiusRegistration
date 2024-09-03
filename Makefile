LIBRARY_SRC = Include/
SRC = ./

SOFT_SOURCE = \
	SOFT1.0/cospmls.cpp \
	SOFT1.0/csecond.cpp \
	SOFT1.0/FFTcode.cpp \
	SOFT1.0/fft_grids.cpp \
	SOFT1.0/fft_grids_so3.cpp \
	SOFT1.0/FST_semi_memo.cpp \
	SOFT1.0/FST_semi_memo_fftw.cpp \
	SOFT1.0/indextables.cpp \
	SOFT1.0/makeWigner.cpp \
	SOFT1.0/naive_synthesis.cpp \
	SOFT1.0/newFCT.cpp \
	SOFT1.0/oddweights.cpp \
	SOFT1.0/OURmods.cpp \
	SOFT1.0/OURperms.cpp \
	SOFT1.0/permroots.cpp \
	SOFT1.0/primitive.cpp \
	SOFT1.0/primitive_FST.cpp \
	SOFT1.0/rotate_so3.cpp \
	SOFT1.0/rotate_so3_mem.cpp \
	SOFT1.0/seminaive.cpp \
	SOFT1.0/seminaive_fftw.cpp \
	SOFT1.0/so3_correlate_fftw.cpp \
	SOFT1.0/so3_correlate_sym.cpp \
	SOFT1.0/soft.cpp \
	SOFT1.0/soft_fftw.cpp \
	SOFT1.0/soft_fftw_pc.cpp \
	SOFT1.0/soft_fftw_wo.cpp \
	SOFT1.0/soft_sym.cpp \
	SOFT1.0/utils_so3.cpp \
	SOFT1.0/utils_vec_cx.cpp \
	SOFT1.0/weights.cpp \
	SOFT1.0/wignerTransforms.cpp \
	SOFT1.0/wignerTransforms_fftw.cpp \
	SOFT1.0/wignerTransforms_sym.cpp

PARAMETRIZATION_SOURCE=Parametrization/Parametrization.cpp
SPHERE_MAP_SOURCE=SphereMap/SphereMap.cpp
REGISTER_SOURCE=Register/Register.cpp




SOFT_TARGET=SOFT1.0
PARAMETRIZATION_TARGET=Parametrization
SPHERE_MAP_TARGET=SphereMap
REGISTER_TARGET=Register





CFLAGS += -fpermissive -fopenmp -Wno-deprecated -std=c++11
LFLAGS += -lgomp -lfftw3 -lfftw3f

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -ffast-math -g
LFLAGS_RELEASE = -O3 -g

BIN = Bin/Linux/
OBJECTS = Bin/Linux/Objects/
INCLUDE = /usr/local/include/ -IInclude

CC=gcc
CXX=g++
MD=mkdir

SOFT_OBJECTS=$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(SOFT_SOURCE))))
PARAMETRIZATION_OBJECTS=$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(PARAMETRIZATION_SOURCE))))
SPHERE_MAP_OBJECTS=$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(SPHERE_MAP_SOURCE))))
REGISTER_OBJECTS=$(addprefix $(OBJECTS), $(addsuffix .o, $(basename $(REGISTER_SOURCE))))




all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: $(BIN)
all: $(BIN)$(SOFT_TARGET)
all: $(BIN)$(PARAMETRIZATION_TARGET)
all: $(BIN)$(SPHERE_MAP_TARGET)
all: $(BIN)$(REGISTER_TARGET)




debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: $(BIN)
debug: $(BIN)$(SOFT_TARGET)
debug: $(BIN)$(PARAMETRIZATION_TARGET)
debug: $(BIN)$(SPHERE_MAP_TARGET)
debug: $(BIN)$(REGISTER_TARGET)




clean:
	rm -r $(BIN)

$(BIN):
	$(MD) -p $(BIN)
	$(MD) -p $(OBJECTS)$(SOFT_TARGET)
	$(MD) -p $(OBJECTS)$(PARAMETRIZATION_TARGET)
	$(MD) -p $(OBJECTS)$(SPHERE_MAP_TARGET)
	$(MD) -p $(OBJECTS)$(REGISTER_TARGET)
	
	
	

$(BIN)$(SOFT_TARGET): $(SOFT_OBJECTS)

$(BIN)$(PARAMETRIZATION_TARGET): $(PARAMETRIZATION_OBJECTS) $(SOFT_OBJECTS) 
	$(CXX) -o $@ $(SOFT_OBJECTS) $(PARAMETRIZATION_OBJECTS) $(LFLAGS)

$(BIN)$(SPHERE_MAP_TARGET): $(SPHERE_MAP_OBJECTS) $(SOFT_OBJECTS) 
	$(CXX) -o $@ $(SOFT_OBJECTS) $(SPHERE_MAP_OBJECTS) $(LFLAGS)

$(BIN)$(REGISTER_TARGET): $(REGISTER_OBJECTS) $(SOFT_OBJECTS) 
	$(CXX) -o $@ $(SOFT_OBJECTS) $(REGISTER_OBJECTS) $(LFLAGS)

$(BIN)%.o: $(LIBRARY_SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(BIN)%.o: $(LIBRARY_SRC)%.cpp
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(BIN)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(OBJECTS)%.o: $(LIBRARY_SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(OBJECTS)%.o: $(LIBRARY_SRC)%.cpp
	$(CC) -c -o $@ $(CFLAGS) $(INCLUDE) $<

$(OBJECTS)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<
