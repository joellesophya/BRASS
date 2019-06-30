
# Set library path for the compiler (ex: -L/usr/local/lib)
LIBS_PATH = -L /usr/lib
LIBS_LNX_D_FORTRAN = -lgfortran -lm

PATH_TO_EIGEN = ./src/Eigen

CPP = g++
CPPFLAGS = -w -O3

OUTPUT = ./BRASS

OBJECTS = $(wildcard ./src/*.cpp) $(wildcard ./src/*.c)
       
all: brass

brass: dsyevr $(OBJECTS)
	$(CPP) $(OBJECTS) ./src/libdsyevr.a $(LIBS_PATH) $(LIBS_LNX_D_FORTRAN) -I$(PATH_TO_EIGEN) $(CPPFLAGS) -o $(OUTPUT)

dsyevr: 
	(cd ./src/lapack_routine;$(MAKE))
	
debug: CPPFLAGS = -O0 -w -g
debug: brass

clean:
	@rm -f $(OUTPUT)
	@rm -rf ./src/*.o
	@rm -f ./src/libdsyevr.a
	@rm -f ./src/lapack_routine/*.o
	
