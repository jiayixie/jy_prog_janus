#
INST_DIR = $(HOME)/bin

BIN  =  temp2 #justCOR_LHEN #tmp_justcor_LHEN_NaN2Ex
#justCOR_LHEN_NAN2


fflags = -O2 -ffixed-line-length-none -m64

# if you don't use shared libraries comment the next line
#shlib = -R/opt/SUNWspro/lib

FFLAGS = $(DBG)  $(fflags) $(shlib)

CFLAGS = $(DBG) -m64

LDLIBS = -lfftw3

FC=gfortran
DBG = 
all : $(BIN)

FOBJS = dcommon.o dmultifft.o temp.o #justcor_LHEN.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS) 

install : $(BIN)
	install  -s $(BIN) $(INST_DIR)

clean :
	rm -f $(BIN) core $(FOBJS)
