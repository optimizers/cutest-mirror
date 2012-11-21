CPPCMD		 = /usr/bin/cpp
CPPINCLUDEDIR	 = 
CPPOTHERFLAGS	 = -P -C -nostdinc -nostdinc++
CPP		 = $(CPPCMD) -I$(CPPINCLUDEDIR) $(CPPOTHERFLAGS)

UMAKE		 = ${CUTEST}/bin/umake.sh

CMP		 = gfo	 
MACHINE		 = Intel-like-PC	 
MCH		 = pc	 
COMPILE77	 = gfortran -c -ffixed-form -fno-second-underscore -Wall
COMPILE		 = gfortran -c -ffree-form -fno-second-underscore -Wall
COMPILE90	 = gfortran -c -fno-second-underscore  -Wall
LOAD		 = gfortran -ffree-form -fno-second-underscore		 
LOAD90		 = gfortran -fno-second-underscore		 
ISF9095		 = yes	
COMPILE9095	 = gfortran -c -ffree-form -fno-second-underscore -Wall
LOAD9095	 = gfortran -ffree-form -fno-second-underscore	 
FFLAGS		 = -O	 
FFLAGS90	 = -O	 
SPECIALFLAGS	= 	 
FDEBUGFLAGS	 = -g	 
FPROFILEFLAGS	= -pg	 
NBYTES		 = 8	 
CCOMPILE	 = gcc -c
CLOAD		 = gcc
CFLAGS		 = -O
CXXFLAGS = -O
SPECIALLIBS = -ldl
R1		 = 1.1920930E-07
R2		 = 5.9604646E-08
R3		 = 1.1754945E-38
R4		 = 1.1754945E-38
R5		 = 3.4028234E+38
D1		 = 2.2204460492503132D-16
D2		 = 1.1102230246251566D-16
D3		 = 2.225073858507202D-308
D4		 = 2.225073858507202D-308
D5		 = 1.797693134862314D+308

PRECISION = double
PRC = CD
PRC90 = \!D

SIZE = large
SIZ = CBIG
SIZ90 = \!BIG

PATHSEP = /
SHELL = /bin/sh

TOP = ../..
CURRENT_DIR = /double/lib

TMP = /tmp

SYSTEM = Linux
KEYSYS = CUNIX

MKDIR = mkdir
CP = cp
RM = rm -f
MV = mv
CAT = cat
CHMOD = chmod
SED = sed
LN = ln
LS = ls
AR = ar
ARFLAGS= -cr
RMDIR = rmdir
GREP = grep
AWK = gawk
HEAD = head
TAIL = tail
WC = wc
MAKE = make
LD = ld

MEXFORTRAN = mex
MEXFFLAGS = -fortran

INSTALL = install

COMMANDLIST = CMP MACHINE MCH COMPILE COMPILE90 LOAD LOAD90 COMPILE9095 LOAD9095 FFLAGS FDEBUGFLAGS FPROFILEFLAGS NBYTES CCOMPILE CLOAD CFLAGS SPECIALLIBS R1 R2 R3 R4 R5 D1 D2 D3 D4 D5 PRECISION PRC PRC90 SIZE SIZ SIZ90 SYSTEM KEYSYS MKDIR CP RM MV CAT CHMOD SED LN LS AR ARFLAGS RMDIR GREP AWK HEAD TAIL WC MAKE LD MEXFORTRAN MEXFFLAGS

.SUFFIXES: .f90 .c .a

DESTDIR =	

INTCONFIG = ./
TOOLSSRCDIR = ./
TOOLSSRCDIRSLASH = ./
TOOLSLIBDIR = ./

VPATH = $(TOOLSSRCDIR)

local.o:
	cd $(MYCUTEST)/$(PRECISION)/config; 	$(MAKE) $(MFLAGS) local.o; 	cd $(MYCUTEST)/$(PRECISION)/lib;

.f.a:
	@printf ' %-12s %-15s\t\t' "Building" $%
	$(COMPILE) $(FFLAGS) -o $(TOOLSLIBDIR)/$% $(TOOLSLIBDIR)/${%:.o=.f}
	$(AR) -crs $@ $(TOOLSLIBDIR)/$%
	@printf '[Ok]\n'
	$(RM) $(TOOLSLIBDIR)/$%

.f90.o:
	@printf ' %-12s %-15s\t\t' "Building" ${@}
	$(COMPILE90) $(FFLAGS90) -o ${@} ${@:.o=.f90}
	@printf '[Ok]\n'

#.f.o:
#	@printf ' %-12s %-15s\t\t' "Building" ${@}
#	$(COMPILE) $(FFLAGS) -o ${@} ${@:.o=.f}
#	@printf '[Ok]\n'

.c.a: $(CUTEST)/common/include/cutest.h
	@printf ' %-12s %-15s\t\t' "Building" $%
	$(CCOMPILE) $(CFLAGS) -I$(CUTEST)/common/include -o $(TOOLSLIBDIR)/$% $(TOOLSSRCDIR)/${%:.o=.c}
	$(AR) -crs $@ $(TOOLSLIBDIR)/$%
	@printf '[Ok]\n'

.c.o:
	@printf ' %-12s %-15s\t\t' "Building" ${@}
	$(CCOMPILE) $(CFLAGS) -I$(CUTEST)/common/include -o ${@} $(TOOLSSRCDIR)/${@:.o=.c}
	@printf '[Ok]\n'

all::

baseconfig:
	cd $(MYCUTEST)/$(PRECISION)/config; 	$(MAKE) $(MFLAGS) all;

all:: baseconfig

.: $(INTCONFIG)/tools.siz $(INTCONFIG)/cast.sed
	touch $(TOOLSSRCDIRSLASH)*.f

all:: .

C_O_FILES = ccutest.o ffilesop.o
CLIBFILES = libcutest.a(ccutest.o) libcutest.a(ffilesop.o)

DotOs =	cutest.o asmbe.o asmbl.o ccfg.o ccfsg.o ccifg.o ccifsg.o cdh.o cdimen.o cdimse.o cdimsh.o cdimsj.o ceh.o cfn.o cgr.o cgrdh.o cidh.o cish.o cjprod.o cnames.o cofg.o cprod.o cprod1.o creprt.o cscfg.o cscifg.o csetup.o csgr.o csgreh.o csgrsh.o csh.o csh1.o cvarty.o elgrd.o connames.o hsprd.o initw.o others.o pbname.o ubandh.o udh.o udimen.o udimse.o udimsh.o ueh.o ufn.o ugr.o ugrdh.o ugreh.o ugrsh.o unames.o uofg.o uprod.o ureprt.o usetup.o ush.o uvarty.o varnames.o

libcutest.so: local.o $(DotOs)
	$(LD) -G $(DotOs) $(MYCUTEST)/$(PRECISION)/config/local.o -o $@

libcutest.a: libcutest.a(cutest.o) libcutest.a(asmbe.o) libcutest.a(asmbl.o) libcutest.a(ccfg.o) libcutest.a(ccfsg.o) libcutest.a(ccifg.o) libcutest.a(ccifsg.o) libcutest.a(cdh.o) libcutest.a(cdimen.o) libcutest.a(cdimse.o) libcutest.a(cdimsh.o) libcutest.a(cdimsj.o) libcutest.a(ceh.o) libcutest.a(cfn.o) libcutest.a(cgr.o) libcutest.a(cgrdh.o) libcutest.a(cidh.o) libcutest.a(cish.o) libcutest.a(cjprod.o) libcutest.a(cnames.o) libcutest.a(cofg.o) libcutest.a(cprod.o) libcutest.a(cprod1.o) libcutest.a(creprt.o) libcutest.a(cscfg.o) libcutest.a(cscifg.o) libcutest.a(csetup.o) libcutest.a(csgr.o) libcutest.a(csgreh.o) libcutest.a(csgrsh.o) libcutest.a(csh.o) libcutest.a(csh1.o) libcutest.a(cvarty.o) libcutest.a(elgrd.o) libcutest.a(connames.o) libcutest.a(hsprd.o) libcutest.a(initw.o) libcutest.a(others.o) libcutest.a(pbname.o) libcutest.a(ubandh.o) libcutest.a(udh.o) libcutest.a(udimen.o) libcutest.a(udimse.o) libcutest.a(udimsh.o) libcutest.a(ueh.o) libcutest.a(ufn.o) libcutest.a(ugr.o) libcutest.a(ugrdh.o) libcutest.a(ugreh.o) libcutest.a(ugrsh.o) libcutest.a(unames.o) libcutest.a(uofg.o) libcutest.a(uprod.o) libcutest.a(ureprt.o) libcutest.a(usetup.o) libcutest.a(ush.o) libcutest.a(uvarty.o) libcutest.a(varnames.o) 
#${CLIBFILES}
	$(AR) -crs libcutest.a
	$(RM) *.o


utest: libcutest.a
	sifdecode ../../sif/ALLINITU
	$(COMPILE77) $(FFLAGS) -o ELFUN.o ELFUN.f
	$(COMPILE77) $(FFLAGS) -o GROUP.o GROUP.f
	$(COMPILE77) $(FFLAGS) -o RANGE.o RANGE.f
	$(COMPILE77) $(FFLAGS) -o EXTER.o EXTER.f
	$(LOAD90) $(FFLAGS90) -o utest ../test/utest.f90 -L./ -lcutest \
         ELFUN.o GROUP.o RANGE.o EXTER.o
	./utest
	$(RM) *.o

ctest: libcutest.a 
	sifdecode ../../sif/ALLINITC
	$(COMPILE77) $(FFLAGS) -o ELFUN.o ELFUN.f
	$(COMPILE77) $(FFLAGS) -o GROUP.o GROUP.f
	$(COMPILE77) $(FFLAGS) -o RANGE.o RANGE.f
	$(COMPILE77) $(FFLAGS) -o EXTER.o EXTER.f
	$(LOAD90) $(FFLAGS90) -o ctest ../test/ctest.f90 -L./ -lcutest \
         ELFUN.o GROUP.o RANGE.o EXTER.o
	./ctest
	$(RM) *.o

cleanlib:
	$(RM) $(DotOs) *.f *.c

all:: local.o
all:: libcutest.a

all:: cleanlib

clean:: 
	$(RM) *.CKP *.ln *.BAK *.bak core errs ,* *~ tags TAGS make.log MakeOut "#"*

purge:: 
	$(RM) libcutest.a

purge:: 
	$(RM) libcutest.so

purge:: 
	$(RM) *.f *.o

MFLAGS =

emptyrule::

clean:: 
	$(RM) *.CKP *.ln *.BAK *.bak core errs ,* *~ tags TAGS make.log MakeOut "#"*

Makefile:: 
	-@if [ -f Makefile ]; then set -x; \
	$(RM) Makefile.bak; $(MV) Makefile Makefile.bak; \
	else exit 0; fi 
	$(UMAKE_CMD) -DTOPDIR=$(TOP) -DCURDIR=$(CURRENT_DIR) -Dlinux -DIsgfo -DLargeSize -DDoublePrecision -DIsa90Compiler -DIsgcc -UNoCcomp -Ptraditional -Pw

Makefiles::

includes::

