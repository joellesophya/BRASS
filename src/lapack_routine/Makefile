FORTRAN = gfortran 
FOPS = -O2 -frecursive

# for dsyevr
DLASRC = \
disnan.o  dlaed3.o  dlaed9.o    dlanst.o  dlarft.o  dlassq.o  dormtr.o  dsytrd.o  lsame.o \
dlacpy.o  dlaed4.o  dlaeda.o    dlansy.o  dlartg.o  dlatrd.o  dstedc.o  ieeeck.o  xerbla.o \
dlae2.o   dlaed5.o  dlaev2.o    dlapy2.o  dlascl.o  dorm2l.o  dsteqr.o  iladlc.o \
dlaed0.o  dlaed6.o  dlaisnan.o  dlarfb.o  dlaset.o  dorm2r.o  dsterf.o  iladlr.o \
dlaed1.o  dlaed7.o  dlamch.o    dlarf.o   dlasr.o   dormql.o  dsyevr.o  ilaenv.o \
dlaed2.o  dlaed8.o  dlamrg.o    dlarfg.o  dlasrt.o  dormqr.o  dsytd2.o  iparmq.o \
daxpy.o   dcopy.o   ddot.o      dgemm.o   dgemv.o   dger.o    dnrm2.o   drot.o   \
dscal.o   dswap.o   dsymv.o     dsyr2.o   dsyr2k.o  dtrmm.o   dtrmv.o   idamax.o \
dstebz.o  dstein.o  dstemr.o    dlaebz.o  dlarnv.o  dlagtf.o  dlagts.o  dlarrc.o \
dlarre.o  dlarrj.o  dlarrr.o    dlarrv.o  dlaruv.o  dlarra.o  dlarrk.o  dlasq2.o \
dlarrd.o  dlarrb.o  dlarrf.o    dlar1v.o  dlasq3.o  dlasq4.o  dlasq5.o  dlasq6.o \
dlaneg.o

all:	$(DLASRC)
	ar cr ../libdsyevr.a $(DLASRC) 
	ranlib ../libdsyevr.a
	

$(DLASRC): $(FRC)

FRC:
	@FRC=$(FRC)

.f.o:
	$(FORTRAN) -static $(FOPS) -c $< -o $@ 
