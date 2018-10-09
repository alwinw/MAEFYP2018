#############################################################################
# Makefile for panel-method codes.
# RCS Information:
# $Id: Makefile,v 1.1.1.1 2008/04/01 01:16:57 hmb Exp $
#############################################################################

OPTIONS  = -m64 -O2 -w

CFLAGS   = $(OPTIONS)
CPPFLAGS = -I$(HOME)/include
LDFLAGS  = -lm


OBJ       = panel2d.o lu.o singular.o
HEADER    = panel.h

panel2d: $(OBJ) $(HEADER)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

OBJ: $(HEADER)

clean:
	-rm -f *.o a.out *~ core




