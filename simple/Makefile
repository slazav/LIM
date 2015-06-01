FFLAGS= -Werror -Wconversion\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -fno-range-check -O
#  -fcheck=all

FC=gfortran

all: lim

lim: lim.o ../tn/tn.a
	$(FC) $(LDFLAGS) -g -o $@ $+  $(LDLIBS)

../tn/tn.a:
	make -C ../tn


####

clean:
	rm -f *.o lim
	make -C tn clean
