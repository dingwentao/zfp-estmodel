CXX = mpic++
#DEFS = -DWITHOUT_COMPRESSION
CXXFLAGS = -O3 -ansi -Wall -I../inc $(DEFS)
SZPATH = /home/dtao/SZ
SZFLAG  = -I$(SZPATH)/include -L$(SZPATH)/lib

#all: zfp diffusion
all: dtzip zfp diffusion

dtzip: dtzip.cpp zfpcompress.h zfpcompress.o
	$(CXX) $(CXXFLAGS) $(SZFLAG) dtzip.cpp -o dtzip zfpcompress.o -lm -lsz -lzlib

zfp: zfp.cpp zfpcompress.h zfpcompress.o
	$(CXX) $(CXXFLAGS) zfp.cpp -o zfp zfpcompress.o

zfpcompress.o: zfpcompress.cpp
	$(CXX) $(CXXFLAGS) -c zfpcompress.cpp

diffusion: diffusion.cpp
	$(CXX) $(CXXFLAGS) diffusion.cpp -o diffusion

clean:
	rm -f dtzip zfp zfpcompress.o diffusion
