#http://clang.llvm.org/docs/AddressSanitizer.html

CXXFLAGS= -O -g -Wextra -Wall -Wshadow -fsanitize=address -Wno-unused-variable -I/opt/local/include -I/usr/local/include/cairo


billard:billard.o graphics.o Graphics.h
	$(CXX) $(CXXFLAGS) graphics.o billard.o -o $@ -lcairo -L/opt/local/lib -lX11 -lm

hs.o: hs.cc Graphics.h
start.o: start.cc Graphics.h
billard.o: billard.cc Graphics.h
graphics.o: Graphics.h graphics.cc
