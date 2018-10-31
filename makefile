CC		 = gcc
CXX		 = g++
DEBUG		 = -g
OPTFLAGS	 = -O2
LIBDIR_LAPACK	 = /usr/lib/lapack-3.8.0
LGFORTRAN_PATH	 = /usr/lib/gfortran
INCLUDE		 = $(LIBDIR_LAPACK)/LAPACKE/include
CFLAGS		 = -std=c++11 -Wall ${DEBUG} ${OPTFLAGS} -I$(INCLUDE)
CXXFLAGS         = -std=c++11 -Wall ${DEBUG} ${OPTFLAGS} -I$(INCLUDE)
LDFLAGS	 	 = -Wall ${DEBUG}
LDLIBS		 = -lbbhutil -L$(LIBDIR_LAPACK) -llapacke -llapack -lblas -L$(LGFORTRAN_PATH) -lgfortran -lm
LOCDIR		 = home/seth/research


h1-etol-slow: h1-etol-slow.o
	-${CXX} -o h1-etol-slow h1-etol-slow.o ${LDLIBS}
	rm -f h1-etol-slow.o

h1-etol-slow.o: solvers.h ekg-proc.h ekg-clean.h jacobian.h ekg-fns.h fda-fns.h fda-io.h h1-etol-slow.cpp
	$(CXX) -c $(CXXFLAGS) h1-etol-slow.cpp

h1-etol-fast: h1-etol-fast.o
	-${CXX} -o h1-etol-fast h1-etol-fast.o ${LDLIBS}
	rm -f h1-etol-fast.o

h1-etol-fast.o: solvers.h ekg-proc.h ekg-clean.h jacobian.h ekg-fns.h fda-fns.h fda-io.h h1-etol-fast.cpp
	$(CXX) -c $(CXXFLAGS) h1-etol-fast.cpp

htol-etol-slow: htol-etol-slow.o
	-${CXX} -o htol-etol-slow htol-etol-slow.o ${LDLIBS}
	rm -f htol-etol-slow.o

htol-etol-slow.o: solvers.h ekg-proc.h ekg-clean.h jacobian.h ekg-fns.h fda-fns.h fda-io.h htol-etol-slow.cpp
	$(CXX) -c $(CXXFLAGS) htol-etol-slow.cpp

htol-etol-fast: htol-etol-fast.o
	-${CXX} -o htol-etol-fast htol-etol-fast.o ${LDLIBS}
	rm -f htol-etol-fast.o

htol-etol-fast.o: solvers.h ekg-proc.h ekg-clean.h jacobian.h ekg-fns.h fda-fns.h fda-io.h htol-etol-fast.cpp
	$(CXX) -c $(CXXFLAGS) htol-etol-fast.cpp

h1-e1-slow: h1-e1-slow.o
	-${CXX} -o h1-e1-slow h1-e1-slow.o ${LDLIBS}
	rm -f h1-e1-slow.o

h1-e1-slow.o: solvers.h ekg-proc.h ekg-clean.h jacobian.h ekg-fns.h fda-fns.h fda-io.h h1-e1-slow.cpp
	$(CXX) -c $(CXXFLAGS) h1-e1-slow.cpp

h1-e1-fast: h1-e1-fast.o
	-${CXX} -o h1-e1-fast h1-e1-fast.o ${LDLIBS}
	rm -f h1-e1-fast.o

h1-e1-fast.o: solvers.h ekg-proc.h ekg-clean.h jacobian.h ekg-fns.h fda-fns.h fda-io.h h1-e1-fast.cpp
	$(CXX) -c $(CXXFLAGS) h1-e1-fast.cpp

htol-indetol-slow: htol-indetol-slow.o
	-${CXX} -o htol-indetol-slow htol-indetol-slow.o ${LDLIBS}
	rm -f htol-indetol-slow.o

htol-indetol-slow.o: solvers.h ekg-proc.h ekg-clean.h jacobian.h ekg-fns.h fda-fns.h fda-io.h htol-indetol-slow.cpp
	$(CXX) -c $(CXXFLAGS) htol-indetol-slow.cpp

htol-indetol-fast: htol-indetol-fast.o
	-${CXX} -o htol-indetol-fast htol-indetol-fast.o ${LDLIBS}
	rm -f htol-indetol-fast.o

htol-indetol-fast.o: solvers.h ekg-proc.h ekg-clean.h jacobian.h ekg-fns.h fda-fns.h fda-io.h htol-indetol-fast.cpp
	$(CXX) -c $(CXXFLAGS) htol-indetol-fast.cpp

ekg-diagind: ekg-diagind.o
	-${CXX} -o ekg-diagind ekg-diagind.o ${LDLIBS}
	rm -f ekg-diagind.o

ekg-diagind.o: solvers.h ekg-proc.h ekg-clean.h jacobian.h ekg-fns.h fda-fns.h fda-io.h ekg-diagind.cpp
	$(CXX) -c $(CXXFLAGS) ekg-diagind.cpp

ekg-conv: fda-io.h ekg-conv.cpp
	$(CXX) -c $(CXXFLAGS) ekg-conv.cpp
	$(CXX) -o ekg-conv ekg-conv.o ${LDLIBS}
	rm -f ekg-conv.o

test: test.o
	-${CXX} -o test test.o ${LDLIBS}
	rm -f test.o

test.o: ekg-proc.h ekg-clean.h ekg-fns.h fda-fns.h fda-io.h test.cpp
	$(CXX) -c $(CXXFLAGS) test.cpp

.PHONY : clean
clean :
	rm -f *.o *~

# "make" uses the following implicit rule to compile
# a .c file into a .o file:
#
#       $(CC) -c $(CFLAGS) <the-.c-file>
#
# or to compile a .cpp file into a .o file:
# 
#       $(CXX) -c $(CXXFLAGS) <the-.c-file>
#
# and the following implicit rule for linking:
#
#       $(CC) $(LDFLAGS) <all-dependent-.o-files> $(LDLIBS)
#
# (regardless of whether .o files from C++ or C source) 
# which means the above commands are inserted automatically by
# compiler if no commands given after "target : dependency"
# (note that commands signalled by tab)
#
# file.o implicitly takes file.cpp (or .c) as a dependency,
# so it is inserted by the compiler if not given explicitly
# 
# "make" uses the following rules to decide what to run:
#
#    If you type make and there's a file called makefile,
# it  runs the commands from first target of that file,
# provided dependent files are more recent than target.
#    If you type make and there's a file called Makefile,
# but none called makefile, it runs the commands from first
# target of that file, provided dependent files are more
# recent than target.
#    If you type make -f <filename>  it runs the commands
# from the first target of <filename>, provided dependent files
# are more recent than target.
#    If you type make <target> it looks for a file called
# makefile (then Makefile) and locates the target. 
#
# In all cases it will also run commands for targets that
# are dependencies of the first or given target.
