# Path to PACS Examples
PACS_ROOT = /home/ebagnati/GitHub/pacs-examples

# Paths to libraries and headers in pacs-Examples repository
EXAMPLES_INCLUDE = ${PACS_ROOT}/Examples/include
EXAMPLES_LIB = ${PACS_ROOT}/Examples/lib
# Paths to libraries and headers of muparser
MUPARSER_INCLUDE = ${PACS_ROOT}/Extras/muparser/include
MUPARSER_LIB = ${PACS_ROOT}/Extras/muparser/lib
# Paths to libraries and headers of muParserInterface
MUPARSERINTERFACE_INCLUDE = ${PACS_ROOT}/Examples/src/muParserInterface
MUPARSERINTERFACE_LIB = ${PACS_ROOT}/Examples/src/muParserInterface

# Compiler and linker variables
CXX      = mpic++
CXXFLAGS ?= -std=c++20
CPPFLAGS ?= -fopenmp -O3 -Wall -pedantic -I. -I../include -I${EXAMPLES_INCLUDE} -I${MUPARSER_INCLUDE} -I${MUPARSERINTERFACE_INCLUDE}
LDFLAGS = -L. -Wl,-rpath=$(EXAMPLES_LIB) -L${EXAMPLES_LIB} -L${MUPARSER_LIB} -L${MUPARSERINTERFACE_LIB}
LDLIBS = -lpacs -lmuparser -lmuParserInterface
export LD_LIBRARY_PATH+=:${MUPARSERINTERFACE_LIB}

LINK.o := $(LINK.cc) # Use C++ linker.
DOXYFILE ?= Doxyfile
DEPEND = make.dep

EXEC = main
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

.PHONY = all $(EXEC) $(OBJS) clean distclean $(DEPEND)

all: $(DEPEND) $(EXEC)

$(EXEC): $(OBJS)

$(OBJS): %.o: %.cpp

clean:
	$(RM) $(DEPEND)
	$(RM) *.o

distclean: clean
	$(RM) $(EXEC)
	$(RM) *.csv *.out *.bak *~

doc:
	doxygen $(DOXYFILE)

$(DEPEND): $(SRCS)
	@$(RM) $(DEPEND)
	@for file in $(SRCS); \
	do \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $${file} >> $(DEPEND); \
	done

-include $(DEPEND)