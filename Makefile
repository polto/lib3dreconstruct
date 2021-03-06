
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#   make - Name
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MAKE_NAME=lib3dreconstruct

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#   make - Directories
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MAKE_BINARY=bin
MAKE_MANUAL=man
MAKE_SOURCE=src
MAKE_OBJECT=obj

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#   make - Files
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MAKE_SRC=$(wildcard $(MAKE_SOURCE)/*.cpp)
MAKE_OBJ=$(addprefix $(MAKE_OBJECT)/, $(addsuffix .o, $(notdir $(basename $(MAKE_SRC)))))

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#   make - Constants
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MAKE_CC=g++
MAKE_LD=ar
MAKE_COPT=-O3 -Wall -c
MAKE_LOPT=

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#   make - All
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

all:directories $(MAKE_NAME)

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#   make - Binaries
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

$(MAKE_NAME):$(MAKE_OBJ)
	$(MAKE_LD) rcs $(MAKE_BINARY)/$(MAKE_NAME).a $^

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#   make - Objects
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

$(MAKE_OBJECT)/%.o:$(MAKE_SOURCE)/%.cpp
	$(MAKE_CC) $(MAKE_COPT) -o $@ $<

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#   make - Directories
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

directories:
	mkdir -p $(MAKE_OBJECT)
	mkdir -p $(MAKE_BINARY)

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#   make - Cleaning
# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

clean:
	rm $(MAKE_BINARY)/* -f
	rm $(MAKE_OBJECT)/*.o -f

# ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
