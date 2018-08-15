# Path to SHAN library
SHAN_DIR=$(HOME)/bsc_workshop/shan
# Path to GPI2 library -No need for Crays
#GPI2_DIR = $(HOME)/GPI2-1.3.0-2018

INCLUDE_DIR += $(SHAN_DIR)/include
LIBRARY_DIR += $(SHAN_DIR)/lib64

TARGET   = shanLBM.exe
CC       = cc
# compiling flags here
CFLAGS   = -h std=c11 -DDEBUG -DUSE_SHAN -DUSE_GASPI
CFLAGS += $(addprefix -I,$(INCLUDE_DIR))

#LINKER   = cc -o
# linking flags here
#LDFLAGS += $(addprefix -L,$(LIBRARY_DIR))
FLAGS= $(LIBRARY_DIR)/libSHAN.a

# change these to set the proper directories where each files shoould be
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
rm       = rm -f

# Do not modify this part of the script

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CC) -o $@ $(LDFLAGS) $(OBJECTS)


$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@


.PHONEY: clean
clean:
	$(rm) $(OBJECTS) 


.PHONEY: remove
remove: clean
	$(rm) $(BINDIR)/$(TARGET)



