CC = gcc
CFLAGS = -I. -I../ -O2
LINKER = gcc
LFLAGS = -Wall -I../ -O2 -lm -lgsl -lgslcblas
TARGET = EvsR_c

BINDIR = ../../bin
#GSL_SRCDIR = ../edited_gsl_src
ENERGY_SRCDIR = ../energy_src
SHARED_SRCDIR = ..
OBJDIR = ../../obj

LOCAL_SRC = EvsR_c.c

ENERGY_SRC := pinvs.c red.c bksub.c energy.c nrutil.c polint.c trapzd.c \
              solvde.c difeq.c finite_differences.c qromb.c shared.c spline.c

SHARED_SRC := utilities.c

#EDITED_GSL_SRC := $(wildcard $(GSL_SRCDIR)/*.c)

INCLUDES := headerfile.h
#EDITEDGSL_INC := $(wildcard $(GSL_SRCDIR)/*.h)


#EDITED_GSL_OBJS := $(EDITED_GSL_SRC:$(GSL_SRCDIR)/%.c=$(OBJDIR)/%.o)
ENERGY_OBJS := $(ENERGY_SRC:%.c=$(OBJDIR)/%.o)
SHARED_OBJS := $(SHARED_SRC:%.c=$(OBJDIR)/%.o)
LOCAL_OBJS := $(LOCAL_SRC:%.c=$(OBJDIR)/%.o)
#OBJECTS := $(ENERGY_OBJS) $(EDITED_GSL_OBJS) $(LOCAL_OBJS)
OBJECTS := $(ENERGY_OBJS) $(SHARED_OBJS) $(LOCAL_OBJS)
rm = rm -f

$(BINDIR)/$(TARGET): $(OBJECTS)
	@$(LINKER) $(OBJECTS) $(LFLAGS) -o $@
	@echo "Linking complete!"

$(LOCAL_OBJS) : $(OBJDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(SHARED_OBJS) : $(OBJDIR)/%.o: $(SHARED_SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

$(ENERGY_OBJS) : $(OBJDIR)/%.o: $(ENERGY_SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@
	@echo "Compiled "$<" successfully!"

#$(EDITED_GSL_OBJS) : $(OBJDIR)/%.o: $(GSL_SRCDIR)/%.c
#	$(CC) $(CFLAGS) -c $< -o $@
#	@echo "Compiled "$<" successfully!"

PHONY: clean
clean:
	$(rm) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	$(rm) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"
