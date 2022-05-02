PROGRAM_NAME := camodel.exe

# Compiler
CC := g++

# Compiler flags and command
CFLAGS := -std=c++11 -O3
DFLAGS := -std=c++11 -g -ggdb -W -Wall -Wpedantic -pedantic-errors
PFLAGS := -std=c++11 -W -Wall -Wpedantic -pedantic-errors
LFLAGS := -lyaml-cpp -lgsl
OBJSDIR = ./objs
COMPILE_COMMAND := $(CC) $(CFLAGS)

# Model objects
CAModel_objs := parameters.o location.o environment.o util.o output.o
CAModel_rd_objs := reactdiff.o
CAModel_vessel_objs := vessel.o
CAModel_epithelial_objs := epithelial.o
CAModel_recepexpr_objs := recepexpr.o
CAModel_recepexpr_pcv2_discon_objs := re_pcv2_discon.o
CAModel_viralrepl_objs := viralrepl.o
CAModel_viralrepl_pcv4_discon_objs := vr_pcv4_discon.o
CAModel_agent_objs := agent.o
CAModel_agent_innate_objs := macrophage.o nk.o neutrophil.o monocyte.o

ALL_objs := $(CAModel_objs) $(CAModel_rd_objs) $(CAModel_vessel_objs) $(CAModel_epithelial_objs) \
$(CAModel_recepexpr_objs) $(CAModel_recepexpr_pcv2_discon_objs) \
$(CAModel_viralrepl_objs) $(CAModel_viralrepl_pcv4_discon_objs) \
$(CAModel_agent_objs) $(CAModel_agent_innate_objs)

# Pre-defined output directory
OUTDIR := ./Output


# compile the project
all: prelim_objsdir prelim_outdir camodel

prelim_objsdir:
	@if [ ! -d $(OBJSDIR) ]; then mkdir -p $(OBJSDIR); fi

prelim_outdir:
	@if [ ! -d $(OUTDIR) ]; then mkdir -p $(OUTDIR); fi

camodel: main.cpp $(ALL_objs)
	$(COMPILE_COMMAND) -o $(PROGRAM_NAME) $(addprefix $(OBJSDIR)/,$(ALL_objs)) main.cpp $(LFLAGS)

debug: COMPILE_COMMAND := $(CC) $(DFLAGS)
debug: all

profile: COMPILE_COMMAND := $(CC) $(PFLAGS) -pg -O3
profile: all


# CAModel Components
parameters.o: ./CAModel/parameters.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

location.o: ./CAModel/location.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

environment.o: ./CAModel/environment.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

util.o: ./CAModel/util.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

output.o: ./CAModel/output.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

# ReactDiff
reactdiff.o: ./CAModel/reactdiff/reactdiff.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

# Vessel
vessel.o: ./CAModel/vessel/vessel.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

# Epithelial
epithelial.o: ./CAModel/epithelial/epithelial.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

# Receptor Expression
recepexpr.o: ./CAModel/recepexpr/recepexpr.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

## pcv2_discon
re_pcv2_discon.o: ./CAModel/recepexpr/pcv2_discon/re_pcv2_discon.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

# Viral Replication
viralrepl.o: ./CAModel/viralrepl/viralrepl.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

## pcv4_discon
vr_pcv4_discon.o: ./CAModel/viralrepl/pcv4_discon/vr_pcv4_discon.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

# Agent
agent.o: ./CAModel/agent/agent.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

## Macrophage
macrophage.o: ./CAModel/agent/innate/macrophage.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

## NK
nk.o: ./CAModel/agent/innate/nk.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

## Neutrophil
neutrophil.o: ./CAModel/agent/innate/neutrophil.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<

## Monocyte
monocyte.o: ./CAModel/agent/innate/monocyte.cpp
	$(COMPILE_COMMAND) -o $(OBJSDIR)/$@ -c $<



# Clean targets
clean:
	rm -f $(OBJSDIR)/*.o
	rm -f $(PROGRAM_NAME)*


# Archive targets
zip:
	zip -r code.zip README.md Makefile main.cpp LICENSE CAModel/
	mv code.zip camodel_code_$$(date +%Y_%m_%d).zip
