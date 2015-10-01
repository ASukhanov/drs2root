# Author: A. Sukhanov, 09/30/2015
include Makefile.arch
#------------------------------------------------------------------------------
dqO         = drs2root.$(ObjSuf) drs2rootDict.$(ObjSuf)
dqS         = drs2root.$(SrcSuf) drs2rootDict.$(SrcSuf)
dqSO        = drs2root.$(DllSuf)

OBJS          = $(dqO)
PROGRAMS      = $(dqSO)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

drs2root:           $(dqSO)
$(dqSO):       $(dqO)
		$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(EXPLLINKLIBS) $(OutPutOpt)$@
		@echo "$@ done"

clean:
		@rm -f $(OBJS) core

distclean:      clean
		-@mv -f linearIO.root linearIO.roott
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d .def so_locations
		@rm -rf cxx_repository
		-@mv -f linearIO.roott linearIO.root
		-@cd RootShower && $(MAKE) distclean

.SUFFIXES: .$(SrcSuf)

###

drs2root.$(ObjSuf): drs2root.h
drs2rootDict.$(SrcSuf): drs2root.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

