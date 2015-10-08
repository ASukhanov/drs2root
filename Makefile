# Author: A. Sukhanov, 09/30/2015
include Makefile.arch
#------------------------------------------------------------------------------
dqO         = drs4root.$(ObjSuf) drs4rootDict.$(ObjSuf) mfilter.$(ObjSuf)
dqS         = drs4root.$(SrcSuf) drs4rootDict.$(SrcSuf)
dqSO        = drs4root.$(DllSuf)

OBJS          = $(dqO)
PROGRAMS      = $(dqSO)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRAMS)

drs4root:           $(dqSO)
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

drs4root.$(ObjSuf): drs4root.h
drs4rootDict.$(SrcSuf): drs4root.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

