# this makefile is not used by WorkingWiki/ProjectEngine directly -
# it is provided for use by projects.  Rather than have to construct
# a complex command line to make a recursive make call into another 
# project's directory like
# $(MAKE) -C $(PROJECT_DIR_X) -f $(RESOURCES)/makefile-before -f Makefile -f $(RESOURCES)/makefile-after target
# they can do
# $(MAKE) -C $(PROJECT_DIR_X) -f $(RESOURCES)/makefile target
# and trust this makefile to figure out which files to include.
# It has to handle at least 3 cases:
# 1. in the "old" WW, as above (where the before and after makefile are
#    called explicitly)
#    In this case its job is to include makefile-before, the project's
#    main makefile (if any) and makefile-after in that order
# 2. in the "old" WW, in an exported project.
#    In this case it should call all 3 in order.
# 3. in the "new" WW (where a GNUmakefile is provided and you should just call
#    that), in a project written for the "old" WW.
#    In this case it should just call GNUmakefile.
# 4. in the "new" WW, in an exported project.
#    In this case it should just call GNUmakefile.
# How to do this?
# I don't think there's a good way to know whether a given other project
# has an auto-GNUmakefile or a regular one.  Maybe I just write makefile-before
# and -after so that they can safely be called twice.
include $(RESOURCES)/makefile-before
ifneq ($(wildcard GNUmakefile),)
	include GNUmakefile
else ifneq ($(wildcard makefile),)
	include makefile
else ifneq ($(wildcard Makefile),)
	include Makefile
endif
include $(RESOURCES)/makefile-after
