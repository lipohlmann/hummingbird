# Adapted from MFEM (https://github.com/mfem/mfem/blob/master/doc/makefile)

SHELL = /bin/bash
MFEM_DIR ?= ..
DOXYGEN_CONF = Doxyfile

# doxygen uses: graphviz, latex
html: $(DOXYGEN_CONF)
	@# Generate the html documentation
	@( cat $(DOXYGEN_CONF) ; printf "$(MFEM_DOXYGEN_FLAGS)\n" ) | doxygen -
	@echo "<meta http-equiv=\"REFRESH\" content=\"0;URL=build/html/index.html\">" > documentation.html
	@cat warnings.log 1>&2
	@# Generate the log of undocumented methods
	@( cat $(DOXYGEN_CONF) ; echo "GENERATE_HTML=NO" ; echo "EXTRACT_ALL=NO" ; echo "WARN_LOGFILE=undoc.log" ; echo "QUIET=YES" ) | doxygen - &> /dev/null

clean:
	rm -rf $(DOXYGEN_CONF) build documentation.html *~
	rm -rf  undoc.log warnings.log

