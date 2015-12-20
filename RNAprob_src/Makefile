##########
## Describe which targets the Makefile can use.
##########

##########
## Define the relative path to the RNAstructure root directory.
## Include all macro, dependency, and variable definitions.
##########

ROOTPATH=.
include ${ROOTPATH}/compiler.h
include ${ROOTPATH}/library_defines.h
include ${ROOTPATH}/dependencies.h

##########
## Define repository build actions.
##########

# Make all the executables at once.
all:
	make serial;
	make SMP;

# Make all the serial executables.
serial:
	@echo "Building of all RNAstructure serial programs started."
	@echo
	
	make RNAprob;
	make scorer;
	@echo
	@echo "Building of the serial RNAstructure programs finished."


# Copy the executables to the /usr/local directory.
install:
	cp -r exe/ /usr/local/RNAstructure

##########
## Define targets.
##########



# Build the Fold text interface.
RNAprob: exe/RNAprob
exe/RNAprob: fold/Fold.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} fold/Fold.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the scorer interface.
scorer: exe/scorer
exe/scorer: scorer/Scorer_Interface.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES}
	${LINK} scorer/Scorer_Interface.o ${CMD_LINE_PARSER} ${STRUCTURE_SCORER} ${RNA_FILES}



##########
## Cleanup.
## Object cleanup removes all temporary build objects.
## Executable cleanup removes all possible executables.
##########

# Remove object files and any temporary files from building.
clean:
	find . -depth -name '*~' -delete
	find . -depth -name '*.o' -delete
	find . -depth -name '*.class' -delete

# Remove object files and executables.
realclean: clean
	rm -f RNA_class/*_class
	find exe -maxdepth 1 -type f ! -name RNAstructureScript ! -name RNAstructure.bat -delete
	find exe -type d ! -name exe ! -name CVS | xargs rm -rf

