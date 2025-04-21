#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Functions.o \
	${OBJECTDIR}/Geometry.o \
	${OBJECTDIR}/main.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/p2test \
	${TESTDIR}/TestFiles/p3base_test \
	${TESTDIR}/TestFiles/quality \
	${TESTDIR}/TestFiles/TestTransport \
	${TESTDIR}/TestFiles/qualityTrig

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64 -fopenmp -mtune=native -msse4.2 -ftree-vectorize -ffunction-sections
CXXFLAGS=-m64 -fopenmp -mtune=native -msse4.2 -ftree-vectorize -ffunction-sections

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=--64

# Link Libraries and Options
LDLIBSOPTIONS=-L/usr/lib64 -L/usr/lib64/atlas -lm -lboost_program_options -lblas -llapack /home/shared/arpackpp/external/libopenblas.a /home/shared/arpackpp/external/libsuperlu.a /home/shared/arpackpp/external/libarpack.a -ldl -lgfortran /home/shared/arpackpp/external/SuiteSparse/UMFPACK/Lib/libumfpack.a /home/shared/arpackpp/external/SuiteSparse/CHOLMOD/Lib/libcholmod.a /home/shared/arpackpp/external/SuiteSparse/COLAMD/Lib/libcolamd.a /home/shared/arpackpp/external/SuiteSparse/CCOLAMD/Lib/libccolamd.a /home/shared/arpackpp/external/SuiteSparse/CAMD/Lib/libcamd.a /home/shared/arpackpp/external/SuiteSparse/AMD/Lib/libamd.a /home/shared/arpackpp/external/SuiteSparse/metis-4.0/libmetis.a /home/shared/arpackpp/external/SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/libopenblas.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/libsuperlu.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/libarpack.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/SuiteSparse/UMFPACK/Lib/libumfpack.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/SuiteSparse/CHOLMOD/Lib/libcholmod.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/SuiteSparse/COLAMD/Lib/libcolamd.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/SuiteSparse/CCOLAMD/Lib/libccolamd.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/SuiteSparse/CAMD/Lib/libcamd.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/SuiteSparse/AMD/Lib/libamd.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/SuiteSparse/metis-4.0/libmetis.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: /home/shared/arpackpp/external/SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition ${OBJECTFILES} ${LDLIBSOPTIONS} -fopenmp

${OBJECTDIR}/Functions.o: nbproject/Makefile-${CND_CONF}.mk Functions.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Functions.o Functions.cpp

${OBJECTDIR}/Geometry.o: nbproject/Makefile-${CND_CONF}.mk Geometry.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Geometry.o Geometry.cpp

${OBJECTDIR}/Mesh.h.gch: nbproject/Makefile-${CND_CONF}.mk Mesh.h 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o "$@" Mesh.h

${OBJECTDIR}/main.o: nbproject/Makefile-${CND_CONF}.mk main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-conf ${TESTFILES}
${TESTDIR}/TestFiles/p2test: ${TESTDIR}/tests/P2_MeshTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/p2test $^ ${LDLIBSOPTIONS} `cppunit-config --libs`   

${TESTDIR}/TestFiles/p3base_test: ${TESTDIR}/tests/P3_MeshTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/p3base_test $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/quality: ${TESTDIR}/tests/QualityCheck.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/quality $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/TestTransport: ${TESTDIR}/tests/TransportCheck.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/TestTransport $^ ${LDLIBSOPTIONS} 

${TESTDIR}/TestFiles/qualityTrig: ${TESTDIR}/tests/TrigQualityCheck.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.cc}   -o ${TESTDIR}/TestFiles/qualityTrig $^ ${LDLIBSOPTIONS} 


${TESTDIR}/tests/P2_MeshTest.o: tests/P2_MeshTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -I. -std=c++11 `cppunit-config --cflags` -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/P2_MeshTest.o tests/P2_MeshTest.cpp


${TESTDIR}/tests/P3_MeshTest.o: tests/P3_MeshTest.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/P3_MeshTest.o tests/P3_MeshTest.cpp


${TESTDIR}/tests/QualityCheck.o: tests/QualityCheck.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -O2 -Wall -s -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/QualityCheck.o tests/QualityCheck.cpp


${TESTDIR}/tests/TransportCheck.o: tests/TransportCheck.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/TransportCheck.o tests/TransportCheck.cpp


${TESTDIR}/tests/TrigQualityCheck.o: tests/TrigQualityCheck.cpp 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -I. -std=c++11 -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/TrigQualityCheck.o tests/TrigQualityCheck.cpp


${OBJECTDIR}/Functions_nomain.o: ${OBJECTDIR}/Functions.o Functions.cpp 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/Functions.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Functions_nomain.o Functions.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/Functions.o ${OBJECTDIR}/Functions_nomain.o;\
	fi

${OBJECTDIR}/Geometry_nomain.o: ${OBJECTDIR}/Geometry.o Geometry.cpp 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/Geometry.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Geometry_nomain.o Geometry.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/Geometry.o ${OBJECTDIR}/Geometry_nomain.o;\
	fi

${OBJECTDIR}/Mesh_nomain.h.gch: ${OBJECTDIR}/Mesh.h.gch Mesh.h 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/Mesh.h.gch`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o "$@" Mesh.h;\
	else  \
	    ${CP} ${OBJECTDIR}/Mesh.h.gch ${OBJECTDIR}/Mesh_nomain.h.gch;\
	fi

${OBJECTDIR}/main_nomain.o: ${OBJECTDIR}/main.o main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	@NMOUTPUT=`${NM} ${OBJECTDIR}/main.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.cc) -g -O2 -Wall -DEIGEN_DONT_PARALLELIZE -DP1BASIS -I../lib -I../../../../../../shared/arpackpp/include -I../../../../../../shared/arpackpp/external/SuiteSparse/SuiteSparse_config -I../../../../../../shared/arpackpp/external/SuiteSparse/AMD/Include -I../../../../../../shared/arpackpp/examples/areig -I../../../../../../shared/arpackpp/external/SuiteSparse/UMFPACK/Include -I../../../../../../shared/arpackpp/external/SuiteSparse/CHOLMOD/Include -I../../../../../../shared/arpackpp/external/SuperLU/SRC -I../../../../../../shared/arpackpp/external/SuiteSparse/metis-4.0/Lib -I. -I/home/shared/arpackpp/include -I/home/shared/arpackpp/examples/matrices/complex -I/home/shared/arpackpp/examples/areig -I/home/shared/arpackpp/examples/areig/sym -I/home/shared/arpackpp/examples/matrices/sym -I. -I. -I. -I. -std=c++11 -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main_nomain.o main.cpp;\
	else  \
	    ${CP} ${OBJECTDIR}/main.o ${OBJECTDIR}/main_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/p2test || true; \
	    ${TESTDIR}/TestFiles/p3base_test || true; \
	    ${TESTDIR}/TestFiles/quality || true; \
	    ${TESTDIR}/TestFiles/TestTransport || true; \
	    ${TESTDIR}/TestFiles/qualityTrig || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fem_eigenvaluedecomposition

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
