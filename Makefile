# -*-Makefile-*-
CC=g++
# CXXFLAGS+=-Wall -O2
CXXFLAGS+=-g

#CXXLAPACKFLAGS is used to link lapack library, if lapack not used ,this can be ommited
CXXLAPACKFLAGS=-L${LAPACKROOT} -llapack -lrefblas -lgfortran -lm

#CXXLASPACKFLAGS=-L${LASPACKROOT} -l

sphxx=sphxx

$(sphxx):main.o
	$(CC) $(CXXFLAGS) *.o -o $(sphxx) -llaspack -lm 

main.o:main.cpp SPHSolver.o
	$(CC) $(CXXFLAGS) -c main.cpp

SPHSolver.o:SPHSolver.cpp Region.o KFile.o Model.o SPHEOS.o NblSch.o GetKnlList.o\
			SPHEqu.o DeltaT.o ControlSPH.o SPHIPCSEqu.o CalculateRange.o\
			ContinuityEqu.o Operator.o UpdatePosition.o SPHInit.o\
			IdentPtLocal.o Remesh.o GetDumProperty.o CalBndNorm.o
	$(CC) $(CXXFLAGS) -c SPHSolver.cpp

BaseFile=BasePt.o Knl.o Force.o Section.o EOS.o ControlSPH.o StatData.o	Matrix.o\
		 Node.o
$(BaseFile):%.o:%.cpp
	$(CC) $(CXXFLAGS) -c $<

Operator.o:Operator.cpp SPHPt.o Region.o Matrix.o
	$(CC) $(CXXFLAGS) -c Operator.cpp

Region.o:Region.cpp BasePt.o PtPair.o Knl.o Force.o Part.o Section.o EOS.o\
		ControlSPH.o StatData.o Matrix.o Mesh.o Node.o PtMshPair.o
	$(CC) $(CXXFLAGS) -c Region.cpp
Part.o:Part.cpp SPHPt.o
	$(CC) $(CXXFLAGS) -c Part.cpp
SectionSPH.o:SectionSPH.cpp Section.o
	$(CC) $(CXXFLAGS) -c SectionSPH.cpp
SectionNULL.o:SectionNULL.cpp Section.o
	$(CC) $(CXXFLAGS) -c SectionNULL.cpp

SPHPt.o:SPHPt.cpp BasePt.o
	$(CC) $(CXXFLAGS) -c SPHPt.cpp

Box.o:Box.cpp Region.o BasePt.o
	$(CC) $(CXXFLAGS) -c Box.cpp
PtPair.o:PtPair.cpp BasePt.o
	$(CC) $(CXXFLAGS) -c PtPair.cpp
Mesh.o:Mesh.cpp Node.o BasePt.o
	$(CC) $(CXXFLAGS) -c Mesh.cpp
PtMshPair.o:PtMshPair.cpp BasePt.o Mesh.o
	$(CC) $(CXXFLAGS) -c PtMshPair.cpp
NblSch.o:NblSch.cpp Region.o BasePt.o PtPair.o Box.o PtMshPair.o PtOperate.o
	$(CC) $(CXXFLAGS) -c NblSch.cpp
GetKnlList.o:GetKnlList.cpp Region.o BasePt.o Knl.o SPHPt.o CSPM.o Operator.o Mesh.o
	$(CC) $(CXXFLAGS) -c GetKnlList.cpp

# input & output & modelling
KFile.o:KFile.cpp Region.o SPHPt.o Part.o IdealGas.o WeaklyCompress.o SectionSPH.o Force.o\
		SectionNULL.o NblSch.o SPHIPCS.o Matrix.o
	$(CC) $(CXXFLAGS) -c KFile.cpp
Model.o:Model.cpp Region.o SPHPt.o
	$(CC) $(CXXFLAGS) -c Model.cpp

#EOS class
IdealGas.o:IdealGas.cpp EOS.o
	$(CC) $(CXXFLAGS) -c IdealGas.cpp
WeaklyCompress.o:WeaklyCompress.cpp EOS.o
	$(CC) $(CXXFLAGS) -c WeaklyCompress.cpp
SPHEOS.o:SPHEOS.cpp IdealGas.o WeaklyCompress.o Region.o SPHPt.o Part.o EOS.o
	$(CC) $(CXXFLAGS) -c SPHEOS.cpp


SPHIPCS.o:SPHIPCS.cpp Region.o SPHPt.o EOS.o
	$(CC) $(CXXFLAGS) -c SPHIPCS.cpp
CSPM.o:CSPM.cpp BasePt.o Region.o Operator.o
	$(CC) $(CXXFLAGS) -c CSPM.cpp

SPHEqu.o:SPHEqu.cpp Region.o SPHShearingForce.o SPHViscoelastic.o SPHViscosity.o SPHPressure.o\
		 ExtForce.o SPHAStress.o SPHAVForce.o SPHSmoothingEqu.o UnsteadyTerm.o EquSolve.o KFile.o\
		 SPHSTForce.o BndForce.o SPHEHD.o
	$(CC) $(CXXFLAGS) -c SPHEqu.cpp

#FCollect1 is collection of same file arch
FCollect1=SPHShearingForce.o SPHViscoelastic.o ExtForce.o SPHAVForce.o SPHSmoothingEqu.o
$(FCollect1)::%.o:%.cpp
	$(CC) $(CXXFLAGS) -c $<

SPHViscosity.o:SPHViscosity.cpp SPHPt.o Region.o SPHIPCS.o WeaklyCompress.o Matrix.o GetKnlList.o CSPM.o
	$(CC) $(CXXFLAGS) -c SPHViscosity.cpp
SPHPressure.o:SPHPressure.cpp Region.o SPHPt.o GetKnlList.o
	$(CC) $(CXXFLAGS) -c SPHPressure.cpp
SPHAStress.o:SPHAStress.cpp Region.o SPHPt.o GetKnlList.o
	$(CC) $(CXXFLAGS) -c SPHAStress.cpp
UnsteadyTerm.o:UnsteadyTerm.cpp SPHPt.o Region.o Matrix.o
	$(CC) $(CXXFLAGS) -c UnsteadyTerm.cpp
EquSolve.o:EquSolve.cpp SPHPt.o Region.o Operator.o Matrix.o
	$(CC) $(CXXFLAGS) -c EquSolve.cpp
SPHSTForce.o:SPHSTForce.cpp SPHPt.o Region.o CSPM.o Operator.o
	$(CC) $(CXXFLAGS) -c SPHSTForce.cpp
BndForce.o:BndForce.cpp Region.o SPHPt.o Operator.o
	$(CC) $(CXXFLAGS) -c BndForce.cpp

DeltaT.o:DeltaT.cpp Region.o SPHPt.o
	$(CC) $(CXXFLAGS) -c DeltaT.cpp
SPHIPCSEqu.o:SPHIPCSEqu.cpp Region.o SPHPt.o
	$(CC) $(CXXFLAGS) -c SPHIPCSEqu.cpp
CalculateRange.o:CalculateRange.cpp Region.o
	$(CC) $(CXXFLAGS) -c CalculateRange.cpp

ContinuityEqu.o:ContinuityEqu.cpp Region.o SPHPt.o BasePt.o CSPM.o GetKnlList.o
	$(CC) $(CXXFLAGS) -c ContinuityEqu.cpp
UpdatePosition.o:UpdatePosition.cpp Region.o DeltaT.o SPHPt.o BasePt.o Knl.o
	$(CC) $(CXXFLAGS) -c UpdatePosition.cpp
SPHInit.o:SPHInit.cpp BasePt.o Region.o
	$(CC) $(CXXFLAGS) -c SPHInit.cpp
IdentPtLocal.o:IdentPtLocal.cpp Region.o BasePt.o Knl.o PtPair.o CSPM.o Operator.o
	$(CC) $(CXXFLAGS) -c IdentPtLocal.cpp
Remesh.o:Remesh.cpp Region.o SPHPt.o Part.o NblSch.o Matrix.o Mesh.o GetKnlList.o
	$(CC) $(CXXFLAGS) -c Remesh.cpp
GetDumProperty.o:GetDumProperty.cpp Region.o SPHPt.o Operator.o
	$(CC) $(CXXFLAGS) -c GetDumProperty.cpp
CalBndNorm.o:CalBndNorm.cpp BasePt.o Region.o NblSch.o Operator.o KFile.o
	$(CC) $(CXXFLAGS) -c CalBndNorm.cpp
SPHEHD.o:SPHEHD.cpp Region.o BasePt.o SPHPt.o Part.o CSPM.o
	$(CC) $(CXXFLAGS) -c SPHEHD.cpp  -llaspack -lm
PtOperate.o:PtOperate.cpp SPHPt.o
	$(CC) $(CXXFLAGS) -c PtOperate.cpp


.PHONY:clean
clean:
	-rm *.o
	-rm $(sphxx)

#the following is used with emacs flymake, flymake doesnot work without that
check-syntax:
	$(CC) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)
