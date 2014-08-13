EXEC = show
CC = mpicc 
OBJS = main.o saveFile.o parameterSetting.o findparam.o boundary.o loadPlasma.o loadLaser.o fieldSolve.o fieldShareY_DSX.o interpolation.o particlePush.o updateCurrent_DSX.o rearrangeParticles.o particleShareY.o removeEdge.o movingDomain.o probe.o dumpData.o clean.o filter.o boostShot.o
#OBJS = main.o saveFile.o parameterSetting.o findparam.o loadLaser.o boundary.o loadPlasma.o fieldSolve.o fieldShareY.o removeEdge.o interpolation.o clean.o particlePush.o rearrangeParticles.o particleShareY.o updateCurrent.o movingDomain.o probe.o dumpData.o filter.o
INCL = constants.h laser.h mesh.h particle.h plasma.h
LIBS = -lm
$(EXEC):$(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
