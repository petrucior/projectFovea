#
#	Makefile projectFovea
#

#
# A variável CC indica que compilador que será utilizado
#
CC= g++

#
# A variável CFLAGS indica que opção de compilação queremos
#
CFLAGS= -Wall -pedantic -ansi -g -c `pkg-config --libs opencv`

#
# A variável LFLAGS indica que opção de compilação queremos
# MACOS: -openmp  Unix: -fopenmp 
LFLAGS= -Wall -g #-openmp #-fopenmp

#
# A variável DIR indica o diretório dos arquivos de cabeçalho do projeto
#
DIR= ~/projectFovea

#
# A variável INCS indica o caminho dos arquivos de cabeçalho
#
INCS= -I../ -I/usr/local/include/opencv2 -I/usr/local/include

#
# A variável LIBS indica o caminho das bibliotecas e as bibliotecas usadas na ligação
#
LIBS= `pkg-config --libs opencv`

#
# A variável OBJS indica os arquivos objetos comuns a todos os executáveis
#
OBJS= testmulti.o #testrectangle.o testpolygons.o testblocks.o testlevel.o testfovea.o testmulti.o testdetect.o teststatistics.o testfunction.o testsurf.o

#
# A variável ALL indica os arquivos executáveis
#
#all: match extract tracking multi
all: testmulti #testrectangle testpolygons testblocks testlevel testfovea testmulti testdetect teststatistics testfunction testsurf

#
# Regra de compilação e ligação do executáveis
# 
testrectangular: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o testrectangle $(LIBS) -O3

#
# Regra de compilação do arquivo objeto testrectangular.o
#
testrectangular.o: testrectangle.cpp $(DIR)/shape.hpp $(DIR)/rectangle.hpp
	$(CC) $(CFLAGS) testrectangle.cpp -o testrectangle.o

#
# Regra de compilação e ligação do executáveis
# 
testpolygons: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o testpolygons $(LIBS) -O3

#
# Regra de compilação do arquivo objeto testrectangular.o
#
testpolygons.o: testpolygons.cpp $(DIR)/shape.hpp $(DIR)/polygons.hpp
	$(CC) $(CFLAGS) testpolygons.cpp -o testpolygons.o

#
# Regra de compilação e ligação do executáveis
# 
testblocks: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o testblocks $(LIBS) -O3

#
# Regra de compilação do arquivo objeto testblocks.o
#
testblocks.o: testblocks.cpp $(DIR)/shape.hpp $(DIR)/blocks.hpp
	$(CC) $(CFLAGS) testblocks.cpp -o testblocks.o


#
# Regra de compilação e ligação do executáveis
# 
testlevel: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o testlevel $(LIBS) -O3

#
# Regra de compilação do arquivo objeto testlevel.o
#
testlevel.o: testlevel.cpp $(DIR)/shape.hpp $(DIR)/rectangle.hpp $(DIR)/polygons.hpp $(DIR)/level.hpp
	$(CC) $(CFLAGS) testlevel.cpp -o testlevel.o

#
# Regra de compilação e ligação do executáveis
# 
testfovea: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o testfovea $(LIBS) -O3

#
# Regra de compilação do arquivo objeto testfovea.o
#
testfovea.o: testfovea.cpp $(DIR)/shape.hpp $(DIR)/rectangle.hpp $(DIR)/polygons.hpp $(DIR)/level.hpp $(DIR)/feature.hpp $(DIR)/fovea.hpp
	$(CC) $(CFLAGS) testfovea.cpp -o testfovea.o

#
# Regra de compilação e ligação do executáveis
# 
testmulti: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o testmulti $(LIBS) -O3

#
# Regra de compilação do arquivo objeto testfovea.o
#
testmulti.o: testmulti.cpp $(DIR)/shape.hpp $(DIR)/rectangle.hpp $(DIR)/polygons.hpp $(DIR)/level.hpp $(DIR)/feature.hpp $(DIR)/fovea.hpp $(DIR)/multifovea.hpp
	$(CC) $(CFLAGS) testmulti.cpp -o testmulti.o

#
# Regra de compilação e ligação do executáveis
# 
testdetect: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o testdetect $(LIBS) -O3

#
# Regra de compilação do arquivo objeto testdetect.o
#
testdetect.o: testdetect.cpp $(DIR)/shape.hpp $(DIR)/rectangle.hpp $(DIR)/polygons.hpp $(DIR)/level.hpp $(DIR)/feature.hpp $(DIR)/fovea.hpp $(DIR)/multifovea.hpp
	$(CC) $(CFLAGS) testdetect.cpp -o testdetect.o

#
# Regra de compilação e ligação do executáveis
# 
teststatistics: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o teststatistics $(LIBS) -O3

#
# Regra de compilação do arquivo objeto teststatistics.o
#
teststatistics.o: teststatistics.cpp $(DIR)/shape.hpp $(DIR)/rectangle.hpp $(DIR)/polygons.hpp $(DIR)/level.hpp $(DIR)/feature.hpp $(DIR)/fovea.hpp $(DIR)/multifovea.hpp $(DIR)/statistics.hpp
	$(CC) $(CFLAGS) teststatistics.cpp -o teststatistics.o


#
# Regra de compilação e ligação do executáveis
# 
testfunction: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o testfunction $(LIBS) -O3

#
# Regra de compilação do arquivo objeto testfunction.o
#
testfunction.o: testfunction.cpp $(DIR)/shape.hpp $(DIR)/rectangle.hpp $(DIR)/polygons.hpp $(DIR)/level.hpp $(DIR)/feature.hpp $(DIR)/fovea.hpp $(DIR)/multifovea.hpp $(DIR)/statistics.hpp
	$(CC) $(CFLAGS) testfunction.cpp -o testfunction.o

#
# Regra de compilação e ligação do executáveis
# 
testsurf: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o testsurf $(LIBS) -O3

#
# Regra de compilação do arquivo objeto testfunction.o
#
testsurf.o: testsurf.cpp $(DIR)/shape.hpp $(DIR)/rectangle.hpp $(DIR)/polygons.hpp $(DIR)/level.hpp $(DIR)/feature.hpp $(DIR)/foveatedHessianDetector.hpp $(DIR)/fovea.hpp 
	$(CC) $(CFLAGS) testsurf.cpp -o testsurf.o


#
# Regra que elimina arquivos intermediários
# 
clean:
	rm -f *.o *.cpp~

#
# Regra que elimina arquivos intermediários e executável
# 
realclean:
	rm -f *.o *.cpp~ testrectangle testpolygons testblocks testlevel testfovea testmulti testdetect teststatistics testfunction testsurf

