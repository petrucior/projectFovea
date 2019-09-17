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
CFLAGS= -Wall -pedantic -ansi -g -c `pkg-config --cflags --libs opencv`

#
# A variável LFLAGS indica que opção de compilação queremos
#
LFLAGS= -Wall -g

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
LIBS=

#
# A variável OBJS indica os arquivos objetos comuns a todos os executáveis
#
OBJS= testrectangle.o #testpolygons.o 

#
# A variável ALL indica os arquivos executáveis
#
#all: match extract tracking multi
all: testrectangle #testpolygons 

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
# Regra que elimina arquivos intermediários
# 
clean:
	rm -f *.o *.cpp~

#
# Regra que elimina arquivos intermediários e executável
# 
realclean:
	rm -f *.o *.cpp~ testrectangular testpolygons

