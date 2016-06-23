all: spinDB.C
	g++ -Wall -Werror -o spinDB spinDB.C -I$(OFFLINE_MAIN)/include -L$(OFFLINE_MAIN)/lib -luspin `root-config --cflags --libs`
