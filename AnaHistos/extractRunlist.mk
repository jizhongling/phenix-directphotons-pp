all: extractRunlist.C
	g++ -ggdb -Wall -Werror -o extractRunlist extractRunlist.C -I$(OFFLINE_MAIN)/include -L$(OFFLINE_MAIN)/lib
