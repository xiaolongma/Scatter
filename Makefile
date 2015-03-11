
SLIBS= -lm
M32= -m32
FLAG= -g -o 
scatterpath: scatterpath.o distazsub.o
	gcc $(FLAG) $@ $^ $(SLIBS)
%.o:%.c
	gcc -c $< -lm
clean: 
	/bin/rm -f *.o
