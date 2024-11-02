CC=			gcc
#CFLAGS=		-g -Wall -O3 -Wextra -Wno-unused-result -Wunused-parameter -fno-strict-aliasing
CFLAGS=		-g -Wall -O3 -Wno-unused-function
CPPFLAGS=
INCLUDES=	
OBJS=
PROG=		pathfinder
PROG_EXTRA=
LIBS=		-lm -lz
DESTDIR=	~/bin

.PHONY:all extra clean depend
.SUFFIXES:.c .o

ifneq ($(asan),)
		CFLAGS+=-fsanitize=address
		LIBS+=-fsanitize=address
endif

all: $(PROG)

extra: all $(PROG_EXTRA)

debug: $(PROG)
debug: CFLAGS += -DDEBUG

pathfinder: pathfinder.c path.c graph.c sstream.c misc.c kalloc.c kopen.c
		$(CC) $(CFLAGS) pathfinder.c path.c graph.c sstream.c misc.c kalloc.c kopen.c -o $@ -L. $(LIBS) $(INCLUDES)

clean:
		rm -fr *.o a.out $(PROG) $(PROG_EXTRA)

install:
		cp $(PROG) $(DESTDIR)

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(CPPFLAGS) -- *.c)

