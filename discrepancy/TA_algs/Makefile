gccflags =  -g  -Wall -O3
GXX=gcc
GPP=g++

PROGRAMS = TA_new_0 TA_new_1 TA_new_2 TA_new_3 TA_split_0 TA_split_1 TA_split_2 TA_split_3 TA_shrink_0 TA_shrink_1 TA_shrink_2 TA_shrink_3 TA_shrinksplit_0 TA_shrinksplit_1 TA_shrinksplit_2 TA_shrinksplit_3
# Not quite implemented.
OTHER_PROGS = TA_parallel
P_REAL = TA_new_0 TA_new_1 TA_new_2 TA_new_3 TA_delta_0 TA_bardelta_0 TA_delta_1 TA_bardelta_1 TA_delta_2 TA_bardelta_2 TA_delta_3 TA_bardelta_3 TA_shrink_0 TA_shrink_1 TA_shrink_2 TA_shrink_3 TA_shrink_delta_0 TA_shrink_bardelta_0 TA_shrink_delta_1 TA_shrink_bardelta_1 TA_shrink_delta_2 TA_shrink_bardelta_2 TA_shrink_delta_3 TA_shrink_bardelta_3


.PHONY : TA_improved
TA_improved : TA_improved_delta TA_improved_bardelta

TA_improved_delta : TA_shrink_delta.c TA_shrink_bardelta.c TA_common.c
	${GXX} $(gccflags) -DVERSION=2 -o TA_improved_delta TA_shrink_delta.c TA_common.c -lm

TA_improved_bardelta : TA_shrink_delta.c TA_shrink_bardelta.c TA_common.c
	${GXX} $(gccflags) -DVERSION=2 -o TA_improved_bardelta TA_shrink_bardelta.c TA_common.c -lm

TA_basic : TA_base.c TA_common.c
	${GXX} $(gccflags) -DVERSION=0 -o TA_basic TA_base.c TA_common.c -lm


TA_parallel : TA_shrink_delta_parallel.cpp TA_shrink_bardelta_parallel.cpp TA_common_parallel.cpp
	${GPP} ${gccflags} -DVERSION=2 -o TA_improved_delta_par TA_shrink_delta_parallel.cpp TA_common_parallel.cpp -lm
	${GPP} ${gccflags} -DVERSION=2 -o TA_improved_bardelta_par TA_shrink_bardelta_parallel.cpp TA_common_parallel.cpp -lm

all : $(PROGRAMS)

clean :
	rm -f *.o $(PROGRAMS)

TA_new_0 : TA_base.c TA_common.c
	${GXX} $(gccflags) -DVERSION=0 -o TA_new_0 TA_base.c TA_common.c -lm

TA_new_1 : TA_base.c TA_common.c
	${GXX} $(gccflags) -DVERSION=1 -o TA_new_1 TA_base.c TA_common.c -lm

TA_new_2 : TA_base.c TA_common.c
	${GXX} $(gccflags) -DVERSION=2 -o TA_new_2 TA_base.c TA_common.c -lm

TA_new_3 : TA_base.c TA_common.c
	${GXX} $(gccflags) -DVERSION=3 -o TA_new_3 TA_base.c TA_common.c -lm


TA_split_0 : TA_base_delta.c TA_common.c TA_base_bardelta.c
	${GXX} $(gccflags) -DVERSION=0 -o TA_delta_0 TA_base_delta.c TA_common.c -lm
	${GXX} $(gccflags) -DVERSION=0 -o TA_bardelta_0 TA_base_bardelta.c TA_common.c -lm

TA_split_1 : TA_base_delta.c TA_common.c TA_base_bardelta.c
	${GXX} $(gccflags) -DVERSION=1 -o TA_delta_1 TA_base_delta.c TA_common.c -lm
	${GXX} $(gccflags) -DVERSION=1 -o TA_bardelta_1 TA_base_bardelta.c TA_common.c -lm

TA_split_2 : TA_base_delta.c TA_common.c TA_base_bardelta.c
	${GXX} $(gccflags) -DVERSION=2 -o TA_delta_2 TA_base_delta.c TA_common.c -lm
	${GXX} $(gccflags) -DVERSION=2 -o TA_bardelta_2 TA_base_bardelta.c TA_common.c -lm

TA_split_3 : TA_base_delta.c TA_common.c TA_base_bardelta.c
	${GXX} $(gccflags) -DVERSION=3 -o TA_delta_3 TA_base_delta.c TA_common.c -lm
	${GXX} $(gccflags) -DVERSION=3 -o TA_bardelta_3 TA_base_bardelta.c TA_common.c -lm


TA_shrink_0 : TA_shrink.c TA_common.c
	${GXX} $(gccflags) -DVERSION=0 -o TA_shrink_0 TA_shrink.c TA_common.c -lm

TA_shrink_1 : TA_shrink.c TA_common.c
	${GXX} $(gccflags) -DVERSION=1 -o TA_shrink_1 TA_shrink.c TA_common.c -lm

TA_shrink_2 : TA_shrink.c TA_common.c
	${GXX} $(gccflags) -DVERSION=2 -o TA_shrink_2 TA_shrink.c TA_common.c -lm

TA_shrink_3 : TA_shrink.c TA_common.c
	${GXX} $(gccflags) -DVERSION=3 -o TA_shrink_3 TA_shrink.c TA_common.c -lm

TA_shrinksplit_0 : TA_shrink_delta.c TA_shrink_bardelta.c TA_common.c
	${GXX} $(gccflags) -DVERSION=0 -o TA_shrink_delta_0 TA_shrink_delta.c TA_common.c -lm
	${GXX} $(gccflags) -DVERSION=0 -o TA_shrink_bardelta_0 TA_shrink_bardelta.c TA_common.c -lm

TA_shrinksplit_1 : TA_shrink_delta.c TA_shrink_bardelta.c TA_common.c
	${GXX} $(gccflags) -DVERSION=1 -o TA_shrink_delta_1 TA_shrink_delta.c TA_common.c -lm
	${GXX} $(gccflags) -DVERSION=1 -o TA_shrink_bardelta_1 TA_shrink_bardelta.c TA_common.c -lm

TA_shrinksplit_2 : TA_shrink_delta.c TA_shrink_bardelta.c TA_common.c
	${GXX} $(gccflags) -DVERSION=2 -o TA_shrink_delta_2 TA_shrink_delta.c TA_common.c -lm
	${GXX} $(gccflags) -DVERSION=2 -o TA_shrink_bardelta_2 TA_shrink_bardelta.c TA_common.c -lm

TA_shrinksplit_3 : TA_shrink_delta.c TA_shrink_bardelta.c TA_common.c
	${GXX} $(gccflags) -DVERSION=3 -o TA_shrink_delta_3 TA_shrink_delta.c TA_common.c -lm
	${GXX} $(gccflags) -DVERSION=3 -o TA_shrink_bardelta_3 TA_shrink_bardelta.c TA_common.c -lm

copy : $(PROGRAMS)
	scp $(P_REAL) twit:~/
