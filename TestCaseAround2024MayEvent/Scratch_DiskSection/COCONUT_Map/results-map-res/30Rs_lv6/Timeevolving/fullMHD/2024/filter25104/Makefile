#makefile for the Project of resultdrawing
#---------by Hp

  f90 = mpif90 
  MKL_FLAG = -lopenblas -lpthread  #-lmkl_lapack95 -lmkl_em64t -lguide -lpthread
  DEBUG_FLAGS = -g -O3 -ffree-line-length-none  #-O2 -r8 -i4 -fpe0 -traceback

objects = result_drawing_COCONUT_Seq.o \
			   	  

./sperical.YY :  $(objects)
#./impmeridianto20Rsaout :  $(objects)
#./rhoVraout :  $(objects)
#./impBxyzout :  $(objects)
	${f90} $(DEBUG_FLAGS)  -o $@ $(objects) $(MKL_FLAG)

#. f90 文件生成 .o 的默认模式
%.o : %.f90
	${f90} $(DEBUG_FLAGS) -c $<  -o $@

%.o : %.for
	${f90} -c $<  -o $@ 

%.o : %.f
	${f90} -c $<  -o $@ 


result_drawing_COCONUT_Seq.o : result_drawing_COCONUT_Seq.f90

clean:
	-rm *.o *.mod

