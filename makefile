x.dat : ejercicio.x
	./ejercicio.x 

ejercicio.x : ejercicio.cpp
	c++ ejercicio.cpp -o ejercicio.x

clean :
	rm ejercicio.x 