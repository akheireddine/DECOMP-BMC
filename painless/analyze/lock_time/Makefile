libmylibc.so: mylibc.c mylibc.h
	g++ -fPIC -shared -o $@ $< -ldl -lpthread

clean:
	rm -rf libmylibc.so 
