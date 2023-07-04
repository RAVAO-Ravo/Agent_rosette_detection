directory_name = "pyx_files"

all : 
	python3 setup.py build_ext --inplace
	rm -rf *.c
	rm -r build
	
mrproper : clean
	mkdir $(directory_name)
	mv *.pyx $(directory_name)/
	zip -r $(directory_name).zip $(directory_name)/
	rm -r $(directory_name)

clean :
	rm -rf *.c
	rm -r build