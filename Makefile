build-rs:
	rustc raytracer.rs

build-cpp:
	c++ -o raytracer -O3 -Wall raytracer.cpp
