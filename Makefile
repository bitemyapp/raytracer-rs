build-rs:
	cargo build

build-cpp:
	c++ -o raytracer -O3 -Wall raytracer.cpp
