run:
	cargo run --bin raytracer && hexdump untitled.ppm && hexdump unt.ppm | head

build-rs:
	cargo build

build-cpp:
	c++ -o raytracer -O3 -Wall raytracer.cpp

# test-rs:
# 	rustc test.rs -o testrs && ./testrs && hexdump test_rs.ppm

test-rs:
	cargo run --bin testrs && hexdump test_rs.ppm && hexdump test_cpp.ppm

test-cpp:
	g++ test.cpp -o testcpp && ./testcpp && hexdump test_cpp.ppm
