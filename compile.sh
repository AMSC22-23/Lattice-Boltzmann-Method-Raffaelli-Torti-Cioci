rm -r build
mkdir -p build
cmake -B build -DCMAKE_BUILD_TYPE=Debug
make -C build -j$(nproc)
mv build/lbm .