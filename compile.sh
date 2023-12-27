cmake -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
make -C build -j$(nproc)