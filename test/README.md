
```bash
rm -rf build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DH5CPP_BUILD_TESTS=ON -DH5CPP_BUILD_EXAMPLES=OFF
cmake --build build -j
ctest --test-dir build --output-on-failure
```