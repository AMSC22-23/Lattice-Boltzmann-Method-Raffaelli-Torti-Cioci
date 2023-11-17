#include <iostream>

int main() {
  // write a simple program that segfaults
  int i = 2;
  int *p = nullptr;
  i = 3;
  *p = 0;
  std::cout << "Hello, World!" << std::endl;
  return 0;
}