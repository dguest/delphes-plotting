#include "tools.hh"
#include <cstdio>
#include <cstdlib>
#include <string>
int main(int argc, char* argv[]) {
  puts("something");
  if (argc <= 1) {
    puts("no file!");
    return 1;
  }
  std::string input = argv[1];
  int testint = dumpcov(input);
  return testint;
}
