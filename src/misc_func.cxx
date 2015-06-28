#include "misc_func.hh"

#include <string>
#include <fstream>

bool exists(const std::string& file_name) {
  std::ifstream file(file_name.c_str(), std::ios::binary);
  if (!file) {
    file.close();
    return false;
  }
  file.close();
  return true;
}

std::string red(const std::string& st) {
  return "\033[31;1m" + st + "\033[m";
}
