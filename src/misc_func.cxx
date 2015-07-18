#include "misc_func.hh"

#include <string>
#include <fstream>
#include <iostream>

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

// === command line interface ===

CLI::CLI(): _err_code(0)
{
}

int CLI::err_code() const {
  return _err_code;
}

void CLI::usage(std::string prname) {
  std::cerr << "usage: " << prname << ": <input> [<output>]" << std::endl;
}
int CLI::read(int argc, char* argv[]) {
  using namespace std;
  if (argc == 1 || argc > 3) {
    usage(argv[0]);
    return error(1);
  }

  in_name = argv[1];
  if (!exists(in_name)) {
    cerr << in_name << " doesn't exist, exiting!" << endl;
    return error(1);
  }

  // name output
  if (argc > 2) {
    out_name = argv[2];
  } else {
    // make up name for output file if not given
    size_t dotpos = in_name.find_last_of(".");
    // if we find an extension, strip it off
    if ( ! (dotpos == 0 || dotpos == string::npos) ) {
      out_name = in_name.substr(0, dotpos);
    } else {
      out_name = in_name;
    }
    out_name.append(".h5");
  }
  if (exists(out_name)) {
    cerr << out_name << " exists, exiting" << endl;
    return error(1);
  }
  return 0;
}

int CLI::error(int error) {
  _err_code = error;
  return error;
}
