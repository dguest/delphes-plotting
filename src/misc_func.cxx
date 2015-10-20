#include "misc_func.hh"

#include <cstring>
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
  std::cerr << "usage: " << prname << ": [-h] <input> [<output>]"
	    << std::endl;
}
int CLI::read(int argc, char* argv[]) {
  using namespace std;
  if (argc == 1 || argc > 3) {
    usage(argv[0]);
    return error(1);
  } else {
    for (int argn = 1; argn < argc; argn++) {
      if (std::strcmp(argv[argn],"--help") == 0 ||
	  std::strcmp(argv[argn],"-h") == 0) {
	usage(argv[0]);
	return error(1);
      }
    }
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


// ________________________________________________________________________
// new CLI
FileCLI::FileCLI(int nargs, char* argc[]):
  m_nargs(nargs),
  m_name(argc[0]),
  m_overwrite(false)
{
  for (int nn = 0; nn < nargs; nn++) {
    m_args.push_back(argc[nn]);
  }

  int argn = 1;
  while (argn < nargs) {
    argn += check_opts(argn);
  }

  // error checking
  if (!m_overwrite && exists(m_output)){
    throw std::runtime_error("can't overwrite " + m_output);
  }
  if (m_output.size() == 0) throw std::runtime_error("no output given");
  if (m_files.size() == 0) throw std::runtime_error("no inputs given");
}

void FileCLI::usage() const
{
  std::cerr << "usage: " << m_name << ": [options] <inputs>... -o <output>"
	    << std::endl;
}

int FileCLI::check_opts(int argn) {
  std::string arg = m_args.at(argn);
  bool compound = arg.size() > 1 && (arg.at(0) == '-');
  if ( (arg == "--help") || (compound && strchr(arg.c_str(), 'h') ) ) {
    usage();
    exit(1);
  }
  if (compound) {
    if (strchr(arg.c_str(), 'f')) m_overwrite = true;
  }
  if (compound && strchr(arg.c_str(), 'o')) {
    if (! (argn + 1 < m_nargs)) {
      throw std::runtime_error("no argument for -o given");
    }
    m_output = m_args.at(argn + 1);
    return 2;
  }
  std::string file = m_args.at(argn);
  if (!exists(file)) throw std::runtime_error(file + " doesn't exist");
  m_files.push_back(file);
  return 1;
}

std::vector<std::string> FileCLI::in_files() const {
  return m_files;
}
std::string FileCLI::out_file() const {
  return m_output;
}
