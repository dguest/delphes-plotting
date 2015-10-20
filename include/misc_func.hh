#ifndef MISC_FUNC_HH
#define MISC_FUNC_HH

#include <string>
#include <vector>

bool exists(const std::string& file_name);
std::string red(const std::string& string);

struct CLI
{
  CLI();
  std::string out_name;
  std::string in_name;
  int err_code() const;
  void usage(std::string prname);
  int read(int argc, char* argv[]);
private:
  int _err_code;
  int error(int);
};

class FileCLI
{
public:
  FileCLI(int argc, char* argv[]);
  std::vector<std::string> in_files() const;
  std::string out_file() const;
private:
  int check_opts(int argn);
  void usage() const;
  std::vector<std::string> m_files;
  std::string m_output;
  const int m_nargs;
  std::string m_name;
  std::vector<std::string> m_args;
  bool m_overwrite;
};

#endif
