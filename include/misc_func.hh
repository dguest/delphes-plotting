#ifndef MISC_FUNC_HH
#define MISC_FUNC_HH

#include <string>

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


#endif
