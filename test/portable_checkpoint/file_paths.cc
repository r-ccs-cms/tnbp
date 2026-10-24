#include "tnbp/framework/fileutility.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <type_traits>
namespace f=tnbp::file_io;
void check(bool good,const char* why) {if(!good) throw std::runtime_error(why);}
int main() {
 try {
  char pattern[]="/tmp/tnbp-paths-XXXXXX";char* temp=::mkdtemp(pattern);
  check(temp!=nullptr,"mkdtemp");const std::string root=temp;
  const auto dir=f::join(root,"directory with spaces");
  check(f::create_directory(dir),"mkdir");check(!f::create_directory(dir),"existing directory");
  check(f::directory(dir+"///")==dir,"trailing separators");
  check(f::join(dir+"/","file")==dir+"/file","join");
  const auto file=f::join(dir,"file");{std::ofstream out(file);out<<"abc";}
  check(f::is_regular_file(file)&&f::file_size(file)==3,"file size");
  check(!f::is_regular_file(dir),"directory not regular");
  const auto renamed=f::join(dir,"renamed");f::rename(file,renamed);
  check(!f::exists(file)&&f::exists(renamed),"rename");
  const auto link=f::join(dir,"dangling");check(::symlink("absent",link.c_str())==0,"symlink");
  check(f::exists(link)&&!f::is_regular_file(link),"dangling destination exists");
  bool rejected=false;try {f::exists(std::string("a\0b",3));}catch(const std::invalid_argument&){rejected=true;}
  check(rejected,"reject NUL");rejected=false;
  try {f::directory("");}catch(const std::invalid_argument&){rejected=true;}
  check(rejected,"reject empty path");
  check(f::remove(link)&&f::remove(renamed),"unlink");check(!f::remove(file),"remove missing");
  check(f::remove(dir)&&f::remove(root),"remove empty directories");
  std::cout<<"PASS POSIX string paths\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
