#pragma once
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

// POSIX file operations with string paths; no std::filesystem/link dependency.
namespace tnbp::file_io {
inline void check_path(const std::string& path) {
  if(path.empty() || path.find('\0')!=std::string::npos)
    throw std::invalid_argument("empty path or embedded NUL in path");
}
inline void fail(const std::string& operation,const std::string& path,int error) {
  throw std::runtime_error(operation+" "+path+": "+std::strerror(error));
}
inline std::string directory(const std::string& path) {
  check_path(path);auto value=path;
  while(value.size()>1 && value.back()=='/') value.pop_back();
  return value;
}
inline std::string join(const std::string& base,const std::string& relative) {
  check_path(relative);auto dir=directory(base);
  return dir+(dir=="/" ? "" : "/")+relative;
}
inline std::string filename(const std::string& path) {
  auto value=directory(path);auto slash=value.find_last_of('/');
  return slash==std::string::npos ? value : value.substr(slash+1);
}
inline std::string parent(const std::string& path) {
  auto value=directory(path);auto slash=value.find_last_of('/');
  return slash==std::string::npos ? "." : (slash==0 ? "/" : value.substr(0,slash));
}
inline bool exists(const std::string& path) {
  check_path(path);struct stat info{};
  // Count dangling links as existing too, so destinations cannot replace them.
  if(::lstat(path.c_str(),&info)==0) return true;
  const int error=errno;
  if(error==ENOENT || error==ENOTDIR) return false;
  fail("lstat",path,error);return false;
}
inline bool is_regular_file(const std::string& path) {
  check_path(path);struct stat info{};
  if(::stat(path.c_str(),&info)==0) return S_ISREG(info.st_mode);
  const int error=errno;
  if(error==ENOENT || error==ENOTDIR) return false;
  fail("stat",path,error);return false;
}
inline std::uint64_t file_size(const std::string& path) {
  check_path(path);struct stat info{};
  if(::stat(path.c_str(),&info)!=0) fail("stat",path,errno);
  if(!S_ISREG(info.st_mode) || info.st_size<0)
    throw std::runtime_error("not a regular file: "+path);
  return static_cast<std::uint64_t>(info.st_size);
}
inline bool create_directory(const std::string& path) {
  check_path(path);
  if(::mkdir(path.c_str(),0777)==0) return true;
  const int error=errno;struct stat info{};
  if(error==EEXIST && ::stat(path.c_str(),&info)==0 && S_ISDIR(info.st_mode)) return false;
  fail("mkdir",path,error);return false;
}
inline bool remove(const std::string& path) {
  check_path(path);
  if(::remove(path.c_str())==0) return true;
  const int error=errno;
  if(error==ENOENT) return false;
  fail("remove",path,error);return false;
}
inline void rename(const std::string& source,const std::string& destination) {
  check_path(source);check_path(destination);
  if(::rename(source.c_str(),destination.c_str())!=0)
    fail("rename",source+" -> "+destination,errno);
}
inline void remove_empty_directory(const std::string& path) noexcept {
  if(!path.empty() && path.find('\0')==std::string::npos) ::rmdir(path.c_str());
}
} // namespace tnbp::file_io
