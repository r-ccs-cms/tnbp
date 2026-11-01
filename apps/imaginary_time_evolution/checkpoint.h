#pragma once
#include "config.h"
#include "tnbp/ptns/function/portable_checkpoint.h"
#include <optional>

namespace ite {
namespace fs=tnbp::file_io;
struct Resume {
  std::string state_digest,plan_digest;
  Config config;
  std::size_t next_step=0;
  double operator_log_scale=0;
};
inline std::string encode_resume(const Resume& r) {
  std::ostringstream out;out.imbue(std::locale::classic());out<<std::setprecision(17);
  out<<"TNBP_IMAGINARY_TIME_RESUME 2\nSTATE "<<std::quoted(r.state_digest)
     <<"\nCOMPRESSION leaf-batches-v1\nPLAN_VERSION 2\nPLAN "<<std::quoted(r.plan_digest)
     <<"\nCONFIG "<<std::quoted(encode_config(r.config))<<"\nNEXT_STEP "<<r.next_step
     <<"\nOPERATOR_LOG_SCALE "<<r.operator_log_scale<<"\nPHASE step-complete\nEND\n";
  auto body=out.str();return body+"SHA256 "+tnbp::portable_detail::sha256(body)+"\n";
}
inline std::optional<Resume> read_resume(const std::string& directory,const std::string& digest) {
  require(!fs::exists(fs::join(directory,"imaginary_time.resume.pending")),"incomplete resume file");
  auto path=fs::join(directory,"imaginary_time.resume");if(!fs::exists(path))return std::nullopt;
  require(fs::file_size(path)<=64*1024*1024,"resume file exceeds 64 MiB");
  std::ifstream file(path,std::ios::binary);require(bool(file),"cannot open resume file");
  std::ostringstream bytes;bytes<<file.rdbuf();require(!file.bad(),"resume read failed");auto text=bytes.str();
  auto checksum=text.rfind("SHA256 ");
  require(checksum!=std::string::npos && text.substr(checksum)=="SHA256 "+tnbp::portable_detail::sha256(text.substr(0,checksum))+"\n","resume checksum mismatch");
  std::istringstream in(text.substr(0,checksum));in.imbue(std::locale::classic());
  auto token=[&]{std::string v;require(bool(in>>v),"truncated resume");return v;};
  auto tag=[&](const std::string& v){require(token()==v,"invalid resume tag/version: "+v);};
  Resume r;tag("TNBP_IMAGINARY_TIME_RESUME");tag("2");tag("STATE");require(bool(in>>std::quoted(r.state_digest)),"missing state digest");
  tag("COMPRESSION");tag("leaf-batches-v1");
  tag("PLAN_VERSION");tag("2");tag("PLAN");require(bool(in>>std::quoted(r.plan_digest)),"missing plan digest");
  std::string config;tag("CONFIG");require(bool(in>>std::quoted(config)),"missing config");r.config=decode_config(config);
  tag("NEXT_STEP");r.next_step=count(token());tag("OPERATOR_LOG_SCALE");r.operator_log_scale=real(token());
  tag("PHASE");tag("step-complete");tag("END");in>>std::ws;require(in.eof(),"trailing resume data");
  require(r.state_digest==digest,"resume belongs to another state");require(r.next_step<=total_steps(r.config.schedule),"resume step exceeds schedule");return r;
}
template<class SaveState>
tnbp::PortableCheckpoint save(const std::string& destination,Resume r,MPI_Comm comm,SaveState save_state) {
  auto target=fs::directory(destination),staging=target+".imaginary-time.pending",state=fs::join(staging,"state");
  int rank;MPI_Comm_rank(comm,&rank);tnbp::portable_detail::agree(target,comm);
  tnbp::portable_detail::phase(comm,[&] {if(rank==0) {
    require(!fs::exists(target),"checkpoint destination already exists");
    require(fs::create_directory(staging),"checkpoint staging already exists");
  }});
  auto result=save_state(state);r.state_digest=result.state_digest;
  tnbp::portable_detail::phase(comm,[&] {if(rank==0) {
    auto pending=fs::join(state,"imaginary_time.resume.pending");
    std::ofstream f;f.exceptions(std::ios::failbit|std::ios::badbit);auto encoded=encode_resume(r);require(encoded.size()<=64*1024*1024,"resume file exceeds 64 MiB");
    f.open(pending,std::ios::binary);f<<encoded;f.close();
    fs::rename(pending,fs::join(state,"imaginary_time.resume"));
    require(!fs::exists(target),"checkpoint destination already exists");fs::rename(state,target);fs::remove_empty_directory(staging);
  }});
  result.directory=target;return result;
}
} // namespace ite
