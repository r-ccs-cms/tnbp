#pragma once
#include "tnbp/ptns/function/portable_checkpoint.h"
#include <charconv>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

// Estimator restart policy belongs to this application, not the TPS format.
namespace estimator_checkpoint {
namespace fs=tnbp::file_io;
struct Resume {
  std::string state_digest, circuit_sha256, settings;
  std::size_t next_step = 0;
};
inline void require(bool condition, const char* message) {
  if(!condition) throw std::runtime_error(std::string("estimator checkpoint: ")+message);
}
inline std::string encode(const Resume& resume) {
  std::ostringstream out;
  out << "TNBP_ESTIMATOR_RESUME 1\nSTATE " << std::quoted(resume.state_digest)
      << "\nCIRCUIT " << std::quoted(resume.circuit_sha256)
      << "\nNEXT_STEP " << resume.next_step
      << "\nSETTINGS " << std::quoted(resume.settings)
      << "\nPHASE cycle-complete\nEND\n";
  auto body=out.str();
  return body+"SHA256 "+tnbp::portable_detail::sha256(body)+"\n";
}
inline std::optional<Resume> read(const std::string& directory,
                                 const std::string& state_digest) {
  require(!fs::exists(fs::join(directory,"estimator.resume.pending")),
          "incomplete estimator resume file");
  const auto path=fs::join(directory,"estimator.resume");
  if(!fs::exists(path)) return std::nullopt;
  require(fs::file_size(path)<=16384, "resume file too large");
  std::ifstream file(path, std::ios::binary);
  require(bool(file), "cannot open resume file");
  std::ostringstream bytes; bytes<<file.rdbuf();
  require(!file.bad(), "cannot read resume file");
  const auto text=bytes.str();
  const auto checksum=text.rfind("SHA256 ");
  require(checksum!=std::string::npos &&
          text.substr(checksum)=="SHA256 "+tnbp::portable_detail::sha256(text.substr(0,checksum))+"\n",
          "resume checksum mismatch");
  std::istringstream in(text.substr(0,checksum));
  auto tag=[&](const char* expected) { std::string value;
    require(bool(in>>value) && value==expected, "unsupported or malformed resume file"); };
  tag("TNBP_ESTIMATOR_RESUME"); tag("1");
  Resume resume;
  tag("STATE"); require(bool(in>>std::quoted(resume.state_digest)), "missing state digest");
  tag("CIRCUIT"); require(bool(in>>std::quoted(resume.circuit_sha256)), "missing circuit digest");
  tag("NEXT_STEP"); std::string step; require(bool(in>>step), "missing next step");
  auto result=std::from_chars(step.data(),step.data()+step.size(),resume.next_step);
  require(result.ec==std::errc{} && result.ptr==step.data()+step.size(), "invalid next step");
  tag("SETTINGS"); require(bool(in>>std::quoted(resume.settings)), "missing settings");
  tag("PHASE"); tag("cycle-complete"); tag("END");
  in>>std::ws; require(in.eof(), "trailing resume data");
  require(resume.state_digest==state_digest, "resume state digest mismatch");
  return resume;
}
inline void write(const std::string& directory, const Resume& resume) {
  const auto pending=fs::join(directory,"estimator.resume.pending");
  std::ofstream file;
  file.exceptions(std::ios::failbit|std::ios::badbit);
  file.open(pending,std::ios::binary); file<<encode(resume); file.close();
  fs::rename(pending,fs::join(directory,"estimator.resume"));
}
// Stage the generic state and this app's resume data together. The state
// library knows nothing about app metadata or this outer publication boundary.
template<class SaveState>
tnbp::PortableCheckpoint save(const std::string& destination,
                              Resume resume, MPI_Comm comm, SaveState save_state) {
  const auto target=fs::directory(destination);
  const std::string staging=target+".estimator.pending";
  const auto state=fs::join(staging,"state");
  int rank=0; MPI_Comm_rank(comm,&rank);
  tnbp::portable_detail::agree(target,comm);
  tnbp::portable_detail::phase(comm,[&] {
    if(rank==0) {
      require(!fs::exists(target), "destination already exists");
      require(fs::create_directory(staging), "estimator staging destination exists");
    }
  });
  // A failure leaves staging for inspection, never a published half-checkpoint.
  auto checkpoint=save_state(state);
  resume.state_digest=checkpoint.state_digest;
  tnbp::portable_detail::phase(comm,[&] {
    if(rank==0) {
      write(state,resume);
      require(!fs::exists(target), "destination already exists");
      fs::rename(state,target);
      // Publication has succeeded; leftover empty staging is harmless.
      fs::remove_empty_directory(staging);
    }
  });
  checkpoint.directory=target;
  return checkpoint;
}
} // namespace estimator_checkpoint
