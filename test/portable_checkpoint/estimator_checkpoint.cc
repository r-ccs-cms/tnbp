#include "../../apps/mcmc_sampling/typedef.h"
#include "../common/device.h"
#include "../../apps/estimator_restart/checkpoint.h"
#include "tnbp/tnbp.h"
#include <chrono>
#include <thread>

// Faults are injected here, not through production app options/environment.
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);
  const auto comm=MPI_COMM_WORLD;
  int rank=0;MPI_Comm_rank(comm,&rank);
  try {
    if(argc!=4) throw std::runtime_error("mode destination circuit");
    initialize_tensor_device(comm);
    ContextHandle ctx;tcapi::create_context(ctx);
    namespace fs=tnbp::file_io;
    const std::string mode=argv[1];const std::string target=argv[2];
    const std::vector<std::pair<int,int>> edges{{0,1},{1,2},{2,3}};
    auto owners=tnbp::ContiguousSiteOwners({0,1,2,3},comm);
    std::vector<Tensor> V,E;std::vector<int> sites,ei;
    tnbp::InitTensorProductState(ctx,edges,{2,2,2,2},V,sites,owners,E,ei,comm);
    std::ostringstream settings;settings<<std::setprecision(17)<<100<<' '<<50<<' '<<1e-8<<' '<<1e-8<<' '<<1e-8;
    estimator_checkpoint::Resume resume{"",tnbp::portable_detail::sha256(tnbp::portable_detail::text(argv[3])),settings.str(),3};
    auto checkpoint=estimator_checkpoint::save(target,resume,comm,[&](const std::string& path) {
      auto state=tnbp::SaveTPSPortable(ctx,path,edges,V,sites,owners,E,ei,comm);
      tnbp::portable_detail::phase(comm,[&] {
        if(rank!=0) return;
        if(mode=="fail-after-state") throw std::runtime_error("injected failure after state save");
        if(mode=="fail-write") fs::create_directory(fs::join(path,"estimator.resume.pending"));
        if(mode=="block-destination") {
          fs::create_directory(target);
          tnbp::portable_detail::write_text(fs::join(target,"sentinel"),"preserve existing destination\n");
        }
        if(mode=="pause") {
          const auto ready=fs::join(fs::parent(path),"ready"),release=fs::join(fs::parent(path),"release");
          tnbp::portable_detail::write_text(ready,"state saved; resume not yet written\n");
          const auto end=std::chrono::steady_clock::now()+std::chrono::seconds(20);
          while(!fs::exists(release)) {
            if(std::chrono::steady_clock::now()>end) throw std::runtime_error("publication test release timeout");
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
          }
          fs::remove(ready);fs::remove(release);
        }
      });
      return state;
    });
    if(checkpoint.directory!=target) throw std::runtime_error("wrong published checkpoint directory");
    if(rank==0) std::cout<<"ESTIMATOR_PUBLISHED "<<checkpoint.state_digest<<std::endl;
  } catch(const std::exception& e) {
    std::cerr<<e.what()<<std::endl;MPI_Abort(comm,1);
  }
  MPI_Finalize();
}
