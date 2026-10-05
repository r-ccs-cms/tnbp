#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <mpi.h>
#include "tnbp/tnbp.h"
#include "qasm/any.h"

using Complex = std::complex<double>;
#ifdef TNBP_TEST_CUDA
using Tensor = tcapi::cuda::Tensor<Complex>;
#else
using Tensor = gqten::tensor<Complex>;
#endif
using Context = tcapi::context_handle_t<Tensor>;
using Shape = tcapi::shape_t<Tensor>;
using Dim = tcapi::bond_dim_t<Tensor>;
int checks = 0;
void check(bool ok, const std::string& label) {
    if (!ok) throw std::runtime_error(label);
    ++checks;
}
void near(Complex a, Complex b, const std::string& label, double tol=1e-9) {
    check(std::isfinite(a.real()) && std::isfinite(a.imag()) && std::abs(a-b)<tol,label);
}
Tensor tensor(Context& ctx, const Shape& shape, const std::vector<Complex>& data) {
    return tcapi::assign_from_range<Tensor>(ctx,shape,data.begin(),[shape](const auto& c) {
        return tnbp::address_from_coor(shape,c);
    });
}
std::vector<Complex> host(Context& ctx, const Tensor& t) {
    std::vector<Complex> data(tcapi::size(ctx,t));
    const auto shape=tcapi::shape(ctx,t);
    tcapi::to_range(ctx,t,data.begin(),[shape](const auto& c) {return tnbp::address_from_coor(shape,c);});
    return data;
}
void gpu_assignment(MPI_Comm comm) {
#ifdef TNBP_TEST_CUDA
    int rank, size, count=0;
    MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
    auto cuda_check=[](cudaError_t e) {if(e!=cudaSuccess) throw std::runtime_error(cudaGetErrorString(e));};
    cuda_check(cudaGetDeviceCount(&count));
    // Require launcher isolation, rather than silently putting every rank on GPU 0.
    check(count==1,"Slurm must expose exactly one GPU per rank");
    cuda_check(cudaSetDevice(0));
    cudaDeviceProp prop{}; cuda_check(cudaGetDeviceProperties(&prop,0));
    char hostname[MPI_MAX_PROCESSOR_NAME]{}; int len=0;
    MPI_Get_processor_name(hostname,&len);
    std::ostringstream uuid;
    for(unsigned char b:prop.uuid.bytes) uuid<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(b);
    std::array<char,33> mine{};
    const auto id=uuid.str(); std::copy(id.begin(),id.end(),mine.begin());
    std::vector<char> ids(size*33), hosts(size*MPI_MAX_PROCESSOR_NAME);
    MPI_Allgather(mine.data(),33,MPI_CHAR,ids.data(),33,MPI_CHAR,comm);
    MPI_Allgather(hostname,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,hosts.data(),MPI_MAX_PROCESSOR_NAME,MPI_CHAR,comm);
    for(int i=0;i<size;++i) for(int j=0;j<i;++j) {
        check(std::string(hosts.data()+i*MPI_MAX_PROCESSOR_NAME)==std::string(hosts.data()+j*MPI_MAX_PROCESSOR_NAME),"same-node GPU test");
        check(std::string(ids.data()+i*33)!=std::string(ids.data()+j*33),"distinct physical GPUs per rank");
    }
    const char* visible=std::getenv("CUDA_VISIBLE_DEVICES");
    std::cout<<"GPU rank="<<rank<<" host="<<hostname<<" visible="<<(visible?visible:"<unset>")<<" uuid="<<id<<'\n';
#else
    (void)comm;
#endif
}
void transfers(Context& ctx, MPI_Comm comm) {
    int rank,size; MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&size);
    for(const Shape shape: {Shape{2,3},Shape{5},Shape{}}) {
        std::size_t count=1; for(auto d:shape) count*=d;
        std::vector<Complex> expected(count);
        for(std::size_t i=0;i<count;++i) expected[i]=Complex(i+1,-double(i)-0.5);
        Tensor a;
        if(rank==0) a=tensor(ctx,shape,expected);
        tnbp::MpiBcast(ctx,a,0,comm);
        check(tcapi::shape(ctx,a)==shape,"broadcast shape");
        check(host(ctx,a)==expected,"broadcast values");
        if(size==2) {
            if(rank==0) {tnbp::MpiSend(ctx,a,1,comm); tnbp::MpiRecv(ctx,a,1,comm);}
            else {tnbp::MpiRecv(ctx,a,0,comm); tcapi::cplx_conj(ctx,a); tnbp::MpiSend(ctx,a,0,comm);}
            for(auto& x:expected)x=std::conj(x);
            check(host(ctx,a)==expected,"send/receive followed by backend operation");
        }
    }
}
// Independently contract a two-site MPS on the host. Site order is explicit,
// and each physical amplitude is compared, rather than gauge-dependent factors.
std::vector<Complex> amplitudes(Context& ctx, const std::vector<Tensor>& v,
    const std::vector<int>& sites, const std::map<int,int>& owners, MPI_Comm comm) {
    std::vector<std::vector<Complex>> h(2); Dim bond=0;
    int rank; MPI_Comm_rank(comm,&rank);
    for(int s=0;s<2;++s) {
        Tensor t;
        if(rank==owners.at(s)) {
            auto p=std::find(sites.begin(),sites.end(),s)-sites.begin();
            t=tcapi::copy(ctx,v.at(p));
        }
        tnbp::MpiBcast(ctx,t,owners.at(s),comm);
        const auto sh=tcapi::shape(ctx,t);
        check(sh.size()==2 && sh[1]==2,"two-site MPS dimensions");
        if(s==0)bond=sh[0]; else check(bond==sh[0],"matching MPS bond");
        h[s]=host(ctx,t);
    }
    std::vector<Complex> result(4);
    for(int b=0;b<2;++b)for(int a=0;a<2;++a)for(Dim k=0;k<bond;++k)
        result[a+2*b]+=h[0][k+bond*a]*h[1][k+bond*b];
    return result;
}
void state_near(std::vector<Complex> actual,const std::vector<Complex>& expected,const std::string& label) {
    double norm=0; Complex overlap=0;
    for(std::size_t i=0;i<actual.size();++i){norm+=std::norm(actual[i]);overlap+=std::conj(expected[i])*actual[i];}
    check(norm>0 && std::isfinite(norm) && std::abs(overlap)>0,label+" nonzero norm/overlap");
    const Complex phase=overlap/std::abs(overlap);
    for(std::size_t i=0;i<actual.size();++i)near(actual[i]/std::sqrt(norm)/phase,expected[i],label);
}
void qasm_bp(Context& ctx,MPI_Comm comm) {
    const auto p=qasm::parse_any("OPENQASM 2.0; include \"qelib1.inc\"; qreg q[2]; h q[0]; cx q[0],q[1];");
    const auto edges=tnbp::EdgesFromQasm(p);
    check(edges==std::vector<std::pair<int,int>>{{0,1}},"QASM edge extraction");
    auto op=tnbp::QasmToTPO<Tensor>(ctx,p,edges);
    auto grouped=tnbp::QasmToTPO<Tensor>(ctx,p,edges,std::vector<int>{2});
    check(op.size()==1 && grouped.size()==1,"QASM TPO layers");
    std::vector<Tensor> v,e; std::vector<int> sites,edgeids; std::map<int,int> owners;
    tnbp::InitTensorProductState(ctx,edges,{2,2},v,sites,owners,e,edgeids,comm);
    tnbp::AbsorbTPO(ctx,op[0],edges,v,sites,owners,e,edgeids,comm);
    const double h=1/std::sqrt(2.0);
    const std::vector<Complex> expected{h,0,0,h};
    auto before=amplitudes(ctx,v,sites,owners,comm);
    for(std::size_t i=0;i<4;++i)near(before[i],expected[i],"QASM Bell amplitudes");
    // The explicit gate-count overload must produce the same physical operator.
    std::vector<Tensor> v2,e2; std::vector<int> s2,ei2; std::map<int,int> o2;
    tnbp::InitTensorProductState(ctx,edges,{2,2},v2,s2,o2,e2,ei2,comm);
    tnbp::OptTPObySVD(ctx,edges,grouped[0],1e-12);
    std::vector<Tensor> local_op;
    for(int site:s2) local_op.push_back(tcapi::copy(ctx,grouped[0].at(site)));
    tnbp::AttachTPO(ctx,edges,local_op,v2,s2,o2,e2,ei2,comm);
    const auto second=amplitudes(ctx,v2,s2,o2,comm);
    for(std::size_t i=0;i<4;++i)near(second[i],expected[i],"grouped QASM Bell amplitudes");
    std::vector<Tensor> updated;
    tnbp::BeliefPropagation(ctx,edges,v,sites,owners,edges,e,edgeids,comm,updated);
    e=std::move(updated);
    double residual=0;
    tnbp::BeliefPropagationCondition(ctx,edges,v,sites,owners,e,edgeids,comm,residual);
    check(std::isfinite(residual) && residual<1e-9,"BP tree fixed point");
    std::vector<Dim> dims; std::vector<double> errors;
    tnbp::Truncation(ctx,edges,v,sites,owners,e,edgeids,comm,Dim(2),0.0,0.0,dims,errors);
    check(!dims.empty(),"BP truncation produced dimensions");
    for(auto d:dims)check(d==2,"Bell retains Schmidt rank two");
    for(auto err:errors)near(err,0,"Bell no-discard error");
    state_near(amplitudes(ctx,v,sites,owners,comm),expected,"BP truncation preserves Bell state");
}
void known_truncation(Context& ctx,MPI_Comm comm) {
    const std::vector<std::pair<int,int>> edges{{0,1}};
    std::vector<Tensor> v,e; std::vector<int> sites,edgeids; std::map<int,int> owners;
    tnbp::InitTensorProductState(ctx,edges,{2,2},v,sites,owners,e,edgeids,comm);
    // A real, nondegenerate Schmidt spectrum inside complex tensors. Its
    // rank-one discarded weight is 0.2; the surviving state is |00>.
    for(int site:sites) {
        const double first=site==0?2.0:1.0;
        v.at(std::find(sites.begin(),sites.end(),site)-sites.begin())=tensor(ctx,{2,2},{first,0,0,1});
    }
    std::vector<Tensor> updated;
    tnbp::BeliefPropagation(ctx,edges,v,sites,owners,edges,e,edgeids,comm,updated);
    e=std::move(updated);
    std::vector<Dim> dims;std::vector<double> errors;
    tnbp::Truncation(ctx,edges,v,sites,owners,e,edgeids,comm,Dim(1),0.0,0.0,dims,errors);
    check(!dims.empty() && dims.size()==errors.size(),"truncation diagnostics");
    for(auto d:dims)check(d==1,"rank-one truncation");
    for(auto err:errors)near(err,0.2,"discarded Schmidt weight");
    for(const auto& t:e)check(tcapi::shape(ctx,t)==Shape({1,1}),"truncated messenger shape");
    state_near(amplitudes(ctx,v,sites,owners,comm),{1,0,0,0},"rank-one truncated state");
}
void pauli_smoke(Context& ctx) {
    std::vector<pauli::Term<Complex>> terms{{Complex(1),"ZX"}};
    std::vector<Tensor> ops; std::vector<std::vector<int>> qubits;
    tnbp::SparsePauliToTensorOp(ctx,terms,ops,qubits);
    check(qubits==std::vector<std::vector<int>>{{0,1}},"Pauli qubit order");
    check(ops.size()==1 && tcapi::shape(ctx,ops[0])==Shape({2,2,2,2}),"Pauli tensor shape");
    auto data=host(ctx,ops[0]);
    for(int in1=0;in1<2;++in1)for(int in0=0;in0<2;++in0)
    for(int out1=0;out1<2;++out1)for(int out0=0;out0<2;++out0) {
        const double value=(out0==1-in0 && out1==in1)?(in1?-1:1):0;
        near(data[out0+2*out1+4*in0+8*in1],value,"Pauli ZX dense oracle");
    }
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
    try {
        check(size==1 || size==2,"run with one or two ranks");
        gpu_assignment(MPI_COMM_WORLD);
        Context ctx; tcapi::create_context(ctx);
        transfers(ctx,MPI_COMM_WORLD);
        qasm_bp(ctx,MPI_COMM_WORLD);
        known_truncation(ctx,MPI_COMM_WORLD);
        pauli_smoke(ctx);
        tcapi::destroy_context(ctx);
        std::cout<<"PASS rank="<<rank<<" checks="<<checks<<'\n';
    }catch(const std::exception& e){std::cerr<<"FAIL rank="<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
    MPI_Finalize();
}
