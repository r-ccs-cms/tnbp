#include "../common/check.h"
#include <sstream>
#include <thread>
template<class E> void test() {
    using T=gqten::tensor<E>;
    check(tcapi::version<T>()=="1.0","version baseline");
    tcapi::gqten_handle ctx;tcapi::create_context(ctx);
    for(const auto shape:{tcapi::shape_t<T>{},{2,2}}) {
        auto a=tcapi::fill<T>(ctx,shape,value<E>(3,1));tcapi::scale(ctx,a,E(2));
        const auto original_scale=a.GetScale();
        std::ostringstream buffer;auto* previous=std::cout.rdbuf(buffer.rdbuf());
        const auto precision=std::cout.precision();
        try { tcapi::show(ctx,a); } catch(...) { std::cout.rdbuf(previous);throw; }
        std::cout.rdbuf(previous);
        check(buffer.str().find("Rank: "+std::to_string(shape.size()))!=std::string::npos,"show rank");
        check(buffer.str().find("6")!=std::string::npos,"show logical scaled value");
        check(std::cout.precision()==precision,"show preserves precision");
        check(a.GetScale()==original_scale,"show preserves tensor");
    }
    tcapi::destroy_context(ctx);T a;
    throws<std::logic_error>([&]{tcapi::show(ctx,a);});
}
int main(int argc,char**) {try {
    if(argc>1) {
        using T=gqten::tensor<double>;
        tcapi::gqten_handle ctx;tcapi::create_context(ctx);
        auto a=tcapi::zeros<T>(ctx,{2}); // fill must not emit a nested log
        tcapi::reshape(ctx,a,{1,2}); // metadata must describe the original shape
        tcapi::for_each(ctx,a,[&](double& e){e=double(tcapi::size(ctx,a));});
        try {tcapi::reshape(ctx,a,{3});}catch(const std::invalid_argument&){}
        try {tcapi::for_each(ctx,a,[](double&){throw std::runtime_error("callback");});}catch(const std::runtime_error&){}
        tcapi::size(ctx,a); // nested depth restored after exception
        std::thread one([&]{tcapi::size(ctx,a);}),two([&]{tcapi::size(ctx,a);});one.join();two.join();
        auto uninitialized=tcapi::allocate<T>(ctx,{2});tcapi::size(ctx,uninitialized);
        tcapi::version<T>();tcapi::destroy_context(ctx);
        std::cout<<"PROBE OK\n";return 0;
    }
    test<float>();test<double>();test<std::complex<float>>();test<std::complex<double>>();
    std::cout<<"PASS "<<checks<<" checks\n";return 0;
} catch(const std::exception& e) {std::cerr<<"FAIL "<<e.what()<<'\n';return 1;} }
