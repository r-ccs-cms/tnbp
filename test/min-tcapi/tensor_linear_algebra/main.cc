#include "../common/check.h"
#include <cmath>
#include "trace_combine.h"
#include "contract.h"
#include "svd.h"
#include "qr.h"
#include "eigh.h"
#include "exp.h"
#include "inverse.h"
template<class E> void test() {
    test_trace_combine<E>();
    test_contract<E>();
    test_svd<E>();
    test_qr_lq<E>();
    test_eigh<E>();
    test_exp<E>();
    test_inverse<E>();
    using T=gqten::tensor<E>; using R=tcapi::real_t<T>;
    tcapi::gqten_handle ctx; tcapi::create_context(ctx);
    const R tolerance=R(64)*std::numeric_limits<R>::epsilon();
    auto near=[&](auto x,auto y) { check(std::abs(x-y)<=tolerance*std::max(R(1),R(std::abs(y))),"numerical tolerance"); };
    auto v=tcapi::zeros<T>(ctx,{3});
    for(int i=0;i<3;++i) tcapi::set_elem(ctx,v,{i},value<E>(i+1,2-i));
    tcapi::scale(ctx,v,value<E>(2,1));
    T d; tcapi::diag(ctx,v,d);
    check(tcapi::shape(ctx,d)==tcapi::shape_t<T>({3,3}),"diag vector to matrix");
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        check(tcapi::get_elem(ctx,d,{i,j})==(i==j?tcapi::get_elem(ctx,v,{i}):E{}),"diag logical values and zeros");
    tcapi::diag(ctx,d,d);
    check(tcapi::close(ctx,v,d,R(0)),"diag alias round trip");
    tcapi::set_elem(ctx,d,{0},E(99));
    check(tcapi::get_elem(ctx,v,{0})!=E(99),"diag independent storage");
    for(const auto& dims:{tcapi::shape_t<T>{2,3},{3,2}}) {
        auto a=tcapi::zeros<T>(ctx,dims);
        for(int j=0;j<dims[1];++j) for(int i=0;i<dims[0];++i) tcapi::set_elem(ctx,a,{i,j},value<E>(i+3*j,1+i-j));
        tcapi::scale(ctx,a,value<E>(1,2));
        auto original=tcapi::copy(ctx,a); tcapi::diag(ctx,a);
        check(tcapi::shape(ctx,a)==tcapi::shape_t<T>({2}),"rectangular diagonal length");
        for(int i=0;i<2;++i) check(tcapi::get_elem(ctx,a,{i})==tcapi::get_elem(ctx,original,{i,i}),"rectangular diagonal stride");
    }
    for(const auto& dims:{tcapi::shape_t<T>{},{2,3}}) {
        auto a=tcapi::fill<T>(ctx,dims,value<E>(3,4));
        tcapi::scale(ctx,a,value<E>(-2,1));
        const E logical=value<E>(3,4)*value<E>(-2,1);
        const auto count=tcapi::size(ctx,a);
        const R expected=std::abs(logical)*std::sqrt(R(count));
        near(tcapi::norm(ctx,a),expected);
        const auto* payload=a.GetRaw();
        auto out=tcapi::zeros<T>(ctx,{4});
        near(tcapi::normalize(ctx,a,out),expected);
        near(tcapi::norm(ctx,out),R(1));
        check(a.GetRaw()==payload,"out normalization leaves input allocation");
        near(tcapi::norm(ctx,a),expected);
        const tcapi::elem_coors_t<T> c(dims.size(),0);
        near(tcapi::get_elem(ctx,out,c),logical/expected);
        near(tcapi::normalize(ctx,a,a),expected);
        near(tcapi::norm(ctx,a),R(1));
        auto saved=tcapi::copy(ctx,a);
        tcapi::scale(ctx,a,E(2),out);
        near(tcapi::norm(ctx,out),R(2));
        check(tcapi::close(ctx,a,saved,tolerance),"out scale leaves input unchanged");
        tcapi::scale(ctx,out,E(-3),out);
        near(tcapi::norm(ctx,out),R(6));
        tcapi::scale(ctx,out,E(0));
        check(tcapi::norm(ctx,out)==R(0),"zero scale norm");
        throws<std::domain_error>([&]{tcapi::normalize(ctx,out);});
        check(tcapi::norm(ctx,out)==R(0),"zero normalization preserves input");
        auto target=tcapi::fill<T>(ctx,{1},E(7));
        throws<std::domain_error>([&]{tcapi::normalize(ctx,out,target);});
        check(tcapi::get_elem(ctx,target,{0})==E(7),"zero normalization preserves output");
    }
    auto scalar=tcapi::fill<T>(ctx,{},E(1));
    throws<std::invalid_argument>([&]{tcapi::diag(ctx,scalar,d);});
    check(tcapi::get_elem(ctx,d,{0})==E(99),"invalid diag preserves output");
    auto rank3=tcapi::zeros<T>(ctx,{1,1,1});
    throws<std::invalid_argument>([&]{tcapi::diag(ctx,rank3);});
    for(R x:{std::numeric_limits<R>::infinity(),std::numeric_limits<R>::quiet_NaN()}) {
        auto bad=tcapi::fill<T>(ctx,{},E(x));
        throws<std::domain_error>([&]{tcapi::normalize(ctx,bad,d);});
        check(tcapi::get_elem(ctx,d,{0})==E(99),"nonfinite norm preserves output");
    }
    auto tiny=tcapi::fill<T>(ctx,{1},E(std::numeric_limits<R>::min()/R(8)));
    const auto n=tcapi::normalize(ctx,tiny);
    check(n>R(0),"tiny positive original norm"); near(tcapi::norm(ctx,tiny),R(1));
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::diag(ctx,v);});
    throws<std::logic_error>([&]{tcapi::norm(ctx,v);});
    throws<std::logic_error>([&]{tcapi::normalize(ctx,v);});
    throws<std::logic_error>([&]{tcapi::scale(ctx,v,E(1));});
}
RUN_FOUR_TYPES(test)
