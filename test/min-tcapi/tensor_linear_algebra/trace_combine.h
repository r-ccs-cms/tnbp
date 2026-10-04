#pragma once
#include "../common/check.h"
template<class E> void test_trace_combine() {
    using T=gqten::tensor<E>; using R=tcapi::real_t<T>;
    tcapi::gqten_handle ctx; tcapi::create_context(ctx);
    auto a=tcapi::zeros<T>(ctx,{2,3,2,3,2});
    for(int m=0;m<2;++m) for(int l=0;l<3;++l) for(int k=0;k<2;++k) for(int j=0;j<3;++j) for(int i=0;i<2;++i)
        tcapi::set_elem(ctx,a,{i,j,k,l,m},value<E>(1+i+2*j+3*k+4*l+5*m,i-j+k-l+m));
    tcapi::scale(ctx,a,value<E>(2,1));
    T out;
    tcapi::trace(ctx,a,{{3,1},{2,0}},out);
    check(tcapi::shape(ctx,out)==tcapi::shape_t<T>({2}),"partial trace remaining shape");
    for(int m=0;m<2;++m) {
        E expected{};
        for(int j=0;j<3;++j) for(int i=0;i<2;++i) expected+=tcapi::get_elem(ctx,a,{i,j,i,j,m});
        check(tcapi::get_elem(ctx,out,{m})==expected,"multiple reversed trace pairs");
    }
    T single; tcapi::trace(ctx,a,{{0,2}},single);
    check(tcapi::shape(ctx,single)==tcapi::shape_t<T>({3,3,2}),"relative order of remaining axes");
    for(int m=0;m<2;++m) for(int l=0;l<3;++l) for(int j=0;j<3;++j) {
        E expected{}; for(int i=0;i<2;++i) expected+=tcapi::get_elem(ctx,a,{i,j,i,l,m});
        check(tcapi::get_elem(ctx,single,{j,l,m})==expected,"single-pair partial trace");
    }
    auto alias=tcapi::copy(ctx,a); tcapi::trace(ctx,alias,{{0,2},{1,3}},alias);
    check(tcapi::close(ctx,alias,out,R(0)),"trace output alias");
    auto matrix=tcapi::eye<T>(ctx,3); tcapi::scale(ctx,matrix,value<E>(2,1));
    tcapi::trace(ctx,matrix,{{1,0}});
    check(tcapi::shape(ctx,matrix).empty() && matrix.GetRaw()==nullptr,"full trace scalar storage");
    check(tcapi::get_elem(ctx,matrix,{})==value<E>(6,3),"full trace scalar value");
    tcapi::trace(ctx,matrix,{},out); check(tcapi::get_elem(ctx,out,{})==value<E>(6,3),"scalar empty trace");
    tcapi::trace(ctx,a,{},alias); tcapi::set_elem(ctx,alias,{0,0,0,0,0},E(99));
    check(tcapi::get_elem(ctx,a,{0,0,0,0,0})!=E(99),"empty trace independent copy");
    throws<std::out_of_range>([&]{tcapi::trace(ctx,a,{{-1,0}},out);});
    throws<std::out_of_range>([&]{tcapi::trace(ctx,a,{{0,5}},out);});
    throws<std::invalid_argument>([&]{tcapi::trace(ctx,a,{{0,0}},out);});
    throws<std::invalid_argument>([&]{tcapi::trace(ctx,a,{{0,1}},out);});
    throws<std::invalid_argument>([&]{tcapi::trace(ctx,a,{{0,2},{2,4}},out);});
    check(tcapi::get_elem(ctx,out,{})==value<E>(6,3),"trace error preserves output");
    for(const auto& dims:{tcapi::shape_t<T>{},{2,3}}) {
        const auto x=tcapi::fill<T>(ctx,dims,value<E>(2,1));
        auto y=tcapi::fill<T>(ctx,dims,value<E>(1,-2)); tcapi::scale(ctx,y,value<E>(2,1));
        const auto* xp=x.GetRaw(); const auto* yp=y.GetRaw();
        const tcapi::List<tcapi::CRef<T>> ins{std::cref(x),std::cref(y),std::cref(x)};
        auto combined=tcapi::linear_combine(ctx,ins,{E(2),value<E>(-1,1),E(3)});
        const E expected=E(5)*value<E>(2,1)+value<E>(-1,1)*value<E>(1,-2)*value<E>(2,1);
        const tcapi::elem_coors_t<T> c(dims.size(),0);
        check(tcapi::get_elem(ctx,combined,c)==expected,"CRef weighted sum and repeated inputs");
        auto sum=tcapi::linear_combine(ctx,ins);
        check(tcapi::get_elem(ctx,sum,c)==E(2)*value<E>(2,1)+value<E>(1,-2)*value<E>(2,1),"default coefficients all one");
        check(x.GetRaw()==xp && y.GetRaw()==yp,"CRef inputs keep their storage");
        tcapi::set_elem(ctx,combined,c,E(99));
        check(tcapi::get_elem(ctx,x,c)==value<E>(2,1),"combination output independent of const input");
        auto zero=tcapi::linear_combine<T>(ctx,{std::cref(x),std::cref(x)},{E(1),E(-1)});
        check(tcapi::get_elem(ctx,zero,c)==E(0),"cancellation");
    }
    throws<std::invalid_argument>([&]{tcapi::linear_combine<T>(ctx,{});});
    throws<std::invalid_argument>([&]{tcapi::linear_combine<T>(ctx,{std::cref(a)},{});});
    throws<std::invalid_argument>([&]{tcapi::linear_combine<T>(ctx,{std::cref(a),std::cref(matrix)});});
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::trace(ctx,a,{},out);});
    throws<std::logic_error>([&]{tcapi::linear_combine<T>(ctx,{std::cref(a)});});
}
