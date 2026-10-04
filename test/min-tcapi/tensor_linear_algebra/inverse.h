#pragma once
#include "../common/check.h"
template<class E> void test_inverse() {
    using T=gqten::tensor<E>;using R=tcapi::real_t<T>;
    tcapi::gqten_handle ctx;tcapi::create_context(ctx);
    const R tol=R(1000)*std::numeric_limits<R>::epsilon();
    for(const auto dims:{tcapi::shape_t<T>{1,1},{3,3},{2,3,3,2},{2,2,4}}) {
        const int rows=dims.size()==2?1:2;
        std::size_t n=1;for(int i=0;i<rows;++i)n*=dims[i];
        auto a=tcapi::detail::construct<T>(dims,[&](std::size_t index) {
            const auto i=index%n,j=index/n;
            return i==j ? value<E>(10+i,1) : value<E>(R(int(2*i+j)%3)/R(5),R(int(i)-int(j))/R(7));
        });
        tcapi::scale(ctx,a,value<E>(-2,1));
        const auto saved=tcapi::copy(ctx,a);T out;
        tcapi::inverse(ctx,a,rows,out);
        check(out.Shape()==dims,"inverse preserves original shape");
        check(tcapi::close(ctx,a,saved,R(0)),"inverse preserves input");
        for(std::size_t j=0;j<n;++j)for(std::size_t i=0;i<n;++i) {
            E left{},right{};
            for(std::size_t k=0;k<n;++k) {
                left+=a.GetElem(i+k*n)*out.GetElem(k+j*n);
                right+=out.GetElem(i+k*n)*a.GetElem(k+j*n);
            }
            check(std::abs(left-E(i==j))<tol,"A inverse(A)=I");
            check(std::abs(right-E(i==j))<tol,"inverse(A) A=I");
        }
        auto alias=tcapi::copy(ctx,a);tcapi::inverse(ctx,alias,rows);
        check(tcapi::close(ctx,alias,out,tol),"in-place inverse");
        alias=tcapi::copy(ctx,a);tcapi::inverse(ctx,alias,rows,alias);
        check(tcapi::close(ctx,alias,out,tol),"out form inverse aliases input");
    }
    auto a=tcapi::zeros<T>(ctx,{2,2});
    tcapi::set_elem(ctx,a,{0,1},E(2));tcapi::set_elem(ctx,a,{1,0},E(4));
    T out;tcapi::inverse(ctx,a,1,out);
    check(out.GetElem(0)==E(0) && out.GetElem(3)==E(0) && out.GetElem(1)==E(0.5) && out.GetElem(2)==E(0.25),"LU row pivoting");
    const auto saved=tcapi::copy(ctx,out);
    for(int rows:{-1,0,2,3})throws<std::invalid_argument>([&]{tcapi::inverse(ctx,a,rows,out);});
    auto rect=tcapi::zeros<T>(ctx,{2,3});throws<std::invalid_argument>([&]{tcapi::inverse(ctx,rect,1,out);});
    auto scalar=tcapi::fill<T>(ctx,{},E(1));throws<std::invalid_argument>([&]{tcapi::inverse(ctx,scalar,1,out);});
    for(R x:{std::numeric_limits<R>::infinity(),std::numeric_limits<R>::quiet_NaN()}) {
        auto bad=tcapi::fill<T>(ctx,{2,2},E(x));throws<std::domain_error>([&]{tcapi::inverse(ctx,bad,1,out);});
    }
    check(tcapi::close(ctx,out,saved,R(0)),"invalid inverse preserves output");
    for(int mode=0;mode<3;++mode) {
        auto singular=tcapi::fill<T>(ctx,{2,2},mode==0?E(0):E(1));
        if(mode==2)tcapi::scale(ctx,singular,E(0));
        const auto original=tcapi::copy(ctx,singular);
        throws<std::domain_error>([&]{tcapi::inverse(ctx,singular,1,out);});
        check(tcapi::close(ctx,out,saved,R(0)),"singular inverse preserves separate output");
        throws<std::domain_error>([&]{tcapi::inverse(ctx,singular,1);});
        check(tcapi::close(ctx,singular,original,R(0)),"singular in-place inverse preserves input");
    }
    // Tiny nonzero pivots are not classified as singular by a heuristic cutoff.
    auto tiny=tcapi::eye<T>(ctx,2);tcapi::set_elem(ctx,tiny,{1,1},E(std::numeric_limits<R>::epsilon()));
    tcapi::inverse(ctx,tiny,1,out);
    check(std::abs(out.GetElem(3)*E(std::numeric_limits<R>::epsilon())-E(1))<tol,"small nonzero pivot");
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::inverse(ctx,a,1,out);});
    throws<std::logic_error>([&]{tcapi::inverse(ctx,a,1);});
}
