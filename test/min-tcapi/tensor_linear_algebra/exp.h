#pragma once
#include "../common/check.h"
template<class E> void test_exp() {
    using T=gqten::tensor<E>;using R=tcapi::real_t<T>;
    tcapi::gqten_handle ctx;tcapi::create_context(ctx);
    const R tol=R(1000)*std::numeric_limits<R>::epsilon();
    auto near=[&](E a,E b){check(std::abs(a-b)<=tol*std::max(R(1),R(std::abs(b))),"matrix exponential analytic value");};
    for(const auto dims:{tcapi::shape_t<T>{1,1},{3,3},{2,3,3,2},{2,2,4}}) {
        const int rows=dims.size()==2?1:2;
        std::size_t n=1;for(int i=0;i<rows;++i)n*=dims[i];
        for(bool zero:{true,false}) {
            auto a=tcapi::detail::construct<T>(dims,[&](std::size_t i){return E(!zero && i%n==i/n?R(i%n)/R(3)-R(0.5):R(0));});
            tcapi::scale(ctx,a,E(-2));
            const auto saved=tcapi::copy(ctx,a);T out;
            tcapi::exp(ctx,a,rows,out);
            check(out.Shape()==dims,"exp preserves original tensor shape");
            for(std::size_t j=0;j<n;++j)for(std::size_t i=0;i<n;++i)
                near(out.GetElem(i+j*n),i==j?E(std::exp(std::real(a.GetElem(i+j*n)))):E(0));
            check(tcapi::close(ctx,a,saved,R(0)),"exp preserves input");
            auto alias=tcapi::copy(ctx,a);tcapi::exp(ctx,alias,rows);
            check(tcapi::close(ctx,alias,out,tol),"in-place exponential");
            alias=tcapi::copy(ctx,a);tcapi::exp(ctx,alias,rows,alias);
            check(tcapi::close(ctx,alias,out,tol),"out-form aliases input");
        }
    }
    // H=[[0,z],[conj(z),0]], H^2=|z|^2 I. exp(H)=cosh(r)I+sinh(r)H/r.
    const E z=value<E>(R(0.3),R(0.4));const R r=std::abs(z);
    auto conj=[](E x){if constexpr(std::is_arithmetic_v<E>)return x;else return std::conj(x);};
    auto h=tcapi::zeros<T>(ctx,{2,2});tcapi::set_elem(ctx,h,{0,1},z);tcapi::set_elem(ctx,h,{1,0},conj(z));
    T out;tcapi::exp(ctx,h,1,out);
    near(out.GetElem(0),E(std::cosh(r)));near(out.GetElem(3),E(std::cosh(r)));
    near(out.GetElem(2),z*(std::sinh(r)/r));near(out.GetElem(1),conj(z)*(std::sinh(r)/r));
    auto negative=tcapi::copy(ctx,h);tcapi::scale(ctx,negative,E(-1));T inverse;
    tcapi::exp(ctx,negative,1,inverse);
    for(int j=0;j<2;++j)for(int i=0;i<2;++i) {
        E sum{};for(int k=0;k<2;++k)sum+=out.GetElem(i+2*k)*inverse.GetElem(k+2*j);
        near(sum,E(i==j));
    }
    if constexpr(!std::is_arithmetic_v<E>) {
        auto rotated=tcapi::detail::construct<T>({2,2},[&](std::size_t i){return h.GetElem(i)*E(0,-1);});
        rotated.SetScale(E(0,1));T result;tcapi::exp(ctx,rotated,1,result);
        check(tcapi::close(ctx,result,out,tol),"complex lazy factor, logical Hermitian input");
    }
    auto saved=tcapi::copy(ctx,out);
    for(int rows:{-1,0,2,3})throws<std::invalid_argument>([&]{tcapi::exp(ctx,h,rows,out);});
    auto rect=tcapi::zeros<T>(ctx,{2,3});throws<std::invalid_argument>([&]{tcapi::exp(ctx,rect,1,out);});
    auto scalar=tcapi::fill<T>(ctx,{},E(1));throws<std::invalid_argument>([&]{tcapi::exp(ctx,scalar,1,out);});
    for(R x:{std::numeric_limits<R>::infinity(),std::numeric_limits<R>::quiet_NaN()}) {
        auto bad=tcapi::fill<T>(ctx,{2,2},E(x));throws<std::domain_error>([&]{tcapi::exp(ctx,bad,1,out);});
    }
    check(tcapi::close(ctx,out,saved,R(0)),"invalid exponential preserves output");
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::exp(ctx,h,1,out);});
    throws<std::logic_error>([&]{tcapi::exp(ctx,h,1);});
}
