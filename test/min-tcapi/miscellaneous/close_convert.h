#pragma once
#include "../common/check.h"
#include <cmath>
template<class E> void test_close() {
    using T=gqten::tensor<E>; using R=tcapi::real_t<T>;
    tcapi::gqten_handle ctx; tcapi::create_context(ctx);
    auto a=tcapi::zeros<T>(ctx,{2}); auto b=tcapi::fill<T>(ctx,{2},E(0.75));
    check(tcapi::close(ctx,a,b,R(1)),"close is not squared-sum comparison");
    check(tcapi::close(ctx,a,b,R(0.75)),"inclusive tolerance boundary");
    check(!tcapi::close(ctx,a,b,std::nextafter(R(0.75),R(0))),"below boundary");
    check(tcapi::close(ctx,a,a,R(-0.0)),"negative zero tolerance");
    auto wrong=tcapi::zeros<T>(ctx,{1,2});
    check(!tcapi::close(ctx,a,wrong,R(10)),"shape mismatch despite equal size");
    auto scalar=tcapi::fill<T>(ctx,{},E(3));
    check(tcapi::close(ctx,scalar,scalar,R(0)),"scalar close");
    b.Scale(E(2)); check(!tcapi::close(ctx,a,b,R(1)),"close respects scale");
    if constexpr (!std::is_arithmetic_v<E>) {
        auto x=tcapi::fill<T>(ctx,{},E(0)); auto y=tcapi::fill<T>(ctx,{},E(3,4));
        check(tcapi::close(ctx,x,y,R(5)),"complex modulus boundary");
        check(!tcapi::close(ctx,x,y,R(4)),"complex modulus not component maximum");
    }
    const R inf=std::numeric_limits<R>::infinity(), nan=std::numeric_limits<R>::quiet_NaN();
    throws<std::invalid_argument>([&]{tcapi::close(ctx,a,b,R(-1));});
    throws<std::invalid_argument>([&]{tcapi::close(ctx,a,wrong,nan);});
    check(tcapi::close(ctx,a,b,inf),"infinite tolerance accepts finite values");
    for(R x:{inf,-inf,nan}) {
        auto nonfinite=tcapi::fill<T>(ctx,{2},E(x));
        check(!tcapi::close(ctx,nonfinite,nonfinite,inf),"nonfinite self comparison is false");
    }
    auto huge=tcapi::fill<T>(ctx,{},E(std::numeric_limits<R>::max()));
    auto minus=tcapi::fill<T>(ctx,{},E(-std::numeric_limits<R>::max()));
    check(!tcapi::close(ctx,huge,minus,std::numeric_limits<R>::max()),"overflowing difference exceeds finite tolerance");
    check(tcapi::close(ctx,huge,minus,inf),"infinite tolerance and difference overflow");
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::close(ctx,a,b,R(0));});
}
template<class S,class D> void conversion_pair() {
    using A=gqten::tensor<S>; using B=gqten::tensor<D>; using R=tcapi::real_t<B>;
    tcapi::gqten_handle src,dst; tcapi::create_context(src); tcapi::create_context(dst);
    for(const auto& dims:{tcapi::shape_t<A>{},{2,3}}) {
        auto a=tcapi::fill<A>(src,dims,value<S>(1.25,-2.5)); a.Scale(value<S>(2,1));
        const auto logical=value<S>(1.25,-2.5)*value<S>(2,1);
        auto out=tcapi::zeros<B>(dst,{4});
        tcapi::convert(src,a,dst,out);
        check(tcapi::shape(dst,out)==dims,"conversion shape");
        D expected;
        if constexpr(std::is_arithmetic_v<D>) expected=static_cast<D>(std::real(logical));
        else expected=D(static_cast<R>(std::real(logical)),static_cast<R>(std::imag(logical)));
        const tcapi::elem_coors_t<B> c(dims.size(),0);
        check(tcapi::get_elem(dst,out,c)==expected,"conversion scaled logical components");
        if(dims.empty()) check(out.GetRaw()==nullptr,"scalar conversion releases dense output");
        tcapi::set_elem(dst,out,c,D(99));
        check(tcapi::get_elem(src,a,c)==logical,"conversion independent output");
        if constexpr(std::is_same_v<S,D>) {
            const auto* previous=a.GetRaw(); tcapi::convert(src,a,dst,a);
            check(tcapi::get_elem(dst,a,c)==logical,"same-object conversion values");
            if(!dims.empty()) check(a.GetRaw()!=previous,"same-object conversion deep copy");
        }
        tcapi::destroy_context(src);
        throws<std::logic_error>([&]{tcapi::convert(src,a,dst,out);});
        check(tcapi::get_elem(dst,out,c)==D(99),"invalid source preserves output");
        tcapi::create_context(src); tcapi::destroy_context(dst);
        throws<std::logic_error>([&]{tcapi::convert(src,a,dst,out);}); tcapi::create_context(dst);
    }
}
template<class S> void test_convert_from() {
    conversion_pair<S,float>(); conversion_pair<S,double>();
    conversion_pair<S,std::complex<float>>(); conversion_pair<S,std::complex<double>>();
}
inline void test_conversion_specials() {
    tcapi::gqten_handle ctx; tcapi::create_context(ctx);
    using C= gqten::tensor<std::complex<double>>;
    using F= gqten::tensor<float>;
    const double inf=std::numeric_limits<double>::infinity();
    auto a=tcapi::fill<C>(ctx,{1},std::complex<double>(-0.0,inf));
    F out;
    tcapi::convert(ctx,a,ctx,out);
    const float* raw=out.GetRaw();
    check(raw[0]==0 && std::signbit(raw[0]),"complex-to-real discards nonfinite imaginary component and keeps signed zero");
    gqten::tensor<std::complex<float>> complex;
    tcapi::convert(ctx,a,ctx,complex);
    check(std::isinf(complex.GetRaw()[0].imag()) && std::signbit(complex.GetRaw()[0].real()),"complex components converted independently");
    auto large=tcapi::fill<gqten::tensor<double>>(ctx,{},std::numeric_limits<double>::max());
    tcapi::convert(ctx,large,ctx,out);
    check(std::isinf(tcapi::get_elem(ctx,out,{})),"narrowing overflow defined as infinity");
    auto nan=tcapi::fill<gqten::tensor<double>>(ctx,{},std::numeric_limits<double>::quiet_NaN());
    tcapi::convert(ctx,nan,ctx,out);
    check(std::isnan(tcapi::get_elem(ctx,out,{})),"NaN conversion");
}
