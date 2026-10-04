#pragma once
#include "../common/check.h"
template<class E> void test_components() {
    using T=gqten::tensor<E>;
    tcapi::context_handle_t<T> ctx;
    tcapi::create_context(ctx);
    for(const auto& dims : {tcapi::shape_t<T>{},{2,3}}) {
        auto a=tcapi::fill<T>(ctx,dims,value<E>(2,3));
        a.Scale(value<E>(1,-2)); // complex lazy scaling mixes both components
        const E expected=value<E>(2,3)*value<E>(1,-2);
        const tcapi::elem_coors_t<T> coors(dims.size(),0);
        T out=tcapi::zeros<T>(ctx,{4});
        tcapi::cplx_conj(ctx,a,out);
        E conjugate;
        if constexpr(std::is_arithmetic_v<E>) conjugate=expected;
        else conjugate=std::conj(expected);
        check(tcapi::get_elem(ctx,out,coors)==conjugate,"conjugated logical value");
        check(tcapi::get_elem(ctx,a,coors)==expected,"conjugate leaves input unchanged");
        if(dims.empty()) check(out.GetRaw()==nullptr,"scalar conjugate replaces dense output");
        tcapi::cplx_conj(ctx,out);
        check(tcapi::get_elem(ctx,out,coors)==expected,"double conjugation");
        tcapi::cplx_conj(ctx,out,out);
        check(tcapi::get_elem(ctx,out,coors)==conjugate,"alias conjugate");
        tcapi::set_elem(ctx,out,coors,E(123));
        check(tcapi::get_elem(ctx,a,coors)==expected,"real/complex conjugate deep copy");
        auto complex=tcapi::to_cplx(ctx,a);
        auto real=tcapi::real(ctx,a);
        auto imag=tcapi::imag(ctx,a);
        static_assert(std::is_same_v<decltype(complex),tcapi::cplx_ten_t<T>>);
        static_assert(std::is_same_v<decltype(real),tcapi::real_ten_t<T>>);
        static_assert(std::is_same_v<decltype(imag),tcapi::real_ten_t<T>>);
        check(tcapi::shape(ctx,complex)==dims && tcapi::shape(ctx,real)==dims && tcapi::shape(ctx,imag)==dims,"component shapes");
        check(tcapi::get_elem(ctx,complex,coors)==tcapi::cplx_t<T>(expected),"complex promotion/copy");
        check(tcapi::get_elem(ctx,real,coors)==std::real(expected),"real logical component");
        check(tcapi::get_elem(ctx,imag,coors)==std::imag(expected),"imag logical component");
        tcapi::set_elem(ctx,complex,coors,tcapi::cplx_t<T>(99,88));
        tcapi::set_elem(ctx,real,coors,tcapi::real_t<T>(77));
        check(tcapi::get_elem(ctx,a,coors)==expected,"component outputs independent");
        a.Scale(E(0));
        auto zero_real=tcapi::real(ctx,a);
        auto zero_imag=tcapi::imag(ctx,a);
        check(tcapi::get_elem(ctx,zero_real,coors)==tcapi::real_t<T>(0),"zero scale real");
        check(tcapi::get_elem(ctx,zero_imag,coors)==tcapi::real_t<T>(0),"zero scale imag");
    }
    auto a=tcapi::fill<T>(ctx,{2},E(1));
    a.SetLabels({12});
    auto real=tcapi::real(ctx,a);
    auto complex=tcapi::to_cplx(ctx,a);
    auto imag=tcapi::imag(ctx,a);
    check(real.GetLabels()==a.GetLabels() && complex.GetLabels()==a.GetLabels() && imag.GetLabels()==a.GetLabels(),"component metadata");
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&] { tcapi::cplx_conj(ctx,a); });
    throws<std::logic_error>([&] { tcapi::to_cplx(ctx,a); });
    throws<std::logic_error>([&] { tcapi::real(ctx,a); });
    throws<std::logic_error>([&] { tcapi::imag(ctx,a); });
}
