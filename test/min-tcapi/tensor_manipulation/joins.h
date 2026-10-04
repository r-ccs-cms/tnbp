#pragma once
#include "../common/check.h"
template<class E> void test_joins() {
    using T=gqten::tensor<E>;
    tcapi::context_handle_t<T> ctx; tcapi::create_context(ctx);
    auto a=tcapi::zeros<T>(ctx,{2,3});
    for(int j=0;j<3;++j) for(int i=0;i<2;++i) tcapi::set_elem(ctx,a,{i,j},value<E>(i+2*j+1,2*i-j));
    a.Scale(value<E>(1,2));
    const auto* payload=a.GetRaw();
    for(int axis=0;axis<2;++axis) {
        auto dims=tcapi::shape(ctx,a); dims[axis]=1;
        auto b=tcapi::fill<T>(ctx,dims,value<E>(10,3)); b.Scale(E(2));
        const tcapi::List<tcapi::CRef<T>> refs{std::cref(a),std::cref(b),std::cref(a)};
        auto joined=tcapi::concatenate(ctx,refs,axis);
        auto expected=tcapi::shape(ctx,a); expected[axis]=2*expected[axis]+1;
        check(tcapi::shape(ctx,joined)==expected,"concatenated shape");
        for(int j=0;j<expected[1];++j) for(int i=0;i<expected[0];++i) {
            tcapi::elem_coors_t<T> c{i,j};
            E e;
            if(c[axis]==tcapi::shape(ctx,a)[axis]) e=value<E>(20,6);
            else { if(c[axis]>tcapi::shape(ctx,a)[axis]) c[axis]-=tcapi::shape(ctx,a)[axis]+1; e=tcapi::get_elem(ctx,a,c); }
            check(tcapi::get_elem(ctx,joined,{i,j})==e,"concatenation segment coordinates");
        }
        tcapi::set_elem(ctx,joined,{0,0},E(99));
        check(tcapi::get_elem(ctx,a,{0,0})!=E(99),"concatenation independent storage");
    }
    for(int axis=0;axis<=2;++axis) {
        auto b=tcapi::fill<T>(ctx,{2,3},value<E>(8,-3));
        auto stacked=tcapi::stack<T>(ctx,{std::cref(a),std::cref(b)},axis);
        auto dims=tcapi::shape(ctx,a); dims.insert(dims.begin()+axis,2);
        check(tcapi::shape(ctx,stacked)==dims,"stack shape at each insertion axis");
        for(int n=0;n<2;++n) for(int j=0;j<3;++j) for(int i=0;i<2;++i) {
            tcapi::elem_coors_t<T> c{i,j}; c.insert(c.begin()+axis,n);
            check(tcapi::get_elem(ctx,stacked,c)==(n?value<E>(8,-3):tcapi::get_elem(ctx,a,{i,j})),"stack coordinates");
        }
    }
    check(a.GetRaw()==payload,"join does not replace inputs");
    auto one=tcapi::concatenate<T>(ctx,{std::cref(a)},0);
    tcapi::set_elem(ctx,one,{0,0},E(77));
    check(tcapi::get_elem(ctx,a,{0,0})!=E(77),"single-input concat is deep copy");
    auto scalar=tcapi::fill<T>(ctx,{},value<E>(4,-2));
    auto vector=tcapi::stack<T>(ctx,{std::cref(scalar),std::cref(scalar)},0);
    check(tcapi::shape(ctx,vector)==tcapi::shape_t<T>({2}),"scalar stack vector");
    check(tcapi::get_elem(ctx,vector,{1})==value<E>(4,-2),"scalar stack values");
    auto wrong=tcapi::zeros<T>(ctx,{3,2});
    throws<std::invalid_argument>([&]{tcapi::concatenate<T>(ctx,{},0);});
    throws<std::invalid_argument>([&]{tcapi::stack<T>(ctx,{},0);});
    throws<std::invalid_argument>([&]{tcapi::concatenate<T>(ctx,{std::cref(a),std::cref(wrong)},0);});
    throws<std::invalid_argument>([&]{tcapi::concatenate<T>(ctx,{std::cref(a),std::cref(scalar)},0);});
    throws<std::invalid_argument>([&]{tcapi::stack<T>(ctx,{std::cref(a),std::cref(wrong)},1);});
    throws<std::out_of_range>([&]{tcapi::concatenate<T>(ctx,{std::cref(scalar)},0);});
    throws<std::out_of_range>([&]{tcapi::concatenate<T>(ctx,{std::cref(a)},-1);});
    throws<std::out_of_range>([&]{tcapi::stack<T>(ctx,{std::cref(a)},3);});
    throws<std::out_of_range>([&]{tcapi::stack<T>(ctx,{std::cref(a)},-1);});
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::concatenate<T>(ctx,{std::cref(a)},0);});
    throws<std::logic_error>([&]{tcapi::stack<T>(ctx,{std::cref(a)},0);});
}
