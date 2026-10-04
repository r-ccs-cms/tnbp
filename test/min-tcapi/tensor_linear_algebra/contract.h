#pragma once
#include "../common/check.h"
template<class E> void test_contract() {
    using T=gqten::tensor<E>; using R=tcapi::real_t<T>;
    tcapi::gqten_handle ctx; tcapi::create_context(ctx);
    auto a=tcapi::zeros<T>(ctx,{2,3,4}), b=tcapi::zeros<T>(ctx,{4,3,5});
    for(int i=0;i<2;++i) for(int j=0;j<3;++j) for(int k=0;k<4;++k)
        tcapi::set_elem(ctx,a,{i,j,k},value<E>(i+j+k+1,i-j));
    for(int k=0;k<4;++k) for(int j=0;j<3;++j) for(int l=0;l<5;++l)
        tcapi::set_elem(ctx,b,{k,j,l},value<E>(k+2*j-l,1+j-l));
    tcapi::scale(ctx,a,value<E>(2,1)); tcapi::scale(ctx,b,value<E>(1,-1));
    T out;
    tcapi::contract(ctx,a,"ijk",b,"kjl",out,"li");
    check(tcapi::shape(ctx,out)==tcapi::shape_t<T>({5,2}),"contract output order");
    for(int i=0;i<2;++i) for(int l=0;l<5;++l) {
        E expected{};
        for(int j=0;j<3;++j) for(int k=0;k<4;++k)
            expected+=tcapi::get_elem(ctx,a,{i,j,k})*tcapi::get_elem(ctx,b,{k,j,l});
        check(std::abs(tcapi::get_elem(ctx,out,{l,i})-expected)<R(0.01),"multiple-axis contraction, no conjugation");
    }
    T ints;
    tcapi::contract(ctx,a,{INT32_MIN,-1,INT32_MAX},b,{INT32_MAX,-1,0},ints,{0,INT32_MIN});
    check(tcapi::close(ctx,ints,out,R(0)),"arbitrary signed integer labels");
    auto aa=tcapi::copy(ctx,a), bb=tcapi::copy(ctx,b);
    tcapi::contract(ctx,aa,"ijk",b,"kjl",aa,"li");
    check(tcapi::close(ctx,aa,out,R(0)),"output aliases first input");
    tcapi::contract(ctx,a,"ijk",bb,"kjl",bb,"li");
    check(tcapi::close(ctx,bb,out,R(0)),"output aliases second input");
    E sum{};
    for(int i=0;i<2;++i) for(int j=0;j<3;++j) for(int k=0;k<4;++k) {
        auto x=tcapi::get_elem(ctx,a,{i,j,k}); sum+=x*x;
    }
    aa=tcapi::copy(ctx,a);
    tcapi::contract(ctx,aa,"ijk",aa,"ijk",aa,"");
    check(tcapi::shape(ctx,aa).empty() && aa.GetRaw()==nullptr,"full contraction scalar storage");
    check(std::abs(tcapi::get_elem(ctx,aa,{})-sum)<R(0.01),"both inputs alias output, bilinear sum");
    auto x=tcapi::fill<T>(ctx,{2},value<E>(2,1)), y=tcapi::fill<T>(ctx,{3},value<E>(3,-1));
    tcapi::contract(ctx,x,"i",y,"j",out,"ji");
    check(tcapi::shape(ctx,out)==tcapi::shape_t<T>({3,2}),"outer product shape");
    for(int i=0;i<2;++i) for(int j=0;j<3;++j)
        check(tcapi::get_elem(ctx,out,{j,i})==value<E>(2,1)*value<E>(3,-1),"outer product value");
    auto scalar=tcapi::fill<T>(ctx,{},value<E>(2,-1));
    tcapi::contract(ctx,scalar,"",a,"ijk",out,"kij");
    for(int i=0;i<2;++i) for(int j=0;j<3;++j) for(int k=0;k<4;++k)
        check(tcapi::get_elem(ctx,out,{k,i,j})==value<E>(2,-1)*tcapi::get_elem(ctx,a,{i,j,k}),"scalar with axis reorder");
    aa=tcapi::copy(ctx,a);
    tcapi::contract(ctx,aa,"ijk",scalar,"",aa,"kij");
    check(tcapi::close(ctx,aa,out,R(0)),"right scalar, output alias");
    tcapi::contract(ctx,scalar,"",scalar,"",scalar,"");
    check(tcapi::get_elem(ctx,scalar,{})==value<E>(2,-1)*value<E>(2,-1),"scalar operands alias");
    const auto saved=tcapi::copy(ctx,out);
    throws<std::invalid_argument>([&]{tcapi::contract(ctx,a,"ij",b,"kjl",out,"li");});
    throws<std::invalid_argument>([&]{tcapi::contract(ctx,a,"iij",b,"kjl",out,"li");});
    throws<std::invalid_argument>([&]{tcapi::contract(ctx,a,"ijk",b,"kjl",out,"lli");});
    throws<std::invalid_argument>([&]{tcapi::contract(ctx,a,"ijk",b,"kjl",out,"l");});
    throws<std::invalid_argument>([&]{tcapi::contract(ctx,a,"ijk",b,"kjl",out,"liz");});
    throws<std::invalid_argument>([&]{tcapi::contract(ctx,a,"ijk",b,"kjl",out,"lij");});
    throws<std::invalid_argument>([&]{tcapi::contract(ctx,x,"i",y,"i",out,"");});
    check(tcapi::close(ctx,out,saved,R(0)),"invalid contraction preserves output");
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::contract(ctx,x,"i",y,"j",out,"ij");});
}
