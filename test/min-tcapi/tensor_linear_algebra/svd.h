#pragma once
#include "../common/check.h"
template<class E> void test_svd() {
    using T=gqten::tensor<E>; using R=tcapi::real_t<T>; using RT=tcapi::real_ten_t<T>;
    tcapi::gqten_handle ctx; tcapi::create_context(ctx);
    const R tol=R(1000)*std::numeric_limits<R>::epsilon();
    for (const auto dims:{tcapi::shape_t<T>{2,3,2},{3,2},{2,3},{2,2}}) {
        const int rows=dims.size()==3?2:1;
        auto a=tcapi::zeros<T>(ctx,dims);
        tcapi::for_each_with_coors(ctx,a,[](E& e,const auto& c) {
            int x=0; for(auto i:c) x=3*x+i;
            e=value<E>((x%5)-2,(x%3)-1);
        });
        tcapi::scale(ctx,a,value<E>(2,1));
        T u,v; RT s;
        tcapi::svd(ctx,a,rows,u,s,v);
        std::size_t m=1,n=1; for(int i=0;i<rows;++i)m*=dims[i];
        for(std::size_t i=rows;i<dims.size();++i)n*=dims[i];
        const auto k=std::min(m,n);
        auto us=tcapi::shape_t<T>(dims.begin(),dims.begin()+rows);us.push_back(k);
        auto vs=tcapi::shape_t<T>{static_cast<int>(k)};vs.insert(vs.end(),dims.begin()+rows,dims.end());
        check(u.Shape()==us && v.Shape()==vs && s.Shape()==tcapi::shape_t<T>({int(k),int(k)}),"SVD folded shapes");
        for(std::size_t j=0;j<k;++j) for(std::size_t i=0;i<k;++i) {
            const auto sij=s.GetElem(i+j*k);
            check(i==j ? sij>=R(0) : sij==R(0),"real diagonal nonnegative sigma");
            if(i==j && i>0) check(sij<=s.GetElem((i-1)*(k+1)),"descending singular values");
            E dotu{},dotv{};
            auto conj=[](E x) { if constexpr(std::is_arithmetic_v<E>) return x; else return std::conj(x); };
            for(std::size_t l=0;l<m;++l) dotu+=conj(u.GetElem(l+i*m))*u.GetElem(l+j*m);
            for(std::size_t l=0;l<n;++l) dotv+=v.GetElem(i+l*k)*conj(v.GetElem(j+l*k));
            check(std::abs(dotu-E(i==j))<tol,"U orthonormal columns");
            check(std::abs(dotv-E(i==j))<tol,"Vdag orthonormal rows");
        }
        for(std::size_t j=0;j<n;++j) for(std::size_t i=0;i<m;++i) {
            E sum{};for(std::size_t l=0;l<k;++l)sum+=u.GetElem(i+l*m)*s.GetElem(l+l*k)*v.GetElem(l+j*k);
            check(std::abs(sum-a.GetElem(i+j*m))<tol*R(20),"SVD reconstruction");
        }
        auto alias=tcapi::copy(ctx,a); T av; RT as;
        tcapi::svd(ctx,alias,rows,alias,as,av);
        check(tcapi::close(ctx,as,s,tol),"SVD input alias");
        throws<std::invalid_argument>([&]{tcapi::svd(ctx,a,0,u,s,v);});
        throws<std::invalid_argument>([&]{tcapi::svd(ctx,a,int(dims.size()),u,s,v);});
        throws<std::invalid_argument>([&]{tcapi::svd(ctx,a,rows,u,s,u);});
    }
    auto zero=tcapi::zeros<T>(ctx,{3,2});T u,v;RT s;
    tcapi::svd(ctx,zero,1,u,s,v);
    check(tcapi::norm(ctx,s)==R(0),"zero matrix SVD");
    auto bad=tcapi::fill<T>(ctx,{2,2},E(std::numeric_limits<R>::infinity()));
    throws<std::domain_error>([&]{tcapi::svd(ctx,bad,1,u,s,v);});
    check(tcapi::norm(ctx,s)==R(0),"invalid SVD preserves output");

    if constexpr(std::is_arithmetic_v<E>) {
        auto a=tcapi::eye<T>(ctx,2);
        tcapi::svd(ctx,a,1,u,a,v);
        check(a.Shape()==tcapi::shape_t<T>({2,2}) && std::abs(a.GetElem(0)-R(1))<tol,"sigma aliases real input");
        throws<std::invalid_argument>([&]{tcapi::svd(ctx,a,1,u,u,v);});
    }
    auto scalar=tcapi::fill<T>(ctx,{},E(1));
    throws<std::invalid_argument>([&]{tcapi::svd(ctx,scalar,1,u,s,v);});
    auto spectrum=tcapi::zeros<T>(ctx,{3,3});
    tcapi::set_elem(ctx,spectrum,{0,0},E(4));tcapi::set_elem(ctx,spectrum,{1,1},E(3));
    R err=R(-1);
    tcapi::trunc_svd(ctx,spectrum,1,u,s,v,err,1,R(0));
    check(s.Shape()==tcapi::shape_t<T>({1,1}) && std::abs(s.GetElem(0)-R(4))<tol,"fixed maximum rank");
    check(std::abs(err-R(9)/R(25))<tol,"relative discarded weight");
    for(int j=0;j<3;++j) for(int i=0;i<3;++i) {
        E reconstructed=u.GetElem(i)*s.GetElem(0)*v.GetElem(j);
        check(std::abs(reconstructed-(i==0 && j==0?E(4):E(0)))<tol,"truncated factors reconstructed");
    }
    tcapi::trunc_svd(ctx,spectrum,1,u,s,v,err,3,3,R(0),R(3));
    check(s.Shape()==tcapi::shape_t<T>({2,2}) && err==R(0),"cutoff equality kept, minimum does not restore discarded values");
    tcapi::trunc_svd(ctx,spectrum,1,u,s,v,err,1,3,R(0),R(0));
    check(s.Shape()==tcapi::shape_t<T>({2,2}),"zero singular value unnecessary at zero target");
    auto equal=tcapi::eye<T>(ctx,2);
    tcapi::trunc_svd(ctx,equal,1,u,s,v,err,1,2,R(0.5),R(0));
    check(s.Shape()==tcapi::shape_t<T>({1,1}) && err==R(0.5),"inclusive target error boundary");
    tcapi::trunc_svd(ctx,equal,1,u,s,v,err,1,2,R(0.49),R(0));
    check(s.Shape()==tcapi::shape_t<T>({2,2}) && err==R(0),"increase rank to meet target");
    tcapi::trunc_svd(ctx,zero,1,u,s,v,err,2,R(0));
    check(s.Shape()==tcapi::shape_t<T>({1,1}) && err==R(0),"zero matrix truncation has zero error");
    const auto saved=tcapi::copy(ctx,s);
    err=R(-7);
    throws<std::domain_error>([&]{tcapi::trunc_svd(ctx,spectrum,1,u,s,v,err,3,R(5));});
    check(tcapi::close(ctx,s,saved,R(0)) && err==R(-7),"empty truncation preserves outputs and error");
    throws<std::invalid_argument>([&]{tcapi::trunc_svd(ctx,spectrum,1,u,s,v,err,0,R(0));});
    throws<std::invalid_argument>([&]{tcapi::trunc_svd(ctx,spectrum,1,u,s,v,err,2,1,R(0),R(0));});
    throws<std::invalid_argument>([&]{tcapi::trunc_svd(ctx,spectrum,1,u,s,v,err,1,3,R(-1),R(0));});
    throws<std::invalid_argument>([&]{tcapi::trunc_svd(ctx,spectrum,1,u,s,v,err,3,R(-1));});
    auto alias=tcapi::copy(ctx,spectrum);
    tcapi::trunc_svd(ctx,alias,1,u,s,alias,err,1,R(0));
    check(alias.Shape()==tcapi::shape_t<T>({1,3}),"truncated Vdag aliases input");
    for(R magnitude:{std::numeric_limits<R>::max()/R(16),std::numeric_limits<R>::min()*R(16)}) {
        auto extreme=tcapi::eye<T>(ctx,2);tcapi::scale(ctx,extreme,E(magnitude));
        tcapi::trunc_svd(ctx,extreme,1,u,s,v,err,1,R(0));
        check(std::abs(err-R(0.5))<tol,"truncation weights avoid overflow and underflow");
    }
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::svd(ctx,zero,1,u,s,v);});
}
