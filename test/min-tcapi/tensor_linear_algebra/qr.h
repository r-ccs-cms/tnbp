#pragma once
#include "../common/check.h"
template<class E> void test_qr_lq() {
    using T=gqten::tensor<E>; using R=tcapi::real_t<T>;
    tcapi::gqten_handle ctx; tcapi::create_context(ctx);
    const R tol=R(1000)*std::numeric_limits<R>::epsilon();
    auto conj=[](E x) { if constexpr(std::is_arithmetic_v<E>) return x; else return std::conj(x); };
    for(const auto dims:{tcapi::shape_t<T>{4,2},{2,4},{3,3},{2,3,2,2},{1,3},{3,1}})
    for(int rows=1;rows<int(dims.size());++rows)
    for(bool zero:{false,true}) {
        auto a=tcapi::zeros<T>(ctx,dims);
        if(!zero) tcapi::for_each_with_coors(ctx,a,[](E& e,const auto& c) {
            int x=0;for(auto i:c)x=3*x+i;
            e=value<E>((x%7)-3,(x%5)-2);
        });
        tcapi::scale(ctx,a,value<E>(2,1));
        auto saved=tcapi::copy(ctx,a);
        std::size_t m=1,n=1;
        for(int i=0;i<rows;++i)m*=dims[i];
        for(std::size_t i=rows;i<dims.size();++i)n*=dims[i];
        const auto k=std::min(m,n);
        auto ls=tcapi::shape_t<T>(dims.begin(),dims.begin()+rows);ls.push_back(k);
        auto rs=tcapi::shape_t<T>{int(k)};rs.insert(rs.end(),dims.begin()+rows,dims.end());
        for(bool is_lq:{false,true}) {
            auto decompose=[&](const T& in,T& left,T& right) {
                if(is_lq)tcapi::lq(ctx,in,rows,left,right);
                else tcapi::qr(ctx,in,rows,left,right);
            };
            T left,right;decompose(a,left,right);
            check(left.Shape()==ls && right.Shape()==rs,"QR/LQ folded thin shapes");
            check(tcapi::close(ctx,a,saved,R(0)),"decomposition preserves input");
            for(std::size_t j=0;j<n;++j)for(std::size_t i=0;i<m;++i) {
                E sum{};for(std::size_t l=0;l<k;++l)sum+=left.GetElem(i+l*m)*right.GetElem(l+j*k);
                check(std::abs(sum-a.GetElem(i+j*m))<tol*R(30),"QR/LQ reconstruction");
            }
            for(std::size_t j=0;j<k;++j)for(std::size_t i=0;i<k;++i) {
                E dot{};
                if(is_lq)for(std::size_t l=0;l<n;++l)dot+=right.GetElem(i+l*k)*conj(right.GetElem(j+l*k));
                else for(std::size_t l=0;l<m;++l)dot+=conj(left.GetElem(l+i*m))*left.GetElem(l+j*m);
                check(std::abs(dot-E(i==j))<tol,"Q orthonormal rows or columns");
            }
            if(is_lq) {
                for(std::size_t j=0;j<k;++j)for(std::size_t i=0;i<std::min(m,j);++i)
                    check(std::abs(left.GetElem(i+j*m))<tol,"L lower trapezoidal");
            } else {
                for(std::size_t j=0;j<n;++j)for(std::size_t i=j+1;i<k;++i)
                    check(std::abs(right.GetElem(i+j*k))<tol,"R upper trapezoidal");
            }
            auto alias=tcapi::copy(ctx,a);T other;
            decompose(alias,alias,other);
            check(tcapi::close(ctx,alias,left,tol) && tcapi::close(ctx,other,right,tol),"left output aliases input");
            alias=tcapi::copy(ctx,a);
            decompose(alias,other,alias);
            check(tcapi::close(ctx,other,left,tol) && tcapi::close(ctx,alias,right,tol),"right output aliases input");
            throws<std::invalid_argument>([&]{decompose(a,left,left);});
        }
    }
    auto a=tcapi::eye<T>(ctx,2);auto left=tcapi::fill<T>(ctx,{},E(7));auto right=tcapi::fill<T>(ctx,{},E(9));
    for(bool is_lq:{false,true}) {
        auto call=[&](const T& in,int rows) { if(is_lq)tcapi::lq(ctx,in,rows,left,right);else tcapi::qr(ctx,in,rows,left,right); };
        for(int rows:{-1,0,2,3})throws<std::invalid_argument>([&]{call(a,rows);});
        auto scalar=tcapi::fill<T>(ctx,{},E(1));throws<std::invalid_argument>([&]{call(scalar,1);});
        auto vector=tcapi::fill<T>(ctx,{3},E(1));throws<std::invalid_argument>([&]{call(vector,1);});
        for(R x:{std::numeric_limits<R>::infinity(),std::numeric_limits<R>::quiet_NaN()}) {
            auto bad=tcapi::fill<T>(ctx,{2,2},E(x));throws<std::domain_error>([&]{call(bad,1);});
        }
        check(left.GetScale()==E(7) && right.GetScale()==E(9),"invalid QR/LQ preserves both outputs");
    }
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::qr(ctx,a,1,left,right);});
    throws<std::logic_error>([&]{tcapi::lq(ctx,a,1,left,right);});
}
