#pragma once
#include "../common/check.h"
template<class E> void test_eigh() {
    using T=gqten::tensor<E>;using R=tcapi::real_t<T>;using RT=tcapi::real_ten_t<T>;
    tcapi::gqten_handle ctx;tcapi::create_context(ctx);
    const R tol=R(2000)*std::numeric_limits<R>::epsilon();
    auto conj=[](E x){if constexpr(std::is_arithmetic_v<E>)return x;else return std::conj(x);};
    for(const auto dims:{tcapi::shape_t<T>{1,1},{3,3},{2,2,4},{2,3,3,2}}) {
        const int rows=dims.size()==2?1:2;
        std::size_t n=1;for(int i=0;i<rows;++i)n*=dims[i];
        for(int mode=0;mode<4;++mode) {
            auto a=tcapi::detail::construct<T>(dims,[&](std::size_t index) {
                const auto i=index%n,j=index/n;
                if(mode==1)return E(0);
                if(mode==2)return E(i==j?2:0);
                return i==j ? E(int(i)-2) : value<E>(int(i+j)%3-1,int(i)-int(j));
            });
            tcapi::scale(ctx,a,E(-2));
            if constexpr(!std::is_arithmetic_v<E>) if(mode==3) {
                // Non-real stored scale, but the logical matrix stays Hermitian.
                auto rotated=tcapi::detail::construct<T>(dims,[&](std::size_t i){return a.GetRaw()[i]*E(0,-1);});
                tcapi::detail::replace(a,std::move(rotated));
                a.SetScale(E(0,-2));
            }
            auto saved=tcapi::copy(ctx,a);T v;RT w,only;
            tcapi::eigh(ctx,a,rows,w,v);tcapi::eigvalsh(ctx,a,rows,only);
            check(w.Shape()==tcapi::shape_t<T>({int(n),int(n)}) && only.Shape()==tcapi::shape_t<T>({int(n)}),"eigenvalue matrix versus list");
            auto vs=tcapi::shape_t<T>(dims.begin(),dims.begin()+rows);vs.push_back(n);
            check(v.Shape()==vs,"folded eigenvector shape");
            check(tcapi::close(ctx,a,saved,R(0)),"eigensolver preserves input");
            for(std::size_t j=0;j<n;++j) {
                if(j)check(w.GetElem(j*(n+1))>=w.GetElem((j-1)*(n+1)),"ascending eigenvalues including negative scale");
                check(std::abs(w.GetElem(j*(n+1))-only.GetElem(j))<tol*R(20),"eigvalsh agrees with eigh");
                for(std::size_t i=0;i<n;++i) {
                    if(i!=j)check(w.GetElem(i+j*n)==R(0),"lambda off-diagonal zero");
                    E av{},dot{},recon{};
                    for(std::size_t k=0;k<n;++k) {
                        av+=a.GetElem(i+k*n)*v.GetElem(k+j*n);
                        dot+=conj(v.GetElem(k+i*n))*v.GetElem(k+j*n);
                        recon+=v.GetElem(i+k*n)*w.GetElem(k*(n+1))*conj(v.GetElem(j+k*n));
                    }
                    check(std::abs(av-v.GetElem(i+j*n)*w.GetElem(j*(n+1)))<tol*R(30),"A v = lambda v");
                    check(std::abs(dot-E(i==j))<tol,"unitary eigenvectors");
                    check(std::abs(recon-a.GetElem(i+j*n))<tol*R(30),"Hermitian reconstruction");
                }
            }
            auto alias=tcapi::copy(ctx,a);RT aw;
            tcapi::eigh(ctx,alias,rows,aw,alias);
            check(tcapi::close(ctx,aw,w,tol*R(20)),"eigenvector output aliases input");
            if constexpr(std::is_arithmetic_v<E>) {
                alias=tcapi::copy(ctx,a);tcapi::eigvalsh(ctx,alias,rows,alias);
                check(tcapi::close(ctx,alias,only,tol*R(20)),"eigvalsh aliases input");
                alias=tcapi::copy(ctx,a);tcapi::eigh(ctx,alias,rows,alias,v);
                check(tcapi::close(ctx,alias,w,tol*R(20)),"lambda output aliases input");
                throws<std::invalid_argument>([&]{tcapi::eigh(ctx,a,rows,v,v);});
            }
        }
    }
    auto a=tcapi::eye<T>(ctx,2);RT w=tcapi::fill<RT>(ctx,{},R(7));T v=tcapi::fill<T>(ctx,{},E(9));
    for(bool vectors:{false,true}) {
        auto call=[&](const T& in,int rows){if(vectors)tcapi::eigh(ctx,in,rows,w,v);else tcapi::eigvalsh(ctx,in,rows,w);};
        for(int rows:{-1,0,2,3})throws<std::invalid_argument>([&]{call(a,rows);});
        auto rectangular=tcapi::zeros<T>(ctx,{2,3});throws<std::invalid_argument>([&]{call(rectangular,1);});
        auto scalar=tcapi::fill<T>(ctx,{},E(1));throws<std::invalid_argument>([&]{call(scalar,1);});
        for(R x:{std::numeric_limits<R>::infinity(),std::numeric_limits<R>::quiet_NaN()}) {
            auto bad=tcapi::fill<T>(ctx,{2,2},E(x));throws<std::domain_error>([&]{call(bad,1);});
        }
        check(w.GetScale()==R(7) && v.GetScale()==E(9),"validation preserves outputs");
    }
    tcapi::destroy_context(ctx);
    throws<std::logic_error>([&]{tcapi::eigh(ctx,a,1,w,v);});
    throws<std::logic_error>([&]{tcapi::eigvalsh(ctx,a,1,w);});
}
