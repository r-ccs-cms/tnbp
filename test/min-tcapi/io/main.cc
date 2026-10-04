#include "../common/check.h"
#include <sstream>
#include <chrono>
struct temp_dir {
    std::filesystem::path path;
    temp_dir() {
        const auto base=std::filesystem::temp_directory_path();
        const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
        for(int i=0;;++i) {
            path=base/("min-tcapi-io-"+std::to_string(stamp)+"-"+std::to_string(i));
            if(std::filesystem::create_directory(path))break;
        }
    }
    ~temp_dir(){std::error_code ec;std::filesystem::remove_all(path,ec);}
};
template<class E> void test() {
    using T=gqten::tensor<E>;using R=tcapi::real_t<T>;
    tcapi::gqten_handle ctx;tcapi::create_context(ctx);
    temp_dir dir;const auto path=dir.path/"tensor.bin";
    for(const auto dims:{tcapi::shape_t<T>{},{2,3}}) {
        auto a=tcapi::fill<T>(ctx,dims,value<E>(3,-1));tcapi::scale(ctx,a,value<E>(2,1));
        if(!dims.empty())tcapi::set_elem(ctx,a,{1,2},value<E>(7,3));
        std::stringstream stream(std::ios::in|std::ios::out|std::ios::binary);
        stream.write("HEAD",4);tcapi::save(ctx,a,stream);
        auto b=tcapi::fill<T>(ctx,{},value<E>(-5,2));tcapi::save(ctx,b,stream);stream.write("TAIL",4);
        stream.seekg(4);auto loaded=tcapi::load<T>(ctx,stream);
        check(tcapi::close(ctx,a,loaded,R(0)),"stream round trip at current position");
        check(a.GetScale()==loaded.GetScale(),"native lazy scale preserved");
        auto second=tcapi::load<T>(ctx,stream);check(tcapi::close(ctx,b,second,R(0)),"second concatenated tensor");
        char tail[4];stream.read(tail,4);check(std::string(tail,4)=="TAIL","trailing data untouched");
        std::ostringstream native(std::ios::binary);a.StreamWrite(native);
        std::ostringstream fresh(std::ios::binary);tcapi::save(ctx,a,fresh);
        check(native.str()==fresh.str(),"byte-compatible with gqten writer");
        std::istringstream old_input(fresh.str(),std::ios::binary);T legacy;legacy.StreamRead(old_input);
        check(tcapi::close(ctx,a,legacy,R(0)),"gqten reads new save");
        std::istringstream new_input(native.str(),std::ios::binary);
        check(tcapi::close(ctx,a,tcapi::load<T>(ctx,new_input),R(0)),"new load reads gqten save");
        tcapi::save(ctx,a,path);check(tcapi::close(ctx,a,tcapi::load<T>(ctx,path),R(0)),"filesystem path");
        const auto name=path.string();tcapi::save(ctx,a,name);check(tcapi::close(ctx,a,tcapi::load<T>(ctx,name),R(0)),"string path");
        tcapi::save(ctx,a,name.c_str());check(tcapi::close(ctx,a,tcapi::load<T>(ctx,name.c_str()),R(0)),"C string path");
        const auto with_suffix=name+"ignored";const std::string_view view(with_suffix.data(),name.size());
        tcapi::save(ctx,a,view);check(tcapi::close(ctx,a,tcapi::load<T>(ctx,view),R(0)),"non-terminated string_view path");
        {std::ofstream file(path,std::ios::binary);tcapi::save(ctx,a,file);tcapi::save(ctx,b,file);check(file.is_open(),"save leaves stream open");}
        {std::ifstream file(path,std::ios::binary);auto x=tcapi::load<T>(ctx,file);auto y=tcapi::load<T>(ctx,file);check(file.is_open() && tcapi::close(ctx,x,a,R(0)) && tcapi::close(ctx,y,b,R(0)),"ifstream consecutive load");}
        for(std::size_t size=0;size<native.str().size();++size) {
            std::istringstream truncated(native.str().substr(0,size),std::ios::binary);
            throws<std::ios_base::failure>([&]{tcapi::load<T>(ctx,truncated);});
        }
        std::istringstream masked("",std::ios::binary);masked.exceptions(std::ios::failbit|std::ios::badbit);
        throws<std::ios_base::failure>([&]{tcapi::load<T>(ctx,masked);});
        std::ostringstream bad;bad.setstate(std::ios::badbit);
        throws<std::ios_base::failure>([&]{tcapi::save(ctx,a,bad);});
    }
    for(std::int32_t dim:{0,-1}) {
        std::ostringstream encoded(std::ios::binary);std::size_t rank=1;
        encoded.write(reinterpret_cast<const char*>(&rank),sizeof(rank));encoded.write(reinterpret_cast<const char*>(&dim),sizeof(dim));
        std::istringstream in(encoded.str(),std::ios::binary);
        throws<std::invalid_argument>([&]{tcapi::load<T>(ctx,in);});
    }
    {
        std::size_t rank=std::numeric_limits<std::size_t>::max();
        std::string bytes(reinterpret_cast<const char*>(&rank),sizeof(rank));
        std::istringstream in(bytes,std::ios::binary);
        throws<std::overflow_error>([&]{tcapi::load<T>(ctx,in);});
    }
    {
        std::size_t rank=3;std::int32_t dim=std::numeric_limits<std::int32_t>::max();
        std::ostringstream encoded(std::ios::binary);encoded.write(reinterpret_cast<const char*>(&rank),sizeof(rank));
        for(int i=0;i<3;++i)encoded.write(reinterpret_cast<const char*>(&dim),sizeof(dim));
        std::istringstream in(encoded.str(),std::ios::binary);
        throws<std::overflow_error>([&]{tcapi::load<T>(ctx,in);});
    }
    auto a=tcapi::zeros<T>(ctx,{2});
    throws<std::ios_base::failure>([&]{tcapi::load<T>(ctx,dir.path/"missing");});
    throws<std::ios_base::failure>([&]{tcapi::save(ctx,a,dir.path/"missing"/"tensor");});
    tcapi::destroy_context(ctx);
    std::stringstream stream;
    throws<std::logic_error>([&]{tcapi::save(ctx,a,stream);});
    throws<std::logic_error>([&]{tcapi::load<T>(ctx,stream);});
}
RUN_FOUR_TYPES(test)
