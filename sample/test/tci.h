#include <string>
#include <vector>

namespace tci {
  template <typename TenT>
  struct tensor_traits {
    using context_handle_t = std::string;
    using shape_t = std::vector<size_t>;
  };

  template <typename TenT>
  using context_handle_t = typename tensor_traits<TenT>::context_handle_t;

  template <typename TenT>
  using shape_t = typename tensor_traits<TenT>::shape_t;

  template <typename TenT>
  void create_context(context_handle_t<TenT> & ctx);

  template <typename TenT,         
	    typename RandNumGen>   
  void random(                     
      context_handle_t<TenT> &ctx,
      const shape_t<TenT> &shape, 
      RandNumGen &gen,
      TenT & a);

  template <typename TenT,         
	    typename RandNumGen>   
  TenT random(                     
      context_handle_t<TenT> &ctx,
      const shape_t<TenT> &shape, 
      RandNumGen &gen);

}
