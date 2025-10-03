// tci/compat.hpp
#pragma once
#include <type_traits>

#if __cplusplus < 202002L
  // remove_cvref_t 代替
  template <class T> struct tci_rm_ref { using type = T; };
  template <class T> struct tci_rm_ref<T&>  { using type = T; };
  template <class T> struct tci_rm_ref<T&&> { using type = T; };
  template <class T> using  tci_remove_reference_t = typename tci_rm_ref<T>::type;

  template <class T> struct tci_rm_cv { using type = T; };
  template <class T> struct tci_rm_cv<const T>          { using type = T; };
  template <class T> struct tci_rm_cv<volatile T>       { using type = T; };
  template <class T> struct tci_rm_cv<const volatile T> { using type = T; };
  template <class T> using  tci_remove_cv_t = typename tci_rm_cv<T>::type;

  template <class T>
  using tci_remove_cvref_t = tci_remove_cv_t<tci_remove_reference_t<T>>;
#else
  #include <utility>
  template <class T>
  using tci_remove_cvref_t = std::remove_cvref_t<T>;
  #include <concepts>
#endif

namespace gqten { template <typename ElemT> class tensor; }

namespace tci {
  // gqten::tensor<...> 判定
  template <typename> struct is_gqten_tensor_impl : std::false_type {};
  template <typename E> struct is_gqten_tensor_impl<gqten::tensor<E>> : std::true_type {};
  template <typename T>
  struct is_gqten_tensor : is_gqten_tensor_impl<tci_remove_cvref_t<T>> {};
  template <typename T>
  constexpr bool is_gqten_tensor_v = is_gqten_tensor<T>::value;

  // 要素型（value_type 想定。ElemT 別名ならここを調整）
  template <typename T, bool = is_gqten_tensor_v<T>>
  struct gqten_elem_type;
  template <typename T>
  struct gqten_elem_type<T, true> { using type = typename tci_remove_cvref_t<T>::value_type; };
  template <typename T>
  using gqten_elem_t = typename gqten_elem_type<T>::type;
}

#if __cplusplus >= 202002L
  #define TCI_REQUIRES(pred)           requires(pred)
  #define TCI_IS_GQTEN_TENSOR(T)       (::tci::is_gqten_tensor_v<T>)
  #define TCI_ELEM_IS_INTEGRAL(T)      (std::integral<::tci::gqten_elem_t<T>>)
  #define TCI_ELEM_IS_FLOATING(T)      (std::floating_point<::tci::gqten_elem_t<T>>)
#else
  #define TCI_REQUIRES(pred)           , typename std::enable_if<(pred), int>::type = 0
  #define TCI_IS_GQTEN_TENSOR(T)       (::tci::is_gqten_tensor<T>::value)
  #define TCI_ELEM_IS_INTEGRAL(T)      (std::is_integral<::tci::gqten_elem_t<T>>::value)
  #define TCI_ELEM_IS_FLOATING(T)      (std::is_floating_point<::tci::gqten_elem_t<T>>::value)
#endif
