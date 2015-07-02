#ifndef ROOT
#define ROOT

#include "TObject.h"
#include <stdexcept>
#include <string>

namespace root {
  template <typename T>
  bool is_class(const TObject& obj) {
    return obj.IsA() == T::Class();
  }
  template <typename T>
  void require_class(const TObject& obj) {
    if (!is_class<T>(obj)) {
      std::string expected = T::Class()->GetName();
      std::string classname = obj.IsA()->GetName();
      std::string err = "expected class of type '" + expected + "' found"
	" '" + classname + "' instead";
      throw std::runtime_error(err);
    }
  }
  template <typename T>
  T* as_a(TObject* obj) {
    require_class<T>(*obj);
    return static_cast<T*>(obj);
  }
}

#endif
