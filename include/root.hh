#ifndef ROOT_HH
#define ROOT_HH

// ROOT utility functions to load things safely
#include "TClass.h"

#include <stdexcept>
#include <string>

namespace root {
  // Most basic wrapper function to check types.
  // Throws an exception if anything goes wrong.
  template <typename T>
  T* as(TObject* obj);

  // check class type
  template <typename T>
  bool is(const TObject& obj);

  // same as above, check for null pointer too
  template <typename T>
  bool is(const TObject* obj);

  // throw exception if wrong class
  template <typename T>
  void require_class(const TObject& obj);
}

// _______________________________________________________________
// implementation

namespace root {
  template <typename T>
  T* as(TObject* obj) {
    if (obj == 0) {
      throw std::runtime_error("TObject is null pointer");
    }
    require_class<T>(*obj);
    return static_cast<T*>(obj);
  }

  template <typename T>
  void require_class(const TObject& obj) {
    if (!is<T>(obj)) {
      std::string expected = T::Class()->GetName();
      std::string classname = obj.IsA()->GetName();
      std::string err = "expected class of type '" + expected + "' found"
	" '" + classname + "' instead";
      throw std::runtime_error(err);
    }
  }

  template <typename T>
  bool is(const TObject* obj) {
    if (obj == 0) return false;
    return is<T>(*obj);
  }
  template <typename T>
  bool is(const TObject& obj) {
    return obj.IsA() == T::Class();
  }
}

#endif
