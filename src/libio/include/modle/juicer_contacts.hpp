#pragma once
#include <cstdint>
#include <fstream>
#include <string>

#include "modle/contacts.hpp"

namespace modle::juicer_contacts {

struct Contact {
  Contact() = default;
  template <typename I, typename R>
  Contact(I b1, I b2, R contacts);
  uint64_t bin1;
  uint64_t bin2;
  double contacts;
};

class ContactsParser {
 public:
  explicit ContactsParser(std::string contact_file);
  ContactMatrix<uint32_t> parse_into_contact_matrix(uint64_t width, std::string_view sep = "\t");

 private:
  std::string _path;
  std::ifstream _f{};

  struct MatrixProperties {
    uint64_t min_bin1{UINT64_MAX};
    uint64_t max_bin1{0};
    uint64_t min_bin2{UINT64_MAX};
    uint64_t max_bin2{0};
    uint64_t bin_size{0};
  };

  [[nodiscard]] MatrixProperties get_matrix_properties(std::string_view sep);
  bool parse_next(std::string& buff, Contact& c, std::string_view sep);
};
}  // namespace modle::juicer_contacts