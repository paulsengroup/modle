#pragma once
#include <cstdint>  // for uint64_t
#include <fstream>
#include <string>
#include <string_view>

#include "modle/contacts.hpp"

namespace modle::juicer_contacts {

struct Contact {
  inline Contact() = default;
  template <typename I, typename R>
  Contact(I b1, I b2, R contacts_num);
  uint64_t bin1;
  uint64_t bin2;
  double contacts;
};

class ContactsParser {
 public:
  inline explicit ContactsParser(std::string contact_file);
  [[nodiscard]] inline ContactMatrix<uint32_t> parse_into_contact_matrix(
      uint64_t width, std::string_view sep = "\t");

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

  [[nodiscard]] inline MatrixProperties get_matrix_properties(std::string_view sep);
  inline bool parse_next(std::string& buff, Contact& c, std::string_view sep);
};

[[nodiscard]] inline ContactMatrix<uint32_t> run_juicer_dump_and_parse_contacts(
    std::string_view path_to_input_matrix, std::string_view path_to_reference_matrix,
    std::string_view chr_name_hic, uint64_t chr_offset_hic, std::string_view path_to_juicer_tools,
    uint64_t juicer_tools_mem);

}  // namespace modle::juicer_contacts

#include "../../juicer_contacts_impl.hpp"