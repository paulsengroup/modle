#include "modle/juicer_contacts.hpp"

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cinttypes>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "fmt/format.h"
#include "modle/contacts.hpp"
#include "modle/utils.hpp"

namespace modle::juicer_contacts {

template <typename I, typename R>
Contact::Contact(I b1, I b2, R contacts) {
  static_assert(std::is_integral<I>());
  static_assert(std::is_arithmetic<R>());
  if (b1 < 0) {
    throw std::runtime_error(
        fmt::format("Exception thrown in Contact constructor: b1 should be a positive integral "
                    "number, is '{}'",
                    b1));
  }
  if (b2 < 0) {
    throw std::runtime_error(
        fmt::format("Exception thrown in Contact constructor: b2 should be a positive integral "
                    "number, is '{}'",
                    b2));
  }
  if (contacts < 0) {
    throw std::runtime_error(fmt::format(
        "Exception thrown in Contact constructor: b2 should be a positive number, is '{}'",
        contacts));
  }
  this->bin1 = static_cast<uint64_t>(b1);
  this->bin2 = static_cast<uint64_t>(b2);
  this->contacts = static_cast<double>(contacts);
}

ContactsParser::ContactsParser(std::string contact_file) : _path(std::move(contact_file)) {
  if (!std::filesystem::exists(this->_path)) {
    throw std::runtime_error(fmt::format("File '{}' does not exist", this->_path));
  }
  if (std::filesystem::is_directory(this->_path)) {
    throw std::runtime_error(fmt::format("'{}' is a directory. Expected a file", this->_path));
  }
}

bool ContactsParser::parse_next(std::string& buff, Contact& c, std::string_view sep) {
  assert(this->_f);
  if (std::getline(this->_f, buff); buff.empty()) {
    if (this->_f.eof()) {
      return false;
    }
    throw std::logic_error(
        "If you see this error, please report it to the developers on GitHub: "
        "ContactsParser::parse_next() "
        "reached a branch that should not be reachable!");
  }
  if (!this->_f && !this->_f.eof()) {
    throw fmt::system_error(errno, "IO error while reading file {}", this->_path);
  }
  std::vector<std::string_view> tmp = absl::StrSplit(buff, sep);
  if (tmp.size() != 3) {
    throw std::runtime_error(
        fmt::format("Malformed record in file '{}': expected 3 fields, got {} in line '{}'",
                    this->_path, tmp.size(), buff));
  }
  utils::parse_numeric_or_throw(tmp[0], c.bin1);
  utils::parse_numeric_or_throw(tmp[1], c.bin2);
  utils::parse_real_or_throw(tmp[2], c.contacts);
  return false;
}

ContactsParser::MatrixProperties ContactsParser::get_matrix_properties(std::string_view sep) {
  Contact c{};
  ContactsParser::MatrixProperties prop{};
  std::string buff;
  while (this->parse_next(buff, c, sep)) {
    prop.min_bin1 = std::min(prop.min_bin1, c.bin1);
    prop.min_bin2 = std::min(prop.min_bin2, c.bin2);
    prop.max_bin1 = std::max(prop.max_bin1, c.bin1);
    prop.max_bin2 = std::max(prop.max_bin2, c.bin2);
    if (prop.bin_size == 0) {
      prop.bin_size = c.bin1 > c.bin2 ? c.bin1 - c.bin2 : c.bin2 - c.bin1;
      continue;
    }
    if (uint64_t new_bin_size = c.bin1 > c.bin2 ? c.bin1 - c.bin2 : c.bin2 - c.bin1;
        new_bin_size != prop.bin_size) {
      throw std::runtime_error(
          fmt::format("Expected bin of size {}, got {}", prop.bin_size, new_bin_size));
    }
  }
  if (!this->_f && !this->_f.eof()) {
    throw fmt::system_error(errno, "IO error while reading file {}", this->_path);
  }
  this->_f.clear();
  this->_f.seekg(0);
  return prop;
}

ContactMatrix<uint32_t> ContactsParser::parse_into_contact_matrix(uint64_t width,
                                                                  std::string_view sep) {
  std::string buff;
  Contact c{};
  const auto matrix_prop = this->get_matrix_properties(sep);
  ContactMatrix<uint32_t> m((matrix_prop.max_bin1 - matrix_prop.min_bin1) / matrix_prop.bin_size,
                            width / matrix_prop.bin_size);
  this->_f = std::ifstream(this->_path);
  if (!this->_f) {
    throw std::runtime_error(fmt::format("Unable to open file '%s' for reading.", this->_path));
  }
  try {
    while (this->parse_next(buff, c, sep)) {
      uint64_t row = (c.bin1 - matrix_prop.min_bin1) / matrix_prop.bin_size;
      uint64_t col = (c.bin2 - matrix_prop.min_bin2) / matrix_prop.bin_size;
      assert(c.contacts >= 0);
      m.set(row, col, static_cast<uint32_t>(std::lround(c.contacts)));
    }
    if (!this->_f.eof()) {
      throw std::runtime_error("An IO error occurred");
    }
  } catch (const std::exception& err) {
    throw std::runtime_error(
        fmt::format("An error occurred while reading file '%s': %s.", this->_path, err.what()));
  }
  return m;
}

}  // namespace modle::juicer_contacts