// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "modle/compressed_io/compressed_io.hpp"

#include <archive.h>
#include <archive_entry.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <exception>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>

#include "modle/common/common.hpp"

namespace modle::compressed_io {

namespace detail {

void ArchiveReaderDeleter::operator()(archive* arch) const noexcept {
  if (arch && archive_read_free(arch) != ARCHIVE_OK) {
    try {
      SPDLOG_WARN("failed to close compressed file opened for reading!");
    } catch (...) {  // NOLINT
    }
  }
}

void ArchiveWriterDeleter::operator()(archive* arch) const noexcept {
  if (!arch) {
    return;
  }

  const bool ok = archive_write_close(arch) == ARCHIVE_OK && archive_write_free(arch) == ARCHIVE_OK;

  if (!ok) {
    try {
      SPDLOG_WARN("failed to close compressed file opened for writing!");
    } catch (...) {  // NOLINT
    }
  }
}

void ArchiveEntryDeleter::operator()(archive_entry* entry) const noexcept {
  if (entry) {
    archive_entry_free(entry);
  }
}

}  // namespace detail

static void handle_errors(archive* arc, int ec) {
  switch (ec) {
    case ARCHIVE_OK:
      return;
    case ARCHIVE_WARN: {
      SPDLOG_WARN("{}", archive_error_string(arc));
      return;
    }
    case ARCHIVE_FAILED:
      throw std::runtime_error(archive_error_string(arc));
    case ARCHIVE_FATAL: {
      const std::string msg{archive_error_string(arc)};
      archive_write_free(arc);
      throw std::runtime_error(msg);
    }
    default:
      throw std::runtime_error(fmt::format("unknown error code {}", ec));
  }
}

static void handle_errors(archive* arc, la_ssize_t ec) {
  if (ec < 0) {
    assert(ec >= std::numeric_limits<int>::min());
    handle_errors(arc, utils::conditional_static_cast<int>(ec));
  }
}

static void try_init_write_archive(Writer::archive_ptr_t& arc) {
  if (arc) {
    handle_errors(arc.get(), archive_write_close(arc.get()));
    return;
  }

  arc.reset(archive_write_new());
  if (!arc) {
    throw std::runtime_error("unable to allocate memory for compressor");
  }
}

static void try_init_read_archive(Reader::archive_ptr_t& arc) {
  if (arc) {
    handle_errors(arc.get(), archive_read_close(arc.get()));
    return;
  }

  arc.reset(archive_read_new());
  if (!arc) {
    throw std::runtime_error("unable to allocate memory for decompressor");
  }
}

template <typename ArchiveEntryPtr>
static void try_init_archive_entry(ArchiveEntryPtr& entry) {
  static_assert(std::is_same_v<ArchiveEntryPtr, Reader::archive_entry_ptr_t> ||
                std::is_same_v<ArchiveEntryPtr, Writer::archive_entry_ptr_t>);
  if (!entry) {
    entry.reset(archive_entry_new());
    if (entry) {
      return;
    }
    if constexpr (std::is_same_v<ArchiveEntryPtr, Reader::archive_entry_ptr_t>) {
      throw std::runtime_error("unable to allocate memory for decompressor");
    } else {
      throw std::runtime_error("unable to allocate memory for compressor");
    }
  }
  archive_entry_clear(entry.get());
}

Reader::Reader(const std::filesystem::path& path, usize buff_capacity) {
  assert(buff_capacity != 0);
  _buff.reserve(buff_capacity);
  open(path);
}

void Reader::open(const std::filesystem::path& path) {
  if (is_open()) {
    close();
  }

  if (path.empty()) {
    throw std::runtime_error("path is empty");
  }
  _path = path;

  try_init_read_archive(_arc);
  try_init_archive_entry(_arc_entry);

  handle_errors(_arc.get(), archive_read_support_filter_all(_arc.get()));
  handle_errors(_arc.get(), archive_read_support_format_empty(_arc.get()));
  handle_errors(_arc.get(), archive_read_support_format_raw(_arc.get()));
  handle_errors(_arc.get(),
                archive_read_open_filename(_arc.get(), _path.c_str(), _buff.capacity()));

  const auto ec = archive_read_next_header2(_arc.get(), _arc_entry.get());
  _eof = ec == ARCHIVE_EOF;

  _idx = 0;
  if (_eof) {
    return;
  }

  if (ec < 0) {
    handle_errors(_arc.get(), ec);
  }
}

Reader::operator bool() const { return is_open() && !eof(); }

bool Reader::operator!() const { return !operator bool(); }

bool Reader::eof() const noexcept {
  assert(is_open());
  return _eof;
}

bool Reader::is_open() const noexcept { return !!_arc; }

void Reader::close() {
  if (is_open()) {
    _path.clear();
    _arc.reset();
    _arc_entry.reset();
    _buff.clear();
    _idx = 0;
    _eof = false;
  }
}

void Reader::reset() {
  const auto path = _path;
  close();
  _buff.clear();
  _tok_tmp_buff.clear();
  open(path);
}

const std::filesystem::path& Reader::path() const noexcept { return _path; }
std::string Reader::path_string() const noexcept { return _path.string(); }
const char* Reader::path_c_str() const noexcept { return _path.c_str(); }

bool Reader::getline(std::string& buff, char sep) {
  assert(is_open());
  buff.clear();
  if (eof()) {
    return false;
  }

  while (!read_next_token(buff, sep)) {
    if (!read_next_chunk()) {
      assert(eof());
      return !buff.empty();
    }
  }
  return true;
}

std::string_view Reader::getline(char sep) {
  assert(is_open());
  if (eof()) {
    return {};
  }

  _tok_tmp_buff.clear();
  while (true) {
    if (const auto tok = read_next_token(sep); !tok.empty()) {
      return tok;
    }
    if (!read_next_chunk()) {
      assert(eof());
      return {};
    }
  }
}

bool Reader::readall(std::string& buff, char sep) {
  assert(is_open());
  buff.clear();
  if (eof()) {
    return false;
  }

  while (!eof()) {
    std::ignore = getline(sep);
    buff.append(_buff.begin(), _buff.end());
    _idx = 0;
    _buff.clear();
  }
  return true;
}

std::string Reader::readall(char sep) {
  std::string buff;
  readall(buff, sep);
  return buff;
}

bool Reader::read_next_chunk() {
  assert(!eof());
  assert(is_open());
  _buff.resize(_buff.capacity());
  const auto bytes_read = archive_read_data(_arc.get(), _buff.data(), _buff.capacity());
  if (bytes_read < 0) {
    handle_errors(_arc.get(), bytes_read);
  } else if (bytes_read == 0) {
    _eof = true;
    _buff.clear();
    _tok_tmp_buff.clear();
    return false;
  }
  _buff.resize(static_cast<usize>(bytes_read));
  _idx = 0;
  return true;
}

bool Reader::read_next_token(std::string& buff, char sep) {
  assert(!eof());
  assert(is_open());
  assert(_idx <= _buff.size());
  if (_idx == _buff.size()) {
    return false;
  }

  const auto pos = _buff.find(sep, _idx);
  const auto i = static_cast<i64>(_idx);
  if (pos == std::string::npos) {
    buff.append(_buff.begin() + i, _buff.end());
    return false;
  }

  assert(pos >= _idx);
  buff.append(_buff.begin() + i, _buff.begin() + static_cast<i64>(pos));
  _idx = pos + 1;
  return true;
}

std::string_view Reader::read_next_token(char sep) {
  assert(!eof());
  assert(is_open());
  assert(_idx <= _buff.size());
  if (_idx == _buff.size()) {
    return {};
  }

  const auto pos = _buff.find(sep, _idx);
  const auto i = static_cast<i64>(_idx);
  if (pos == std::string::npos) {
    _tok_tmp_buff.append(_buff.begin() + i, _buff.end());
    return {};
  }

  assert(pos >= _idx);
  _idx = pos + 1;
  if (_tok_tmp_buff.empty()) {
    return std::string_view{_buff.data() + static_cast<usize>(i), pos - static_cast<usize>(i)};
  }

  _tok_tmp_buff.append(_buff.begin() + i, _buff.begin() + static_cast<i64>(pos));
  return std::string_view{_tok_tmp_buff};
}

Writer::Writer(const std::filesystem::path& path, Compression compression)
    : _compression(compression) {
  open(path);
}

static void setup_compression(archive* arc, Writer::Compression compression) {
  assert(arc);

  // using enum breaks GCC 10 when used in this context
  using C = Writer::Compression;
  int ec{};
  switch (compression) {
    case C::GZIP: {
      ec = archive_write_add_filter_gzip(arc);
      break;
    }
    case C::BZIP2: {
      ec = archive_write_add_filter_bzip2(arc);
      break;
    }
    case C::LZMA: {
      ec = archive_write_add_filter_lzma(arc);
      break;
    }
    case C::XZ: {
      ec = archive_write_add_filter_xz(arc);
      break;
    }
    case C::ZSTD: {
      ec = archive_write_add_filter_zstd(arc);
      break;
    }
    case C::NONE: {
      ec = archive_write_add_filter_none(arc);
      break;
    }
    default:
      if constexpr (utils::ndebug_not_defined()) {
        throw std::logic_error("unsupported compression algorithm");
      } else {
        MODLE_UNREACHABLE_CODE;
      }
  }
  handle_errors(arc, ec);

  ec = archive_write_set_format_raw(arc);
  handle_errors(arc, ec);
}

void Writer::open(const std::filesystem::path& path) {
  try {
    if (path.empty()) {
      throw std::runtime_error("path is empty!");
    }
    if (!is_open()) {
      close();
    }

    try_init_write_archive(_arc);
    try_init_archive_entry(_arc_entry);

    const auto compression = _compression == AUTO ? infer_compression_from_ext(path) : _compression;
    setup_compression(_arc.get(), compression);

    _path = path;
    auto ec = archive_write_open_filename(_arc.get(), _path.c_str());
    handle_errors(_arc.get(), ec);

    archive_entry_set_pathname(_arc_entry.get(), "x");  // placeholder
    archive_entry_set_filetype(_arc_entry.get(), AE_IFREG);

    ec = archive_write_header(_arc.get(), _arc_entry.get());
    handle_errors(_arc.get(), ec);

  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format("failed to open \"{}\" for writing: {}", path, e.what()));
  }
}

bool Writer::is_open() const noexcept { return !!_arc && !!_arc_entry; }

void Writer::close() {
  _path.clear();
  _arc.reset();
  _arc_entry.reset();
}

Writer::operator bool() const { return is_open(); }

bool Writer::operator!() const { return !operator bool(); }

const std::filesystem::path& Writer::path() const noexcept { return _path; }
std::string Writer::path_string() const noexcept { return _path.string(); }
const char* Writer::path_c_str() const noexcept { return _path.c_str(); }

Writer::Compression Writer::infer_compression_from_ext(const std::filesystem::path& p) {
  // clang-format off
  static constexpr std::array ext_mappings{
    std::make_pair(".gz"sv, GZIP),
    std::make_pair(".bz2"sv, BZIP2),
    std::make_pair(".xz"sv, XZ),
    std::make_pair(".lzma"sv, LZMA),
    std::make_pair(".zst"sv, ZSTD),
    std::make_pair(".zstd"sv, ZSTD)
  };
  // clang-format on

  auto ext = p.extension().string();
  std::ranges::transform(ext, ext.begin(), [](char c) { return std::tolower(c); });
  const auto* const match =
      std::ranges::find_if(ext_mappings, [&](const auto& mapping) { return mapping.first == ext; });
  if (match == ext_mappings.end()) {
    SPDLOG_WARN(
        "unable to infer compression algorithm for file \"{}\": no compression will be applied!",
        p);
    return NONE;
  }
  return match->second;
}

void Writer::write(std::string_view buff) {
  if (buff.empty()) {
    return;
  }

  try {
    assert(is_open());
    const auto ec = archive_write_data(_arc.get(), buff.data(), buff.size());
    handle_errors(_arc.get(), ec);
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format("failed to write {} bytes of data to file \"{}\": {}",
                                         buff.size(), path(), e.what()));
  }
}

}  // namespace modle::compressed_io
