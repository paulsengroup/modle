#pragma once

#include <array>
#include <initializer_list>
#include <utility>

#include "modle/common/common.hpp"

namespace modle::utils {
// Source: https://www.youtube.com/watch?v=INn3xa4pMfg
template <class Key, class Value, usize Size>
class ConstMap {
  // std::array<std::pair<Key, Value>, Size> _buff{};
  std::array<Key, Size> _keys{};
  std::array<Value, Size> _vals{};

 public:
  // NOLINTNEXTLINE(hicpp-explicit-conversions)
  constexpr ConstMap(const std::initializer_list<std::pair<Key, Value>>& args);

  template <class KeyInputIt, class ValueInputIt>
  constexpr ConstMap(KeyInputIt first_key, KeyInputIt last_key, ValueInputIt first_value);

  class iterator;
  using const_iterator = const iterator;
  using key_type = Key;
  using mapped_type = Value;
  using value_type = std::pair<const Key, Value>;
  using first_type = const Key;
  using second_type = Value;

  [[nodiscard]] constexpr bool empty() const noexcept;
  [[nodiscard]] constexpr usize size() const noexcept;

  [[nodiscard]] constexpr const Value& at(const Key& key) const;
  [[nodiscard]] constexpr const Value& operator[](const Key& key) const;
  [[nodiscard]] constexpr auto find(const Key& key) const noexcept -> const_iterator;
  [[nodiscard]] constexpr bool contains(const Key& key) const noexcept;

  [[nodiscard]] constexpr auto begin() const noexcept -> const_iterator;
  [[nodiscard]] constexpr auto end() const noexcept -> const_iterator;

  [[nodiscard]] constexpr const std::array<Key, Size>& keys() const noexcept;
  [[nodiscard]] constexpr const std::array<Value, Size>& values() const noexcept;

  class iterator {
    friend ConstMap;

    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    constexpr iterator(const ConstMap& map, usize i = 0) noexcept;

   public:
    constexpr iterator() = delete;
    using value_type = std::pair<const Key, Value>;
    using difference_type = isize;
    using pointer = value_type*;
    using reference = value_type&;
    using iterator_category = std::input_iterator_tag;

    const Key* first;     // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)
    const Value* second;  // NOLINT(cppcoreguidelines-non-private-member-variables-in-classes)

    [[nodiscard]] constexpr auto operator*() const -> const value_type;
    constexpr bool operator==(const iterator& other) const noexcept;
    constexpr bool operator!=(const iterator& other) const noexcept;

    constexpr auto operator++() noexcept -> iterator&;
    constexpr auto operator++(int) const noexcept -> iterator;
  };

 private:
  constexpr void detect_duplicate_keys() const;
};
}  // namespace modle::utils

#include "../../../const_map_impl.hpp"
