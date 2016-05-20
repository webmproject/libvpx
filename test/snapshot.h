/*
 *  Copyright (c) 2016 The WebM project authors. All Rights Reserved.
 *
 *  Use of this source code is governed by a BSD-style license
 *  that can be found in the LICENSE file in the root of the source
 *  tree. An additional intellectual property rights grant can be found
 *  in the file PATENTS.  All contributing project authors may
 *  be found in the AUTHORS file in the root of the source tree.
 */

#ifndef TEST_SNAPSHOT_H_
#define TEST_SNAPSHOT_H_

#include <map>

namespace libvpx_test {

/**
 * Allows capturing and retrieving snapshots of arbitrary blobs of memory,
 * blob size is based on compile time type information.
 *
 * Usage:
 * void example() {
 *   Snapshot snapshot;
 *
 *   int foo = 4;
 *
 *   snapshot(foo);
 *
 *   foo = 10;
 *
 *   assert(snapshot.get(foo) == 4);     // Pass
 *   assert(snapshot.get(foo) == foo);   // Fail (4 != 10)
 *
 *   char bar[10][10];
 *   memset(bar, 3, sizeof(bar));
 *
 *   snapshot(bar);
 *
 *   memset(bar, 8, sizeof(bar));
 *
 *   assert(sum(bar) == 800);                 // Pass
 *   assert(sum(snapshot.get(bar)) == 300);   // Pass
 * }
 */
class Snapshot {
 public:
  virtual ~Snapshot() {
    for (snapshot_map_t::iterator it = snapshots_.begin();
         it != snapshots_.end(); it++) {
      delete[] it->second;
    }
  }

  /**
   * Take new snapshot for object
   */
  template<typename E>
  void take(const E &e) {
    const void *const key = reinterpret_cast<const void*>(&e);

    snapshot_map_t::iterator it = snapshots_.find(key);

    if (it != snapshots_.end())
      delete[] it->second;

    char *const buf = new char[sizeof(E)];

    memcpy(buf, &e, sizeof(E));

    snapshots_[key] = buf;
  }

  /**
   * Same as 'take'
   */
  template<typename E>
  void operator()(const E &e) {
    take(e);
  }

  /**
   * Retrieve last snapshot for object
   */
  template<typename E>
  const E& get(const E &e) const {
    const void *const key = reinterpret_cast<const void*>(&e);

    snapshot_map_t::const_iterator it = snapshots_.find(key);

    assert(it != snapshots_.end());

    return *reinterpret_cast<const E*>(it->second);
  }

 private:
  typedef std::map<const void*, const char*> snapshot_map_t;

  snapshot_map_t snapshots_;
};

}   // namespace libvpx_test

#endif  // TEST_SNAPSHOT_H_
