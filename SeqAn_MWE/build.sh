#!/usr/bin/env bash
# Build the SeqAn v2 standalone smoke test.
#
# SeqAn is vendored as the full upstream repo (submodule at includes/seqan), so
# the include root is includes/seqan/include (one level into the repo), NOT
# includes/. The phased plan's shorthand "-I includes/" predates the submodule
# layout and does not resolve <seqan/...>; this is the corrected flag.
set -euo pipefail
cd "$(dirname "$0")"
g++ -std=gnu++17 -I ../includes/seqan/include smoke.cpp -o smoke
echo "built ./smoke"
