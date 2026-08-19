# Dead Code Audit

Source: `cargo build --message-format=json` (rustc's own `dead_code`/`unused_*` lints — not a heuristic tool, unlike the R-to-Python skill's `vulture` wrapper).

## 1. Definitely-dead candidates (no Python analogue, not `pub`)

_None._

## 2. Keep for Python parity (mirrors Python source even if currently unused)

_None._

## 3. Public-API surface flagged despite being `pub`

_A hit here usually means the item is `pub` but its enclosing module is private — genuinely unreachable from outside the crate despite the `pub` keyword. Either widen the module's visibility or remove the now-misleading `pub`._

_None._

## 4. Stale TODO / FIXME / XXX / HACK comments

_None._

## 5. Summary

| Category | Count |
|---|---:|
| Definitely dead (recommend delete) | 0 |
| Keep for Python parity | 0 |
| Public API (flagged anyway) | 0 |
| Stale comments | 0 |
| Python function/class/method symbols (total) | 184 |
| Rust `pub` symbols (total) | 320 |
