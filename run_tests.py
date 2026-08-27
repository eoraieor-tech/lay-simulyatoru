"""Ehtiyat test runner — pytest quraşdırılmayıbsa istifadə olunur.

    python run_tests.py            bütün testlər
    python run_tests.py -q         yalnız yekun
    IMEX_SKIP_SLOW=1 python run_tests.py    yavaş testləri keç

pytest varsa, ondan istifadə etmək daha yaxşıdır:
    pytest -v
"""

from __future__ import annotations

import importlib
import os
import sys
import time
import traceback

ROOT = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.join(ROOT, "tests")
sys.path.insert(0, ROOT)
sys.path.insert(0, TESTS)

GREEN, RED, DIM, RESET = "\033[92m", "\033[91m", "\033[90m", "\033[0m"


def collect_modules():
    return sorted(name[:-3] for name in os.listdir(TESTS)
                  if name.startswith("test_") and name.endswith(".py"))


def main():
    quiet = "-q" in sys.argv
    passed, failed = 0, []
    started = time.time()

    for module_name in collect_modules():
        module = importlib.import_module(module_name)
        functions = sorted(n for n in dir(module) if n.startswith("test_"))
        if not quiet:
            print(f"\n{DIM}{module_name}{RESET}")
        for function_name in functions:
            test = getattr(module, function_name)
            t0 = time.time()
            try:
                test()
                passed += 1
                if not quiet:
                    print(f"  {GREEN}PASS{RESET}  {function_name}"
                          f"  {DIM}{time.time() - t0:.2f}s{RESET}")
            except Exception:
                failed.append((module_name, function_name, traceback.format_exc()))
                if not quiet:
                    print(f"  {RED}FAIL{RESET}  {function_name}")

    print(f"\n{'=' * 62}")
    for module_name, function_name, tb in failed:
        print(f"\n{RED}FAIL{RESET} {module_name}.{function_name}\n{tb}")
    status = f"{GREEN}HAMISI KEÇDİ{RESET}" if not failed else f"{RED}{len(failed)} XƏTA{RESET}"
    print(f"{status}   keçdi: {passed}   xəta: {len(failed)}   "
          f"vaxt: {time.time() - started:.1f} san")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
