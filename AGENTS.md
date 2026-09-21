# AGENTS.md

## Structure
- A library for geometric sets with a core in C++
- The python bindings are generated through `nanobind` in `bindings`
- unary operations are implemented in the set classes
- binary operations are implemented as functions
- status of implementations is depicted in `documentation/table.md`

## Setup commands
- Use `uv` for python environment
- Activate with `source .venv/bin/activate`
- Reinstall the python bindings with `uv sync --reinstall`
- Run python tests with `pytest tests`

## C++ Test Setup
- Build with tests enabled: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DGS_BUILD_TESTS=ON`
- Build test executable: `cmake --build build --target geosets_tests`
- Run all tests: `./build/tests/cpp/geosets_tests`
- Run specific tests: `./build/tests/cpp/geosets_tests "[interval]"`
- Run with verbose output: `./build/tests/cpp/geosets_tests --success`

## Code style
- Minimize duplicate code
- Avoid fluff during implementations!
- Use clean Eigen vector syntax: `Eigen::Vector2d(x, y)` instead of `Eigen::VectorXd lb(2); lb << x, y;`

## Test structure
- There are some common tests that apply to all classes in `tests/cpp/common_tests.cpp` and `tests/python/common_tests.py`
- Try to implement common tests whenever possible
- Additional set representation specific tests are only in C++
- Keep specific tests lean and focused - avoid redundancy with common tests
- Add cases with negative values to specific tests, as these are often not covered by the common tests
