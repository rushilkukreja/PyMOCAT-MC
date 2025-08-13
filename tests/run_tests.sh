#!/bin/bash
# Test runner script for PyMOCAT-MC
# This script runs all tests without requiring pytest

echo "========================================"
echo "PyMOCAT-MC Test Suite"
echo "========================================"
echo ""

# Store the original directory
ORIGINAL_DIR=$(pwd)

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Change to the repository root (parent of tests directory)
cd "$SCRIPT_DIR/.."

# Check Python version
echo "Python version:"
python3 --version
echo ""

# Test 1: Basic functionality tests
echo "1. Running basic functionality tests..."
echo "----------------------------------------"
python3 tests/test_basic_functionality_simple.py
TEST1_RESULT=$?
echo ""

# Test 2: Import test
echo "2. Running import test..."
echo "----------------------------------------"
python3 tests/test_import.py
TEST2_RESULT=$?
echo ""

# Test 3: Minimal simulation test
echo "3. Running minimal simulation test..."
echo "----------------------------------------"
python3 tests/minimal_test.py
TEST3_RESULT=$?
echo ""

# Test 4: Simple run test
echo "4. Running simple run test..."
echo "----------------------------------------"
python3 tests/test_simple_run.py
TEST4_RESULT=$?
echo ""

# Test 5: Quick Start example
echo "5. Running Quick Start example..."
echo "----------------------------------------"
python3 python_implementation/examples/Quick_Start/quick_start.py
TEST5_RESULT=$?
echo ""

# Summary
echo ""
echo "========================================"
echo "TEST SUMMARY"
echo "========================================"

TOTAL_FAILED=0

if [ $TEST1_RESULT -eq 0 ]; then
    echo "+ Basic functionality tests: PASSED"
else
    echo "X Basic functionality tests: FAILED"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

if [ $TEST2_RESULT -eq 0 ]; then
    echo "+ Import test: PASSED"
else
    echo "X Import test: FAILED"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

if [ $TEST3_RESULT -eq 0 ]; then
    echo "+ Minimal simulation test: PASSED"
else
    echo "X Minimal simulation test: FAILED"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

if [ $TEST4_RESULT -eq 0 ]; then
    echo "+ Simple run test: PASSED"
else
    echo "X Simple run test: FAILED"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

if [ $TEST5_RESULT -eq 0 ]; then
    echo "+ Quick Start example: PASSED"
else
    echo "X Quick Start example: FAILED"
    TOTAL_FAILED=$((TOTAL_FAILED + 1))
fi

echo ""
if [ $TOTAL_FAILED -eq 0 ]; then
    echo "All tests passed!"
    EXIT_CODE=0
else
    echo "$TOTAL_FAILED test(s) failed"
    EXIT_CODE=1
fi

# Return to original directory
cd "$ORIGINAL_DIR"

exit $EXIT_CODE