#!/bin/bash
# Master test runner for Chapter 1: Linear Algebra
# Runs all exercises and reports comprehensive results

set -e  # Exit on error

CHAPTER_DIR="/Users/seancollins/Desktop/Pen&Paper/exercises/Chapter 1"
cd "$CHAPTER_DIR"

echo "╔════════════════════════════════════════════════════════════╗"
echo "║      Chapter 1: Linear Algebra - Test Suite Runner        ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

run_test() {
    local exe=$1
    local name=$2
    
    echo "──────────────────────────────────────────────────────────"
    echo "Testing: $name"
    echo "──────────────────────────────────────────────────────────"
    
    if [ ! -f "./$exe" ]; then
        echo -e "${RED}✗ MISSING${NC}: Executable not found"
        ((FAILED_TESTS++))
        return 1
    fi
    
    # Run test and capture output
    if output=$(./"$exe" 2>&1); then
        # Check for success indicators
        if echo "$output" | grep -q "SUCCESS\|All tests completed"; then
            echo -e "${GREEN}✓ PASSED${NC}"
            ((PASSED_TESTS++))
            
            # Extract test count if available
            if test_count=$(echo "$output" | grep -o "[0-9]\+/[0-9]\+ tests passed" | head -1); then
                echo "  $test_count"
                # Add to total
                passed=$(echo "$test_count" | cut -d'/' -f1)
                total=$(echo "$test_count" | cut -d'/' -f2 | cut -d' ' -f1)
                ((TOTAL_TESTS += total))
            else
                ((TOTAL_TESTS++))
            fi
        else
            echo -e "${RED}✗ FAILED${NC}"
            ((FAILED_TESTS++))
            ((TOTAL_TESTS++))
        fi
    else
        echo -e "${RED}✗ CRASHED${NC}: Exit code $?"
        ((FAILED_TESTS++))
        ((TOTAL_TESTS++))
    fi
    
    echo ""
}

# Run all exercises
run_test "ex_1_1" "1.1 Gram-Schmidt Orthogonalization"
run_test "ex_1_2" "1.2 Linear Transforms"
run_test "ex_1_3" "1.3 Eigenvalue Decomposition"
run_test "ex_1_4" "1.4 Trace, Determinants and Eigenvalues"
run_test "ex_1_5_symmetric" "1.5 Symmetric Matrices"
run_test "ex_1_6_power_method" "1.6 Power Method"

# Final summary
echo "════════════════════════════════════════════════════════════"
echo "                      FINAL SUMMARY                          "
echo "════════════════════════════════════════════════════════════"
echo ""
echo "Total Exercises: 6"
echo "Passed: ${GREEN}$PASSED_TESTS${NC}"
echo "Failed: ${RED}$FAILED_TESTS${NC}"
echo ""

if [ $FAILED_TESTS -eq 0 ]; then
    echo -e "${GREEN}╔════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║  ✓ ALL TESTS PASSED - CHAPTER 1 COMPLETE!         ║${NC}"
    echo -e "${GREEN}╚════════════════════════════════════════════════════╝${NC}"
    exit 0
else
    echo -e "${RED}╔════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║  ✗ SOME TESTS FAILED - REVIEW REQUIRED            ║${NC}"
    echo -e "${RED}╚════════════════════════════════════════════════════╝${NC}"
    exit 1
fi
