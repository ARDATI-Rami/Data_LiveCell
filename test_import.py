#!/usr/bin/env python3
"""Test that the package can be imported."""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

try:
    import wing_disc_analysis as wda
    print("✓ Package imported successfully!")
    print(f"Version: {wda.__version__}")
    print(f"\nAvailable modules:")
    for module in ['io', 'geometry', 'topology', 'roi', 'visualization', 'utils']:
        print(f"  - {module}")

    # Test a simple function
    page_w, page_h = wda.utils.get_page_dimensions()
    print(f"\nPage dimensions: {page_w:.2f} x {page_h:.2f} inches")

    print("\n✓ All imports successful!")

except Exception as e:
    print(f"✗ Import failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

