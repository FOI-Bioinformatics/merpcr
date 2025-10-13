"""
Test PCR size range handling to match me-PCR behavior.

me-PCR handles PCR size ranges (e.g., "200-220") by:
1. Taking the average as the PCR size: (200 + 220) / 2 = 210
2. Adding margin adjustment: (220 - 210 + 1) = 11
3. Total search range: 210 ± (base_margin + 11)

This is critical for compatibility as it affects the search window
for the second primer.

Reference: me-PCR stsmatch.cpp lines 1007-1031
"""

import pytest
from merpcr.core.engine import MerPCR


class TestPCRRangeHandling:
    """Test PCR size range parsing and margin adjustment."""

    def setup_method(self):
        """Set up test instance."""
        self.engine = MerPCR(margin=50)

    def test_simple_pcr_size(self):
        """Test parsing of simple PCR size values."""
        # Single integer value
        pcr_size, margin_adj = self.engine._parse_pcr_size("240")
        assert pcr_size == 240, f"Expected pcr_size=240, got {pcr_size}"
        assert margin_adj == 0, f"Expected margin_adj=0, got {margin_adj}"

        # Different values
        for value in ["100", "500", "1000", "5000"]:
            pcr_size, margin_adj = self.engine._parse_pcr_size(value)
            assert pcr_size == int(value), f"Expected pcr_size={value}, got {pcr_size}"
            assert margin_adj == 0, f"Expected margin_adj=0 for single value, got {margin_adj}"

    def test_pcr_size_range_average(self):
        """Test that PCR size ranges are averaged correctly."""
        # Range "200-220" should give average of 210
        pcr_size, margin_adj = self.engine._parse_pcr_size("200-220")
        assert pcr_size == 210, f"Expected pcr_size=210 for range 200-220, got {pcr_size}"

        # Range "300-400" should give average of 350
        pcr_size, margin_adj = self.engine._parse_pcr_size("300-400")
        assert pcr_size == 350, f"Expected pcr_size=350 for range 300-400, got {pcr_size}"

        # Odd ranges (result in integer division)
        pcr_size, margin_adj = self.engine._parse_pcr_size("100-199")
        # (100 + 199) // 2 = 149
        assert pcr_size == 149, f"Expected pcr_size=149 for range 100-199, got {pcr_size}"

    def test_pcr_size_range_margin_adjustment(self):
        """
        Test margin adjustment calculation for PCR size ranges.

        From me-PCR stsmatch.cpp:1021:
        margin_to_use += (right - pcr_size + 1)

        For range "200-220":
        - left = 200, right = 220
        - pcr_size = (200 + 220) / 2 = 210
        - margin_adjustment = (220 - 210 + 1) = 11
        """
        # Test "200-220"
        pcr_size, margin_adj = self.engine._parse_pcr_size("200-220")
        assert pcr_size == 210, f"Expected pcr_size=210, got {pcr_size}"
        assert margin_adj == 11, f"Expected margin_adj=11 for range 200-220, got {margin_adj}"

        # Test "300-400"
        pcr_size, margin_adj = self.engine._parse_pcr_size("300-400")
        assert pcr_size == 350, f"Expected pcr_size=350, got {pcr_size}"
        # margin_adj = (400 - 350 + 1) = 51
        assert margin_adj == 51, f"Expected margin_adj=51 for range 300-400, got {margin_adj}"

        # Test "100-120"
        pcr_size, margin_adj = self.engine._parse_pcr_size("100-120")
        assert pcr_size == 110, f"Expected pcr_size=110, got {pcr_size}"
        # margin_adj = (120 - 110 + 1) = 11
        assert margin_adj == 11, f"Expected margin_adj=11 for range 100-120, got {margin_adj}"

    def test_pcr_size_zero_uses_default(self):
        """Test that PCR size of 0 uses default value."""
        engine = MerPCR(default_pcr_size=300)

        pcr_size, margin_adj = engine._parse_pcr_size("0")
        assert pcr_size == 300, f"Expected default pcr_size=300 for '0', got {pcr_size}"
        assert margin_adj == 0, f"Expected margin_adj=0 for zero value, got {margin_adj}"

    def test_pcr_size_dash_uses_default(self):
        """Test that PCR size of '-' uses default value."""
        engine = MerPCR(default_pcr_size=240)

        # Just a dash (me-PCR allows this to mean "size unknown")
        pcr_size, margin_adj = engine._parse_pcr_size("-")
        assert pcr_size == 240, f"Expected default pcr_size=240 for '-', got {pcr_size}"
        assert margin_adj == 0, f"Expected margin_adj=0 for dash, got {margin_adj}"

    def test_invalid_pcr_size_uses_default(self):
        """Test that invalid PCR sizes use default value."""
        engine = MerPCR(default_pcr_size=240)

        # Invalid formats should fall back to default
        invalid_formats = ["abc", "12.5", "", "  ", "100-", "-200", "100--200", "x-y"]

        for invalid in invalid_formats:
            pcr_size, margin_adj = engine._parse_pcr_size(invalid)
            assert (
                pcr_size == 240
            ), f"Expected default pcr_size=240 for invalid '{invalid}', got {pcr_size}"
            assert margin_adj == 0, f"Expected margin_adj=0 for invalid format, got {margin_adj}"

    def test_sts_record_stores_adjusted_margin(self):
        """Test that STS records correctly store the adjusted margin."""
        # Create a test STS file with a range
        import tempfile
        import os

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as f:
            # Write STS with PCR size range
            f.write("STS001\tATCGATCGATCG\tGCTAGCTAGCTA\t200-220\tTest STS\n")
            sts_file = f.name

        try:
            engine = MerPCR(margin=50)
            engine.load_sts_file(sts_file)

            # Check that STS records have the adjusted margin
            # Base margin = 50, adjustment for 200-220 = 11
            # Expected margin = 50 + 11 = 61
            found_sts_with_range = False
            for sts in engine.sts_records:
                if sts.id == "STS001":
                    assert sts.pcr_size == 210, f"Expected pcr_size=210, got {sts.pcr_size}"
                    assert (
                        sts.margin == 61
                    ), f"Expected margin=61 (50+11) for range 200-220, got {sts.margin}"
                    found_sts_with_range = True

            assert found_sts_with_range, "Failed to find STS001 in loaded records"

        finally:
            os.unlink(sts_file)

    def test_different_sts_different_margins(self):
        """Test that different STSs can have different margins."""
        import tempfile
        import os

        with tempfile.NamedTemporaryFile(mode="w", suffix=".sts", delete=False) as f:
            # Write multiple STSs with different PCR sizes
            f.write("STS001\tATCGATCGATCG\tGCTAGCTAGCTA\t200\tNo Range\n")
            f.write("STS002\tATCGATCGATCG\tGCTAGCTAGCTA\t200-220\tSmall Range\n")
            f.write("STS003\tATCGATCGATCG\tGCTAGCTAGCTA\t200-300\tLarge Range\n")
            sts_file = f.name

        try:
            engine = MerPCR(margin=50)
            engine.load_sts_file(sts_file)

            sts_margins = {}
            for sts in engine.sts_records:
                if sts.id in ["STS001", "STS002", "STS003"]:
                    sts_margins[sts.id] = sts.margin

            # STS001: no range, margin = 50 + 0 = 50
            assert sts_margins.get("STS001") == 50, "STS001 should have base margin of 50"

            # STS002: range 200-220, margin = 50 + 11 = 61
            assert sts_margins.get("STS002") == 61, "STS002 should have margin of 61 (50+11)"

            # STS003: range 200-300, margin = 50 + 51 = 101
            # avg = 250, adj = (300 - 250 + 1) = 51
            assert sts_margins.get("STS003") == 101, "STS003 should have margin of 101 (50+51)"

        finally:
            os.unlink(sts_file)

    def test_margin_adjustment_mathematical_correctness(self):
        """
        Verify the mathematical correctness of margin adjustment.

        For a range [low, high]:
        - Average PCR size = (low + high) / 2
        - To search the full range, we need:
          - Lower bound: average - (average - low) = low
          - Upper bound: average + (high - average) = high
        - Margin adjustment should cover (high - average)
        - But me-PCR adds 1, so: margin_adj = (high - average + 1)
        """
        test_cases = [
            # (low, high, expected_avg, expected_margin_adj)
            (200, 220, 210, 11),  # (220 - 210 + 1) = 11
            (100, 200, 150, 51),  # (200 - 150 + 1) = 51
            (300, 400, 350, 51),  # (400 - 350 + 1) = 51
            (150, 151, 150, 2),  # (151 - 150 + 1) = 2 (minimal range)
            (1000, 1100, 1050, 51),  # (1100 - 1050 + 1) = 51
        ]

        for low, high, expected_avg, expected_margin_adj in test_cases:
            pcr_size_str = f"{low}-{high}"
            pcr_size, margin_adj = self.engine._parse_pcr_size(pcr_size_str)

            assert (
                pcr_size == expected_avg
            ), f"Range {pcr_size_str}: Expected avg={expected_avg}, got {pcr_size}"
            assert (
                margin_adj == expected_margin_adj
            ), f"Range {pcr_size_str}: Expected margin_adj={expected_margin_adj}, got {margin_adj}"
