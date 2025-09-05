import pytest

from wam2layers.utils.calendar import fix_deprecated_frequency


@pytest.mark.parametrize(
    "freq,expected",
    [
        # Deprecated aliases (should be fixed)
        ("H", "h"),
        ("5H", "5h"),
        ("BH", "bh"),
        ("2BH", "2bh"),
        ("CBH", "cbh"),
        ("3CBH", "3cbh"),
        ("T", "min"),
        ("10T", "10min"),
        ("S", "s"),
        ("15S", "15s"),
        ("L", "ms"),
        ("20L", "20ms"),
        ("U", "us"),
        ("30U", "30us"),
        ("N", "ns"),
        ("40N", "40ns"),
        # Valid aliases (should stay unchanged)
        ("MS", "MS"),  # month start
        ("BMS", "BMS"),  # business month start
        ("QS", "QS"),  # quarter start
        ("AS", "AS"),  # year start
        ("W", "W"),  # week
        ("D", "D"),  # day
        ("min", "min"),  # already valid
        ("ms", "ms"),  # already valid
        ("us", "us"),  # already valid
        ("ns", "ns"),  # already valid
        # Edge cases
        ("5min", "5min"),  # valid new alias
        ("", ""),  # empty string
        ("foobar", "foobar"),  # unknown alias, passthrough
    ],
)
def test_fix_freq(freq, expected):
    assert fix_deprecated_frequency(freq) == expected


@pytest.mark.parametrize(
    "freq,not_expected",
    [
        # Regression guards (make sure these do NOT change incorrectly)
        ("MS", "ms"),  # month start must NOT become millisecond
        ("BMS", "bms"),  # business month start must NOT be lowercased
        ("QS", "qs"),  # quarter start must NOT be lowercased
        ("AS", "as"),  # year start must NOT be lowercased
    ],
)
def test_fix_freq_not_expected(freq, not_expected):
    """Test that fix_freq does not rewrite valid aliases into wrong ones."""
    assert fix_deprecated_frequency(freq) != not_expected
